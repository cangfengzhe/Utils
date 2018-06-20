#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author:	xieshang@grandomics.com
# Date & Time:	2018-01-24 19:15:42
# Description:
# Standard library imports
"""
计算每个read 每个信号值，上游下游的波动情况(采用滑窗，窗口为100，步长为1)，引入了ttest、std、mean以及cv
"""
import argparse
import os
import sys

# Third party imports
try:
    import h5py
    import numpy as np
    from scipy.stats import ttest_ind
except (NameError, ImportError) as E:
    print(E)
    print('A third party package is missing. Please verify your dependencies')
    sys.exit()

# Local lib import


#extract option
def Get_Opt():
    parser = argparse.ArgumentParser(
        prog='detect_ont_str.py',
        usage='%(prog)s [options]',
        description='Tools for analysing ont fast5 data.',
    )
    parser.add_argument(
        '-i',
        '--input',
        type=str,
        required=True,
        dest='input',
        help='input file is a fastq with nanopolish index')
    #parser.add_argument('-a', '--standard-array', type=str, required=False,
    #    dest='standard_array', default='', help='a standard array')
    parser.add_argument(
        '-w',
        '--signal-window',
        type=int,
        required=False,
        dest='window_size',
        default=100,
        help='window size in calculate signal')
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        required=True,
        dest='output',
        help='output file')
    parser.add_argument(
        '-v', '--version', action='version', version='%(prog)s 0.0.1')
    options = parser.parse_args()
    return options


def Read_Index(input):
    Read_Length, Read_Fast5 = dict(), dict()

    index_fai = input + '.index.fai'
    with open(index_fai, 'r') as open_fai:
        for eachline in open_fai.readlines():
            sp = eachline.strip().split()
            read_name = sp[0]
            read_length = int(sp[1])
            Read_Length.update({read_name: read_length})

    index_readdb = input + '.index.readdb'
    with open(index_readdb, 'r') as open_readdb:
        for eachline in open_readdb.readlines():
            sp = eachline.strip().split()
            read_name = sp[0]
            fast5_path = sp[1]
            Read_Fast5.update({read_name: fast5_path})

    return Read_Length, Read_Fast5


def Read_Fast5_File(Read_Fast5):
    """
    Read_Signal return read_name and signal list
    Read_Info return read_name and signal of every position
    """
    Read_Signal, Read_Info = dict(), dict()
    for read_name, fast5_path in Read_Fast5.items():
        Electrical_Signal = list()
        fast5_file = h5py.File(fast5_path, 'r')
        Reads = fast5_file['Raw']['Reads']
        read = Reads.keys()[0]
        Signal = Reads[read]['Signal']
        signal_id = 0
        Read_Info.update({read_name: dict()})
        for signal in Signal:
            signal_id += 1
            Read_Info[read_name].update({signal_id: [str(signal)]})
            Electrical_Signal.append(signal)
            #print read_name, signal
        Read_Signal.update({read_name: Electrical_Signal})
    return Read_Signal, Read_Info


def Ttest_Signal(Read_Length, Read_Signal, Read_Info, window_size):
    """
    左右窗口信号值分别与中间窗口做ttest，步长为1
    """
    for read_name, Signal in Read_Signal.items():
        #read_length = Read_Length[read_name]
        signal_length = len(Signal)

        # left window  vs mid window； right window vs mid window
        for window_id in range(signal_length - window_size * 3):
            signal_array = np.array(
                Signal[window_id + window_size:window_id + window_size * 2])
            left_array = np.array(Signal[window_id:window_id + window_size])
            lt, lp = ttest_ind(left_array, signal_array, equal_var=False)
            right_array = np.array(Signal[window_id + window_size * 2:
                                          window_id + window_size * 3])
            rt, rp = ttest_ind(right_array, signal_array, equal_var=False)
            signal_id = window_id + window_size * 1.5 + 1
            Read_Info[read_name][signal_id].extend([lt, lp, rt, rp])
            #print read_name, signal_id, lt, lp, rt, rp
    return Read_Info


def Read_Standard_Array(standard_array):
    Standard = list()
    with open(standard_array, 'r') as open_standard:
        for eachline in open_standard.readlines():
            st = float(eachline.strip())
            Standard.append(st)
    return Standard


def Calculate_CV(Read_Signal, Read_Info, window_size):
    """
    为了衡量数据间波动情况，我们引入了变异系数coefficient of variation（CV）
    其是概率分布离散程度的一个归一化量度，其定义为标准差与平均值之比
    公式1：
      cv = SD/mean
    """
    Read_CV = dict()
    for read_name, Signal in Read_Signal.items():
        Read_CV.update({read_name: list()})
        #read_length = Read_Length[read_name]
        signal_length = len(Signal)
        #last_CV = 100
        for window_id in range(signal_length - window_size):
            signal_array = np.array(Signal[window_id:window_id + window_size])
            #max_index = np.argmax(signal_array)
            #standard = Standard[100-max_index:]+Standard[:100-max_index]
            #for signal_id in range(100):
            #    signal_array[signal_id] = signal_array[signal_id] / standard[signal_id]
            mean = np.mean(signal_array)
            SD = np.std(signal_array)
            CV = SD / mean * 100   # 标准差/均值 * 100
            Read_CV[read_name].append(CV)
            #pos = float(window_id+window_size/2)/signal_length*read_length
            signal_id = window_id + window_size * 0.5 + 1
            if signal_id > window_size * 1.5:
                Read_Info[read_name][signal_id].extend([mean, SD, CV])
            #print read_name, window_id+51, pos, mean, SD, CV
            #if CV < 10:
            #    if last_CV < 10:
            #        region = [region[0], pos]
            #    else:
            #        region = [pos, pos]
            #else:
            #    if last_CV < 10:
            #        if region[1] - region[0] > 100:
            #            region = [int(round(region[0])), int(round(region[1]))]
            #            Read_Region[read_name].append(region)
            #last_CV = CV
    return Read_CV, Read_Info


def Ttest_CV(Read_Length, Read_CV, Read_Info, window_size):
    for read_name, CV in Read_CV.items():
        #read_length = Read_Length[read_name]
        cv_length = len(CV)
        for window_id in range(cv_length - window_size * 3):
            cv_array = np.array(
                CV[window_id + window_size:window_id + window_size * 2])
            left_array = np.array(CV[window_id:window_id + window_size])
            lt, lp = ttest_ind(left_array, cv_array, equal_var=False)
            right_array = np.array(
                CV[window_id + window_size * 2:window_id + window_size * 3])
            rt, rp = ttest_ind(right_array, cv_array, equal_var=False)
            signal_id = window_id + window_size * 2 + 1
            Read_Info[read_name][signal_id].extend([lt, lp, rt, rp]) # 每条read 每个信号都会保存很多参数
            #print read_name, signal_id, lt, lp, rt, rp
    return Read_Info


def CV_CV(Read_CV, Read_Info, window_size):
    for read_name, CV in Read_CV.items():
        cv_length = len(CV)
        for window_id in range(cv_length - window_size):
            cv_array = np.array(CV[window_id:window_id + window_size])
            mean = np.mean(cv_array)
            SD = np.std(cv_array)
            cv_cv = SD / mean * 100
            signal_id = window_id + window_size + 1
            if signal_id > window_size * 2:
                Read_Info[read_name][signal_id].extend([mean, SD, cv_cv])
    return Read_Info


def Print_Info(Read_Info, window_size):
    print(
        'read_name\tsignal_id\tsignal\tsig_lt\tsig_lp\tsig_rt\tsig_rp\tsig_mean\tsig_sd\tsig_cv\tcv_lt\tcv_lp\tcv_rt\tcv_rp\tcv_mean\tcv_sd\tcv_cv'
    )
    for read_name, Signal in Read_Info.items():
        signal_length = len(Signal)
        for signal_id, Info in Signal.items():
            if signal_id > window_size * 2 and signal_id <= signal_length - window_size * 2:
                Info = [read_name, signal_id] + Info
                Info = [str(info) for info in Info]
                print('\t'.join(Info))


def Write_Output(Read_Region, output):
    with open(output, 'w') as open_out:
        for read_name, Region in Read_Region.items():
            for region in Region:
                open_out.write(read_name + '\t' + str(region[0]) + '\t' + str(
                    region[1]) + '\t' + str(region[1] - region[0] + 1) + '\n')


if __name__ == '__main__':
    opt = Get_Opt()
    Read_Length, Read_Fast5 = Read_Index(opt.input)
    Read_Signal, Read_Info = Read_Fast5_File(Read_Fast5)
    Read_Info = Ttest_Signal(Read_Length, Read_Signal, Read_Info,
                             opt.window_size)
    #Standard = Read_Standard_Array( opt.standard_array )
    Read_CV, Read_Info = Calculate_CV(Read_Signal, Read_Info, opt.window_size)
    Read_Info = Ttest_CV(Read_Length, Read_CV, Read_Info, opt.window_size)
    Read_Info = CV_CV(Read_CV, Read_Info, opt.window_size)
    Print_Info(Read_Info, opt.window_size)
    #Write_Output( Read_Region, opt.output )
