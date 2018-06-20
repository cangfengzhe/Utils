#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
合并两个survior结果
"""

# ---------
# Change Logs:
#
# ---------

__author__ = 'Li Pidong'
__email__ = 'lipidong@126.com'
__version__ = '0.0.1'
__status__ = 'Dev'

import sys
# reload(sys)
# sys.setdefaultencoding('utf-8')
import argparse
import logging
import pdb


def is_na(info, index):
    if info.split(':')[index] == 'NaN':
        return True
    else:
        return False

def log(file_name, logger_name='lipidong', quite=False):
    logger = logging.getLogger(logger_name)
    formatter = logging.Formatter("%(asctime)s-%(levelname)s-%(message)s",
                                 datefmt="%Y-%m-%d %H:%M:%S")
    handler = logging.FileHandler(file_name)
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    if not quite:
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        console.setFormatter(formatter)
        logger.addHandler(console)
    return logger

def get_args():
    parser = argparse.ArgumentParser(prog='合并两个survior结果')
    parser.add_argument('--input1_file', help='')
    parser.add_argument('--input2_file', help='')
    parser.add_argument('--output_file', help='')
    parser.add_argument('--log', help='log file, default=.log', default='.log')
    parser.add_argument("--quite", help="increase output verbosity",
                         action="store_true")
    if len(sys.argv) == 1:
        parser.print_help()
        exit()
    return parser.parse_args()


def read_vcf(f):
    for line in f:
        if line.startswith('#'):
            continue
        dt = line.rstrip().split('\t')
        yield dt


def main():
    args = get_args()
    input1_file = args.input1_file
    input2_file = args.input2_file
    output_file = args.output_file
    log_file = args.log
    quite = args.quite
    global logger
    logger = log(log_file, quite=quite)
    # import pandas as pd
    # input1_df = pd.read_table(input1_file, sep='\t', comment='#',
                              # prefix='x')
    # input2_df = pd.read_table(input2_file, sep='\t', comment='#',
                              # prefix='x')

    # input1_df.merge(input2_df, how='outer', on='')
    template = ['']
# input1_file 是 多个软件合并的vcf 结果
# input2_file  是多个软件与bionano合并的结果
    with open(input1_file) as f1, \
            open(input2_file) as f2,\
            open(output_file, 'w') as fo:
        input1_dict={}
        for dt in read_vcf(f1):
            chrom = dt[0]
            pos = dt[1]
            sv_type = dt[2][:3]
# 根据前几列的信息走key，最后样本的信息做value，去input2_file中去匹配
            key = '%s_%s_%s' % (chrom, pos, sv_type)
            input1_dict[key] = dt[9:]
        pre_sample_num = len(dt) - 9
        output_list = []

        for dt in read_vcf(f2):
            format_type = dt[8].split(':')
            format_len = len(format_type)
            co_index = format_type.index('CO')
            sample1_info = dt[9] # 这一个信息代表input1_file中所有样本信息
            if is_na(sample1_info, co_index):
                co_template = [':'.join(['NaN'] * format_len)] * pre_sample_num
                dt += co_template
            else:
                ty_index = format_type.index('TY')
                sample1_info_list = sample1_info.split(':')
                chrom, pos = sample1_info_list[co_index].split('-')[0].split('_')
                sv_type = sample1_info_list[ty_index].split(',')[0]
                key = '%s_%s_%s' % (chrom, pos, sv_type)
                if key not in input1_dict:
                    pdb.set_trace()
                dt += input1_dict[key]
                # pdb.set_trace()
            dt.pop(9)
            output_list.append('\t'.join(dt) + '\n')

        fo.write('\t'.join(['#CHROM', 'POS',\
                           'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',\
                           'bionano', 'mummer', 'sniffles', 'pbhoney',  'pbsv']) + '\n')
        for line in output_list:
            fo.write(line)


if __name__ == '__main__':
    main()
