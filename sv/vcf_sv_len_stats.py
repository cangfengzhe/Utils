#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
desc
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
import matplotlib
matplotlib.use('Agg') # 必须使用，否则会报错tkinter.TclError: couldn't connect to display
# 在plt之前加载 agg
import seaborn as sns
from matplotlib import pyplot as plt
sns.set_style('whitegrid')
import math
import pdb


def read_vcf(f):
    for line in f:
        if line.startswith('#'):
            continue
        dt = line.rstrip().split('\t')
        yield dt



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
    parser = argparse.ArgumentParser(prog='desc')
    parser.add_argument('--input_file', help='')
    parser.add_argument('--output_file', help='')
    parser.add_argument('--log', help='log file, default=.log', default='.log')
    parser.add_argument("--quite", help="increase output verbosity",
                         action="store_true")
    if len(sys.argv) == 1:
        parser.print_help()
        exit()
    return parser.parse_args()

def main():
    args = get_args()
    input_file = args.input_file
    output_file = args.output_file
    log_file = args.log
    quite = args.quite
    global logger
    logger = log(log_file, quite=quite)
    type_list = ['INS', 'DEL', 'DUP', 'TRA', 'INV']
    len_dict = {}
    len_dict['INS'] = []
    len_dict['DEL'] = []
    len_dict['DUP'] = []
    len_dict['TRA'] = []
    len_dict['INV'] = []

    depth_dict = {}
    for type_name in type_list:
        depth_dict[type_name] = []

    with open(input_file) as f:
        for line in read_vcf(f):
            info = line[7]
            type_info = info.split(';')
            type_dict = {}

            for xx in type_info:
                if '=' in xx:
                    type_name, type_value = xx.split('=')
                    type_dict[type_name] = type_value
            # pdb.set_trace()
            if type_dict['SVLEN'] != 'NA':
                len_dict[type_dict['SVTYPE']].append(math.log(abs(int(type_dict['SVLEN'] )), 2 ))
    # import pandas as pd
    fig, ax = plt.subplots()
    for type_name in type_list:
        sns.distplot(len_dict[type_name], hist=False, rug=False, label=type_name)
    # pdb.set_trace()
# the size of A4 paper
    # ax.set_xscale('log2')
    ax.set(xlabel='log2(Length(bp))', ylabel='Ratio')
    fig.set_size_inches(12, 12) # 保存图片尺寸
    plt.savefig('aa.png')



if __name__ == '__main__':
    main()
