#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
根据reads名称提取bam文件
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
    parser = argparse.ArgumentParser(prog='根据reads名称提取bam文件')
    parser.add_argument('--bam', help='输入bam文件')
    parser.add_argument('--pattern', help='匹配模式-in')
    parser.add_argument('--output_bam', help='输出bam文件')
    parser.add_argument('--log', help='log file, default=.log', default='.log')
    parser.add_argument("--quite", help="increase output verbosity",
                         action="store_true")
    if len(sys.argv) == 1:
        parser.print_help()
        exit()
    return parser.parse_args()

def main():
    args = get_args()
    bam = args.bam
    output_bam = args.output_bam
    pattern = args.pattern
    log_file = args.log
    quite = args.quite
    global logger
    logger = log(log_file, quite=quite)
    with open(bam, 'rb') as f, \
            open(output_bam, 'w') as fo:
        for line in f:
            if line.startswith('@'):
                fo.write(line)
                continue
            read_name = line.split('\t')[0]
            if pattern in read_name:
                print(read_name)
                fo.write(line)



if __name__ == '__main__':
    main()
