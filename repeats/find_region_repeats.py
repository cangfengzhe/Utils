#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
查找指定区域的重复序列
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


UCSC_SIMPLE_REPEAT_FILE='/data/lipidong/database/ucsc/hg19/simpleRepeat.txt'
UCSC_NESTED_REPEAT_FILE='/data/lipidong/database/ucsc/hg19/nestedRepeats.txt'


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
    parser = argparse.ArgumentParser(prog='查找指定区域的重复序列')
    parser.add_argument('--input_region', help='chrXX:start-end')
    parser.add_argument('--output_prefix', help='结果文件, {output_prefix}.{input_region}.simple.xls'
                        '{output_prefix}.{input_region}.nested.xls')
    parser.add_argument('--log', help='log file, default=.log', default='.log')
    parser.add_argument("--quite", help="increase output verbosity",
                         action="store_true")
    if len(sys.argv) == 1:
        parser.print_help()
        exit()
    return parser.parse_args()

def main():

    args = get_args()
    input_region = args.input_region
    output_prefix = args.output_prefix
    log_file = args.log
    quite = args.quite
    global logger
    logger = log(log_file, quite=quite)

    chrom, start_end = input_region.split(':')
    if 'chr' not in chrom:
        chrom = 'chr' + chrom
    import pandas as pd
    start, end = start_end.split('-')
    start, end = int(start), int(end)
    simple_repeat = pd.read_table(UCSC_SIMPLE_REPEAT_FILE, header=None, prefix='x')
    nested_repeat = pd.read_table(UCSC_NESTED_REPEAT_FILE, header=None, prefix='x')

    simple_repeat_result = simple_repeat[(simple_repeat['x1'] == chrom) & \
                                         (((start <= simple_repeat['x2']) & (simple_repeat['x2']  <= end)) | \
                                          ((start <= simple_repeat['x3']) & (simple_repeat['x3']  <= end))
                                          ) ]
    nested_repeat_result = nested_repeat[(nested_repeat['x1'] == chrom) & \
                                         (((start <= nested_repeat['x2']) & (nested_repeat['x2']  <= end)) | \
                                          ((start <= nested_repeat['x3']) & (nested_repeat['x3']  <= end))
                                          ) ]

    # pdb.set_trace()
    # nested_repeat_result = nested_repeat[(nested_repeat['x1'] == chrom) & ((start <= nested_repeat['x2'] <= end) | (start <= nested_repeat['x3'] <= end))]
    simple_repeat_result.to_csv('{output_prefix}.{input_region}.simple.xls'.format(
        output_prefix=output_prefix,
        input_region=args.input_region
    ), index=False, header=None, sep='\t')

    nested_repeat_result.to_csv('{output_prefix}.{input_region}.nested.xls'.format(
        output_prefix=output_prefix,
        input_region=args.input_region
    ), index=False, header=None, sep='\t')


if __name__ == '__main__':
    main()
