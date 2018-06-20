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
    parser.add_argument('--prefix', help='')
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
    prefix = args.prefix
    log_file = args.log
    quite = args.quite
    global logger
    logger = log(log_file, quite=quite)
    snp = '%s.SNP.vcf' % prefix
    indel = '%s.InDel.vcf' % prefix
    with open(input_file) as fi, open(snp, 'w') as f_snp, open(indel, 'w') as f_indel:
        for xx in fi:
            if xx.startswith('#'):
                f_snp.write(xx)
                f_indel.write(xx)
                continue
            line = xx.strip().split('\t')
            ref_len = len(line[3])
            alt_len = max(map(len, line[4].split(',')))
            if ref_len == alt_len == 1:
                f_snp.write(xx)
            else:
                f_indel.write(xx)

if __name__ == '__main__':
    main()
