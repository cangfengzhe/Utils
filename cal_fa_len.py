#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
calculate the length of fasta
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
    parser = argparse.ArgumentParser(prog='calculate the length of fasta')
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
    with open(input_file) as f, open(output_file, 'w') as fo:
        for read in f:
            length = len(next(f).rstrip())
            line = '\t'.join([read.rstrip(), str(length)]) + '\n'
            fo.write(line)




if __name__ == '__main__':
    main()
