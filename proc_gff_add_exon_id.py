#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
为gff文件添加exonid
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
    parser = argparse.ArgumentParser(prog='为gff文件添加exonid')
    parser.add_argument('--input_file', help='')
    parser.add_argument('--output_file', help='')
    parser.add_argument('--log', help='log file, default=.log', default='.log')
    parser.add_argument("--quite", help="increase output verbosity",
                         action="store_true")
    if len(sys.argv) == 1:
        parser.print_help()
        exit()
    return parser.parse_args()

def get_transcript(info):
    for xx in info.split(';'):
        if 'transcript_id' in xx:
            return xx.replace('transcript_id=', '')

def main():
    args = get_args()
    input_file = args.input_file
    output_file = args.output_file
    log_file = args.log
    quite = args.quite
    global logger
    logger = log(log_file, quite=quite)
    exon_id = 0
    with open(input_file) as fi,\
            open(output_file, 'w') as fo:
        for line in fi:
            cols = line.rstrip().split('\t')
            feature = cols[2]
            if feature == 'mRNA':
                exon_id = 0
                # cols[8] += ';exon_id=%s' %(get_transcript(cols[8]))
            if feature == 'exon':
                exon_id += 1
                cols[8] += ';exon_id=%s:%s' %(get_transcript(cols[8]), str(exon_id))
                fo.write('\t'.join(cols) + '\n')

if __name__ == '__main__':
    main()
