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
    parser = argparse.ArgumentParser(prog='desc')
    parser.add_argument('--input_file', help='')
    parser.add_argument('--sample', help='')
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
    sample = args.sample
    log_file = args.log
    quite = args.quite
    global logger
    logger = log(log_file, quite=quite)
    RE=20
    with open(input_file) as f, \
        open(output_file, 'w') as fo:
        for line in f:
            if line.startswith('#'):
                if '#CHROM' in line:
                    header = line.rstrip().split('\t')
                    header[9] = sample + '.Sniffles'
                    fo.write('\t'.join(header) + '\n')
                    continue
                if '_' in line or '-' in line:
                    continue

                fo.write(line)
                continue
            line = line.rstrip().split('\t')
            chrom = line[0]
            if '_' in chrom:
                continue
            sv_type = line[4]

            info = line[7]
            type_info = info.split(';')
            type_dict = {}

            for xx in type_info:
                if '=' in xx:
                    type_name, type_value = xx.split('=')
                    type_dict[type_name] = type_value
            # pdb.set_trace()
            if ('RE' in type_dict) and ('SVLEN' in type_dict) :
                if '_' in type_dict['CHR2']:
                    continue

                if 'NA' in type_dict['RE'] or 'NA' in type_dict['SVLEN']:
                    fo.write('\t'.join(line) + '\n')
                    continue
                if int(type_dict['RE']) >= RE and abs(int(type_dict['SVLEN'])) < 1000000:
                    fo.write('\t'.join(line) + '\n')
                    continue
                if (sv_type =='<TRA>' ) and (int(type_dict['RE']) >= RE):
                    fo.write('\t'.join(line) + '\n')
if __name__ == '__main__':
    main()
