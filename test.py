#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test
"""

# ---------
# Change Logs:
#
# ---------

__author__ = 'Li Pidong'
__email__ = 'pidong.li@genetronhealth.com'
__version__ = '1.0.1'
__status__ = 'Production'

import sys
# reload(sys)
# sys.setdefaultencoding('utf-8')
import argparse
import logging
import pandas

@profile
def log(file_name, logger_name='lipidong', verbose=False):
    logger = logging.getLogger(logger_name)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    handler = logging.FileHandler(file_name)
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    if verbose:
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        logger.addHandler(console)
        logger.setLevel(logging.INFO)
    return logger

def get_args():
    parser = argparse.ArgumentParser(prog='test')
    parser.add_argument('--input_file', help='')
    parser.add_argument('--log', help='log file, default=log.log', default='log.log')
    parser.add_argument("--verbose", help="increase output verbosity",
                         action="store_true")
    if len(sys.argv) == 1:
        parser.print_help()
        exit()
    return parser.parse_args()

def main():
    args = get_args()
    input_file = args.input_file
    log_file = args.log
    verbose = args.verbose
    global logger
    logger = log(log_file, verbose=verbose)
    import pandas

if __name__ == '__main__':
    main()
