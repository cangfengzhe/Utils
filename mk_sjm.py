#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
生成sjm模版
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


class SJM(object):

    """SJM Template"""

    def __init__(self, job_name, after_job, memory, sched_options, cmd):
        """TODO: to be defined1.

        :job_name: TODO
        :after_job: TODO
        :memory: TODO
        :sched_options: TODO
        :cmd: TODO

        """
        self._job_name = job_name
        self._after_job = after_job
        self._memory = memory
        self._sched_options = sched_options
        self._cmd = cmd



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
    parser = argparse.ArgumentParser(prog='生成sjm模版')
    parser.add_argument('--input_file', help='')
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
    log_file = args.log
    quite = args.quite
    global logger
    logger = log(log_file, quite=quite)


if __name__ == '__main__':
    main()
