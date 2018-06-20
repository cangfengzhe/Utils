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
sys.path.append('/data/lipidong/.local/lib/python2.7/site-packages/')
# reload(sys)
# sys.setdefaultencoding('utf-8')
import argparse
import logging
import pdb
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

import venn

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
    parser.add_argument('--output_prefix', help='')
    parser.add_argument('--tools_list', help='default bionano,mummer,sniffles,pbhoney,pbsv',
                        default='bionano,mummer,sniffles,pbhoney,pbsv')
    parser.add_argument('--log', help='log file, default=.log', default='.log')

    parser.add_argument("--quite", help="increase output verbosity",
                         action="store_true")
    if len(sys.argv) == 1:
        parser.print_help()
        exit()
    return parser.parse_args()

def is_na(info, index):
    if info.split(':')[index] == 'NaN':
        return True
    else:
        return False

def main():
    args = get_args()
    input_file = args.input_file
    tools_list = args.tools_list
    output_prefix = args.output_prefix
    log_file = args.log
    quite = args.quite
    global logger
    logger = log(log_file, quite=quite)
    fo_dict = {}
    result_dict = {}
    tools_result_dict = {}
    if tools_list:
        tools_list = tools_list.split(',')
    else:
        tools_list = ['bionano', 'mummer', 'sniffles', 'pbhoney','pbsv']
    # tools_list = ['sniffle', 'pbhoney']
    sv_type_list = ['TRA', 'INS', 'DEL', 'INV', 'DUP']
    for t in tools_list:
        for st in sv_type_list:
            result_dict['%s_%s' %(t, st)] =  []# open('%s_%s_%s.txt' %(output_prefix, t, st), 'w')
            fo_dict['%s_%s' %(t, st)] = open('%s_%s_%s.txt' %(output_prefix, t, st), 'w')
        tools_result_dict[t] = []
    with open(input_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            dt = line.rstrip().split('\t')
            sv_id = dt[2]

            sv_type = sv_id[:3]
            if sv_type not in sv_type_list:
                print('sv type not in sv type list')
            # pdb.set_trace()
            vcf_format = dt[8]
            if 'CO' not in vcf_format.split(':'):
                continue
            co_index = vcf_format.split(':').index('CO')
            # tools_info = dt[9:]
            for index, tools_info in enumerate(dt[9:]):
                if not is_na(tools_info, co_index):
                    result_dict['%s_%s' % (tools_list[index], sv_type)].append(sv_id)
                    tools_result_dict[tools_list[index]].append(sv_id)

    for xx in result_dict:
        line = '\n'.join(result_dict[xx])
        fo_dict[xx].write(line)
        fo_dict[xx].close()

    # venn plot
    exec('venn_plot=venn.venn%s' % len(tools_list))
    for sv_type in sv_type_list:
        labels = venn.get_labels([result_dict['%s_%s' % (tool, sv_type)] for tool in tools_list ]
                                 # fill=['number', 'logic']
                                 )
        fig, ax = venn_plot(labels, names=tools_list)
        fig.set_size_inches(12, 12)
        fig.savefig('%s.%s.png' % (output_prefix, sv_type), bbox_inches='tight')
        plt.close()

    # 全部类型整合
    labels = venn.get_labels([tools_result_dict[tool] for tool in tools_list ])
    fig, ax = venn_plot(labels, names=tools_list)
    fig.savefig('%s.%s.png' % (output_prefix, 'all'), bbox_inches='tight')
    fig.set_size_inches(12, 12)
    plt.close()

if __name__ == '__main__':
    main()
