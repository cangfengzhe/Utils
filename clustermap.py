#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
clustermap
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

import matplotlib
matplotlib.use('Agg')
import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt


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
    parser = argparse.ArgumentParser(prog='clustermap')
    parser.add_argument('--input_file', help='tsv file including values')
    parser.add_argument('--output_prefix', help='the prefix of jpg file')
    # parser.add_argument('--log', help='log file, default=.log', default='.log')
    # parser.add_argument("--quite", help="increase output verbosity",
                         # action="store_true")
    if len(sys.argv) == 1:
        parser.print_help()
        exit()
    return parser.parse_args()

def main():
    args = get_args()
    input_file = args.input_file
    output_prefix = args.output_prefix
    # log_file = args.log
    # quite = args.quite
    # global logger
    # logger = log(log_file, quite=quite)
    df = pd.read_table(input_file)
    # pdb.set_trace()
    # print(df.shape)
    import scipy.spatial as sp, scipy.cluster.hierarchy as hc
    from scipy.cluster.hierarchy import dendrogram
    row_linkage, col_linkage = (hc.linkage(sp.distance.pdist(x), method='average', metric='euclidean')
      for x in (df.values, df.values.T))
    # linkage = hc.linkage(sp.distance.squareform(df), method='average', metric='euclidean')
    # print(col_linkage)
    g = sns.clustermap(df,row_linkage=row_linkage, col_linkage=col_linkage,
        method='average', metric='euclidean', cmap="vlag", figsize=(40, 20))
    # g = sns.clustermap(df, cmap="vlag",method='average',metric='euclidean',
                # row_colors=label_num_color,
                # row_cluster = False,
                    # yticklabels = yticklabels,
                # linewidths=.3,
                # figsize=(13, 13)
                       # )
    # g.ax_col_dendrogram.legend(loc="center", ncol=6)
    g.savefig('%s_clustermap.jpg' % output_prefix, dpi=300)
    fig = plt.figure(figsize=(25, 10))
    # pdb.set_trace()
    labels = df.columns
    dn = dendrogram(col_linkage, labels=labels)
    plt.savefig('%s_tree.jpg' % output_prefix, dpi=300)
if __name__ == '__main__':
    main()
