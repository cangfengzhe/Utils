#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
heatmap for discrete data
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

import matplotlib
matplotlib.use('Agg') # 必须使用，否则会报错tkinter.TclError: couldn't connect to display
# 在plt之前加载 agg
import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd
from matplotlib.colors import ListedColormap


def get_args():
    parser = argparse.ArgumentParser(prog='heatmap for discrete data')
    parser.add_argument('--input_file', help='输入文件')
    parser.add_argument('--output_file', help='保存图片')
    if len(sys.argv) == 1:
        parser.print_help()
        exit()
    return parser.parse_args()


def main():
    args = get_args()
    input_file = args.input_file
    output_file = args.output_file
    df = pd.read_table(input_file, index_col='#Sam')
    new_df = df.replace('T', 1).replace('F', 0).replace('NG', -1)
    new_df.index.name = 'Sample'
    sns.set(rc={'figure.figsize':(8, 8)})
    ax = sns.heatmap(new_df, cmap=ListedColormap(['grey', 'yellow', 'white']))
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks([-0.667, 0, 0.667])
    colorbar.set_ticklabels(['NG', 'F', 'T'])
    figure = ax.get_figure()
    figure.savefig(output_file, dpi=300)



if __name__ == '__main__':
    main()
