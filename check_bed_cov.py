#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
统计bed区间上的覆盖情况
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
import os
import argparse
import logging

import pandas as pd
import ipdb

import utils

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
    parser = argparse.ArgumentParser(prog='统计bed区间上的覆盖情况')
    parser.add_argument('--proj_path', help='项目目录')
    parser.add_argument('--proj_id', help='项目目录')
    parser.add_argument('--bed_file', help='bed文件,第4列可以作为bed名称')
    parser.add_argument('--bam_file', help='bam文件')
    parser.add_argument('--depth_file', help='深度文件，如果给出深度文件，可以不给bam文件;'
                                             '如果不给深度文件，会采用samtools depth进行计算',
                        default=None)
    parser.add_argument('--cover_list', help='统计不同深度的覆盖情况，逗号（英文）隔开，比如 1,2,3', default=None)
    parser.add_argument('--log', help='log file, default=.log', default='.log')
    parser.add_argument("--quite", help="increase output verbosity",
                         action="store_true")
    if len(sys.argv) == 1:
        parser.print_help()
        exit()
    return parser.parse_args()


def run_samtools_depth(bed_file, bam_file, output):
    cmd = 'samtools depth -b {bed_file} {bam_file} > {output}'.format(
        bed_file=bed_file, bam_file=bam_file, output=output)
    utils.safe_run(cmd, retry=2)


def cover_sum(value, depth=[10]):
#     out_list = []
    return pd.Series([sum(value>=int(ii)) for ii in depth], index=depth)
    #  return pd.DataFrame([sum(value>=int(ii)) for ii in depth])

def get_coverage_depth(df, depth_list):
    return df.iloc[:, 2:].apply(cover_sum, depth=depth_list)


def proc_depth_info(bed_file, depth_file, depth_list):
    depth_df = pd.read_table(depth_file, header=None, prefix='x', dtype={'x0':str})
    depth_df.iloc[:, 0] = depth_df.iloc[:, 0].astype(str)
    bed_dt = utils.read_table(bed_file)
    region_df_list = []
    for xx in bed_dt:
        logger.info(' '.join(xx))
        if len(xx) > 3:
            chrom, start, end, region_name = xx[:4]
        elif len(xx) == 3:
            chrom, start, end = xx
        else:
            logger.error('bed文件输入错误，请核实')
            exit()
        #  ipdb.set_trace()
        region_df = depth_df[(depth_df.iloc[:, 0] == chrom) & (depth_df.iloc[:, 1] >= int(start)) & (depth_df.iloc[:, 1] <= int(end))]
        region_cov_df = get_coverage_depth(region_df, depth_list)
        length = int(end) - int(start) + 1
        region_cov_df = region_cov_df / length
        region_mean_df = region_df.iloc[:, 2:].sum()/length

        region_cov_df.loc['mean_depth', :] = region_mean_df

        region_cov_df.loc['region', :] = '%s:%s-%s' %(chrom, start, end)

        if len(xx) > 3:
            region_cov_df.loc['region_name', :] = region_name
        region_cov_df = region_cov_df.transpose()
        region_df_list.append(region_cov_df)
    #  ipdb.set_trace()
    return pd.concat(region_df_list)


def main():
    args = get_args()
    proj_path = args.proj_path
    proj_id = args.proj_id
    bed_file = args.bed_file
    bam_file = args.bam_file
    depth_file = args.depth_file
    cover_list = args.cover_list
    log_file = args.log
    quite = args.quite
    global logger
    logger = log(log_file, quite=quite)
    if not cover_list:
        depth_list = list(range(1, 40))
    else:
        depth_list = cover_list.split(',')

    if not depth_file:
        depth_file = os.path.join(proj_path, '%s.depth.xls' % proj_id)
        run_samtools_depth(bed_file, bam_file, depth_file)
    depth_df = proc_depth_info(bed_file, depth_file, depth_list)
    output_file = os.path.join(proj_path, '%s_cover.xls' % proj_id)
    #  with open(bed_file) as f: ncol = next(f).split('\t')
    depth_df = depth_df.iloc[:, range(len(depth_list)+1, depth_df.shape[1]) + range(0, len(depth_list) + 1)]
    depth_df = depth_df.round(4)
    depth_df.to_csv(output_file, sep='\t', index=False)
    logger.info('结果保存在%s' % output_file)

if __name__ == '__main__':
    main()
