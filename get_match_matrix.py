#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
根据已有vcf文件，分析样本（bam）在vcf突变位点上的情况。
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
from collections import defaultdict

import pysam
import pandas as pd
from scipy.special import comb

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
    parser = argparse.ArgumentParser(prog='根据已有vcf文件，分析样本（bam）在vcf突变位点上的情况。')
    parser.add_argument('--vcf_file', help='')
    parser.add_argument('--bam_file', help='')
    parser.add_argument('--output_prefix', help='')

    if len(sys.argv) == 1:
        parser.print_help()
        exit()
    return parser.parse_args()


def get_match_list_by_chrom(vcf_dict, samfile, chrom, start=None, end=None):
    result_df = pd.DataFrame()
    mat = pd.DataFrame()
    result_list = {}
    vcf_pos = 0

    vcf_ref = {}
    vcf_ref['query_name'] = 'ref'
    vcf_ref['chrom'] = chrom
    vcf_alt = {}
    vcf_alt['query_name'] = 'alt'
    vcf_alt['chrom'] = chrom
    for pileupcolumn in samfile.pileup(chrom, start=None, end=None):
        pos = pileupcolumn.pos
        if pos in vcf_dict[chrom]:
            vcf_pos += 1

            ref, alt = vcf_dict[chrom][pos]
            vcf_alt[vcf_pos] = alt
            vcf_ref[vcf_pos] = ref
    #         print("\ncoverage at base %s = %s, ref=%s, alt=%s" %
    #            (pileupcolumn.pos, pileupcolumn.n, ref, alt))
            for pileupread in pileupcolumn.pileups:
                query_name = pileupread.alignment.query_name
                if query_name not in result_list:
                    result_list[query_name] = {}
                    result_list[query_name]['base'] = defaultdict(dict)
#                     result_list[query_name]['base']['read_name'] = query_name
                    result_list[query_name]['snp'] = {}
                    result_list[query_name]['ins'] = {}
                    result_list[query_name]['del'] = {}
                    result_list[query_name]['snp']['T'] = 0
                    result_list[query_name]['ins']['T'] = 0
                    result_list[query_name]['del']['T'] = 0
                    result_list[query_name]['snp']['F'] = 0
                    result_list[query_name]['ins']['F'] = 0
                    result_list[query_name]['del']['F'] = 0
                indel = pileupread.indel
                # 重要注释，针对InDel问题，考虑ALT的最终状态
                if indel > 0: #ins
                    match = 'ins'
                    ins_len = len(alt)
                    query_pos = pileupread.query_position
                    bases = pileupread.alignment.query_sequence[query_pos: query_pos + ins_len ]

                elif indel < 0: #del
                    match = 'del'
                    del_len = len(ref) - len(alt)
                    pairs = pileupread.alignment.get_aligned_pairs(with_seq=True)
                    bases = ''.join([xx[2] for xx in pairs if xx[1] in range(pos, pos + len(ref)) and xx[0] is not None]).upper()
    #                 print ('\tdel base in read %s = %s, query pos is %s ' %
    #                       (pileupread.alignment.query_name,
    #                        bases, pileupread.query_position))

                elif not pileupread.is_del :
                    match = 'snp'
                    # query position is None if is_del or is_refskip is set.
                    bases = pileupread.alignment.query_sequence[pileupread.query_position]
    #                 print ('\tbase in read %s = %s, query pos is %s ' %
    #                       (pileupread.alignment.query_name,
    #                        pileupread.alignment.query_sequence[pileupread.query_position], pileupread.query_position))
                else:
                    match = 'error'
                    bases = 'N'
                if match not in result_list[query_name]:
                    result_list[query_name][match] = {}
                    result_list[query_name][match]['T'] = 0
                    result_list[query_name][match]['F'] = 0

                if alt == bases:
                    result_list[query_name][match]['T'] += 1
                else:
                    result_list[query_name][match]['F'] += 1

                result_list[query_name]['base'][vcf_pos] = bases # '%s_%s' % (bases, match)
    mat = mat.append(pd.Series(vcf_ref), ignore_index=True)
    mat = mat.append(pd.Series(vcf_alt), ignore_index=True)
    for query_name in result_list:
        result_list[query_name]['base']['query_name'] = query_name
        result_list[query_name]['base']['chrom'] = chrom
        mat = mat.append(pd.Series(result_list[query_name]['base']), ignore_index=True)
        result_df = result_df.append(pd.Series(data=[query_name, chrom,
                      result_list[query_name]['snp']['T'],
                      result_list[query_name]['snp']['F'],
                      result_list[query_name]['ins']['T'],
                      result_list[query_name]['ins']['F'],
                      result_list[query_name]['del']['T'],
                      result_list[query_name]['del']['F']],
#                       dtype='int32',
                      ), ignore_index=True)
    result_df.columns=['read_name','chrom', 'snp_t', 'snp_f', 'ins_t', 'ins_f', 'del_t', 'del_f']
    return result_df, mat, result_list


def proc_vcf(vcf_file):
    vcf_dict = {}
    with open(vcf_file) as f:
        for xx in f:
            if xx.startswith('#'):
                continue
            chrom, pos, id_v, ref, alt = xx.rstrip('\n').split('\t')[:5]
            pos = int(pos) -1 # 1-based 转为 0-based
            if chrom not in vcf_dict:
                vcf_dict[chrom] = {}
            if pos not in vcf_dict[chrom]:
                vcf_dict[chrom][pos] = []
            vcf_dict[chrom][pos] += [ref, alt]
    return vcf_dict


def cal_prob(s):
    snp_prob = 0.85
    indel_prob = 0.85
    s = s[-6:]
    n =  sum(s)
    match = sum(s[::2])
    p = comb(n,match) * \
        snp_prob ** s[0] * (1-snp_prob) ** s[1] * \
        indel_prob ** (s[2] + s[4]) * (1-indel_prob) ** (s[3] + s[5])
    return round(p, 5)




def main():
    args = get_args()
    vcf_file = args.vcf_file
    bam_file = args.bam_file
    output_prefix = args.output_prefix

    vcf_dict = proc_vcf(vcf_file)
    import pysam
    samfile = pysam.AlignmentFile(bam_file, "rb" )
    match_df = []
    matrix_df = []
    for chrom in vcf_dict.keys():
        result = get_match_list_by_chrom(vcf_dict, samfile, chrom, start=None, end=None)
        match_df.append(result[0])
        matrix_df.append(result[1].fillna('.'))
    match_out_df = pd.concat(match_df)
    match_out_df['prob'] = match_out_df.apply(cal_prob, axis=1)
    matrix_out_df = pd.concat(matrix_df)
    match_out_df.to_csv('%s_match_stat.xls' % output_prefix, sep='\t', index=False)
    matrix_out_df.to_csv('%s_matrix_stat.xls' % output_prefix, sep='\t', index=False)
    samfile.close()


if __name__ == '__main__':
    main()
