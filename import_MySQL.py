#!/usr/bin/env python
# coding: utf-8

"""
import SNV/INDEL/CNV/SV info to Mysql

"""

# ---------
# Change Logs:
# 2017-01-06 Friday 14:31 - add sv
# 2017-01-20 Friday 14:33 - revise cnvtable
# ---------

__author__ = 'Li Pidong'
__email__ = 'pidong.li@126.com'
__version__ = '1.0.1'
__status__ = 'production'

import sys
# reload(sys)
# sys.setdefaultencoding('utf-8')
import argparse
import logging
import os

import xlrd
import MySQLdb


def log(file_name, logger_name='root', verbose=False):
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


def open_excel(file):
    try:
        data = xlrd.open_workbook(file)
        return data
    except Exception, e:
        print str(e)


def get_sample_id(smp_id):
    cursor.execute('select id from sample_info where sample_id="%s"' % smp_id)
    dt = cursor.fetchone()
    if dt:
        return dt[0]

def insert_null(table, sample_id, panel):
    cursor.execute('insert into %s(\
                            sample_id, panel) values (\
                            "%s", "%s") ' %
                           (table, sample_id, panel))
    conn.commit()


def get_snp_indel(file_name, panel, tsample):
    data = open_excel(file_name)
    table = data.sheets()[0]
    nrows = table.nrows
    sample_id = get_sample_id(tsample)
    if not sample_id:
        print('The sample id does not exist in MySQL db')
        return None

    # cursor.execute('delete from sample_snp_indel_info\
    # where sample_id="%s" and panel="%s"' % (sample_id, panel))

    for ii in range(1, nrows):
        gene_name = table.cell(ii, 5).value
        refseq_id = table.cell(ii, 12).value
        chrome = str(table.cell(ii, 7).value).split('.')[0]
        start = table.cell(ii, 8).value
        end = table.cell(ii, 9).value
        cDNA_change = table.cell(ii, 13).value
        aa_change = table.cell(ii, 14).value
        mut_type = table.cell(ii, 15).value
        t_freq = table.cell(ii, 17).value

        try:

            cursor.execute('insert into sample_snp_indel_info(\
                            sample_id, panel, gene_name, refseq_id, chrome,\
                            start,  end, cDNA_change, aa_change, mut_type,\
                            t_freq) values (\
                            "%s", "%s", "%s", "%s", "%s", "%s", "%s",\
                            "%s", "%s", "%s", "%s") ' %
                           (sample_id, panel, gene_name, refseq_id,
                            str(chrome), start, end, cDNA_change,
                            aa_change, mut_type, t_freq))

            conn.commit()
        except:
            conn.rollback()
    if nrows == 1:
         insert_null('sample_snp_indel_info', sample_id, panel)

    print('Import {panel} {nrows} SNV/INDEL to DataBase'.
          format(panel=panel, nrows=nrows - 1))


def get_cnv(file_name, panel, tsample):
    data = open_excel(file_name)
    table = data.sheets()[0]
    nrows = table.nrows
    sample_id = get_sample_id(tsample)
    if not sample_id:
        print('The sample id does not exist in MySQL db')
        return None
    # 如果该样本、panel 之前存在， 先删掉
    # cursor.execute('delete from sample_cnv_info \
    #                 where sample_id="%s" and panel="%s"' %
    #                (sample_id, panel))
    for ii in range(1, nrows):
        gene_name = table.cell(ii, 0).value
        chrome = str(table.cell(ii, 1).value).split('.')[0]
        start = table.cell(ii, 2).value
        end = table.cell(ii, 3).value
        fold = table.cell(ii, 4).value
        cnv_type = table.cell(ii, 5).value

        try:
            cursor.execute('insert into sample_cnv_info(\
            sample_id, panel, gene_name, chrome, start, end, fold,\
            cnv_type) values("%s", "%s", "%s", "%s", "%s", "%s",\
            "%s", "%s") ' %
                           (sample_id, panel, gene_name,
                            str(chrome), start, end, fold, cnv_type))

            conn.commit()

        except:
            conn.rollback()

    print('Import {panel} {nrows} CNV to DataBase'.
          format(panel=panel, nrows=nrows - 1))


def get_sv(file_name, panel, tsample):
    data = open_excel(file_name)
    table = data.sheets()[0]
    nrows = table.nrows
    sample_id = get_sample_id(tsample)
    if not sample_id:
        print('The sample id does not exist in MySQL db')
        return None
    # 如果该样本、panel 之前存在， 先删掉

    for ii in range(1, nrows):
        gene_name = table.cell(ii, 0).value.strip()

        if gene_name == '':
            continue

        break_pos = table.cell(ii, 1).value
        extron_pos = table.cell(ii, 2).value
        freq = table.cell(ii, 3).value
        try:

            cursor.execute(
                'insert into sample_sv_info(\
                sample_id, panel, gene_name, break_pos, extron_pos, freq)\
                values("%s", "%s", "%s", "%s", "%s", "%s") ' %
                (sample_id, panel, gene_name, break_pos, extron_pos, freq))
            conn.commit()
        except Exception as e:
            print('ERROR', e)
            conn.rollback()
    print('Import {panel} {nrows} SV to DataBase'.
          format(panel=panel, nrows=nrows - 1))

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        prog='import SNV/INDEL/CNV/SV info to MySQL DataBase ')

    parser.add_argument('--tsample', help='sample name')
    parser.add_argument("--panel", help="panel name/ panel203、 panel51、 panel509、panel88、exome")
    # parser.add_argument("--muttable", help="mutation info file")

    # parser.addHandler('--log', help='log file, default=log.log')
    args = parser.parse_args()
    tsample = args.tsample
    panel = args.panel
    if len(sys.argv) == 1:
        parser.print_help()
        exit()

    if ('ct' in panel) or ('CT' in panel):

        panel = 'CT_DNA'
    # mut_file = args.muttable
    # log_file = args.log
    conn = MySQLdb.connect(host='172.16.34.2',
                           user='lipidong',
                           passwd='lipidong',
                           db='genetron',
                           charset="utf8")
    cursor = conn.cursor()

    sample_id = get_sample_id(tsample)
    if panel == 'exome':
        cursor.execute('delete from sample_snp_indel_info\
            where sample_id="%s" and panel="%s" ' % (
                sample_id, 'WES'))

        cursor.execute('delete from sample_cnv_info\
            where sample_id="%s" and panel="%s" ' % (
                sample_id, 'WES'))
        cursor.execute('delete from sample_sv_info\
            where sample_id="%s" and panel="%s" ' % (
                sample_id, 'WES'))

        exome_mut_table = '{}.exome.mutTable.xls'.format(tsample)
        get_snp_indel(exome_mut_table, 'WES', tsample)
        exome_cnv_table = '{}.exome.cnvTable.xls'.format(tsample)
        if os.path.exists(exome_cnv_table):
            get_cnv(exome_cnv_table, 'WES', tsample)
        panel = 'WES_panel509'

    cursor.execute('delete from sample_snp_indel_info\
        where sample_id="%s" and panel="%s" ' % (
        sample_id, panel))

    cursor.execute('delete from sample_cnv_info\
        where sample_id="%s" and panel="%s" ' % (
        sample_id, panel))
    cursor.execute('delete from sample_sv_info\
        where sample_id="%s" and panel="%s" ' % (
        sample_id, panel))

    mut_table = '{}.mutTable.xls'.format(tsample)
    get_snp_indel(mut_table, panel, tsample)
    cnv_table = '{}.cnvTable.xls'.format(tsample)
    if os.path.exists(cnv_table):
        get_cnv(cnv_table, panel, tsample)
    sv_table = '{}.svTable.xls'.format(tsample)
    if os.path.exists(sv_table):
        get_sv(sv_table, panel, tsample)

    cursor.close()
    conn.close()
