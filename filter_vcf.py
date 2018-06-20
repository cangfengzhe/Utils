#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author:	xieshang@grandomics.com
# Date & Time:	2018-01-03 16:46:57
# Description:
import os
import sys
import argparse
import pdb
# Third party imports
try:
    import pandas as pd
except (NameError, ImportError) as E:
    print (E)
    print ('A third party package is missing. Please verify your dependencies')
    sys.exit()

# Local lib import

reload(sys)
sys.setdefaultencoding('utf8')
#extract option
def Get_Opt():
    parser = argparse.ArgumentParser(
        prog='filter_vcf.py',
        usage='%(prog)s [options]',
        description='''Tools for filtering annovar result
        with grandfreq and hgmd.''',
        )
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='input file is a multianno.xls', dest='input')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='output file name', dest='output')
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s 0.0.1')
    options = parser.parse_args()
    return options

def Read_Vcf( input ):
    with open( input, 'r' )as open_vcf:
        Func, Freq, line_id = list(), list(), 0
        for eachline in open_vcf.readlines():
            line_id += 1
            sp = eachline.strip().split('\t')
            if line_id == 1:
                if 'Func.refGene' in sp:
                    col_refgene_f = sp.index('Func.refGene')
                    Func.append(col_refgene_f)
                else:
                    print ('You have not annotated refGene database yet.')
                    sys.exit()
                if 'Func.wgEncodeGencodeBasicV19' in sp:
                    col_genecode_f = sp.index('Func.wgEncodeGencodeBasicV19')
                    Func.append(col_genecode_f)
                else:
                    print ('You have not annotated Gencode database yet.')
                    sys.exit()
                if 'CLINSIG' in sp :
                    col_clinsig = sp.index('CLINSIG')
                else:
                    print ('You have not annotated clinvar database yet.')
                    sys.exit()
                if ( '1000g2015aug_eas' in sp and
                    '1000g2015aug_sas' in sp and
                    '1000g2015aug_eur' in sp and
                    '1000g2015aug_afr' in sp and
                    '1000g2015aug_amr' in sp and
                    '1000g2015aug_all' in sp ):
                    col_kg_eas = sp.index('1000g2015aug_eas')
                    col_kg_sas = sp.index('1000g2015aug_sas')
                    col_kg_eur = sp.index('1000g2015aug_eur')
                    col_kg_afr = sp.index('1000g2015aug_afr')
                    col_kg_amr = sp.index('1000g2015aug_amr')
                    col_kg_all = sp.index('1000g2015aug_all')
                    Freq.extend([col_kg_eas, col_kg_sas, col_kg_eur, col_kg_afr,
                        col_kg_amr, col_kg_all])
                else:
                    print ('You have not annotated 1000g2015aug database yet.')
                    sys.exit()
                if ( 'ExAC_ALL' in sp and
                    'ExAC_AFR' in sp and
                    'ExAC_AMR' in sp and
                    'ExAC_EAS' in sp and
                    'ExAC_FIN' in sp and
                    'ExAC_NFE' in sp and
                    'ExAC_OTH' in sp and
                    'ExAC_SAS' in sp ):
                    col_exac_all = sp.index('ExAC_ALL')
                    col_exac_afr = sp.index('ExAC_AFR')
                    col_exac_amr = sp.index('ExAC_AMR')
                    col_exac_eas = sp.index('ExAC_EAS')
                    col_exac_fin = sp.index('ExAC_FIN')
                    col_exac_nfe = sp.index('ExAC_NFE')
                    col_exac_oth = sp.index('ExAC_OTH')
                    col_exac_sas = sp.index('ExAC_SAS')
                    Freq.extend([col_exac_all, col_exac_afr, col_exac_amr,
                        col_exac_eas, col_exac_fin, col_exac_nfe, col_exac_oth,
                        col_exac_sas])
                else:
                    print ('You have not annotated ExAC database yet.')
                    sys.exit()
                if 'gnomAD_exome_ALL' in sp:
                    if ( 'gnomAD_exome_ALL' in sp and
                        'gnomAD_exome_AFR' in sp and
                        'gnomAD_exome_AMR' in sp and
                        'gnomAD_exome_ASJ' in sp and
                        'gnomAD_exome_EAS' in sp and
                        'gnomAD_exome_FIN' in sp and
                        'gnomAD_exome_NFE' in sp and
                        'gnomAD_exome_OTH' in sp and
                        'gnomAD_exome_SAS' in sp ):
                        col_gnomad_all = sp.index('gnomAD_exome_ALL')
                        col_gnomad_afr = sp.index('gnomAD_exome_AFR')
                        col_gnomad_amr = sp.index('gnomAD_exome_AMR')
                        col_gnomad_asj = sp.index('gnomAD_exome_ASJ')
                        col_gnomad_eas = sp.index('gnomAD_exome_EAS')
                        col_gnomad_fin = sp.index('gnomAD_exome_FIN')
                        col_gnomad_nfe = sp.index('gnomAD_exome_NFE')
                        col_gnomad_oth = sp.index('gnomAD_exome_OTH')
                        col_gnomad_sas = sp.index('gnomAD_exome_SAS')
                        Freq.extend([col_gnomad_all, col_gnomad_afr,
                            col_gnomad_amr, col_gnomad_asj, col_gnomad_eas,
                            col_gnomad_fin, col_gnomad_nfe,
                            col_gnomad_oth, col_gnomad_sas])
                    else:
                        print ('You have not annotated gnomAD database yet.')
                        sys.exit()
                elif 'gnomAD_genome_ALL' in sp:
                    if ( 'gnomAD_genome_ALL' in sp and
                        'gnomAD_genome_AFR' in sp and
                        'gnomAD_genome_AMR' in sp and
                        'gnomAD_genome_ASJ' in sp and
                        'gnomAD_genome_EAS' in sp and
                        'gnomAD_genome_FIN' in sp and
                        'gnomAD_genome_NFE' in sp and
                        'gnomAD_genome_OTH' in sp ):
                        col_gnomad_all = sp.index('gnomAD_genome_ALL')
                        col_gnomad_afr = sp.index('gnomAD_genome_AFR')
                        col_gnomad_amr = sp.index('gnomAD_genome_AMR')
                        col_gnomad_asj = sp.index('gnomAD_genome_ASJ')
                        col_gnomad_eas = sp.index('gnomAD_genome_EAS')
                        col_gnomad_fin = sp.index('gnomAD_genome_FIN')
                        col_gnomad_nfe = sp.index('gnomAD_genome_NFE')
                        col_gnomad_oth = sp.index('gnomAD_genome_OTH')
                        Freq.extend([col_gnomad_all, col_gnomad_afr,
                            col_gnomad_amr, col_gnomad_asj,
                            col_gnomad_eas, col_gnomad_fin,
                            col_gnomad_nfe, col_gnomad_oth])
                    else:
                        print ('You have not annotated gnomAD database yet.')
                        sys.exit()
                else:
                    print ('You have not annotated gnomAD database yet.')
                    sys.exit()
                if 'cg69' in sp:
                    col_cg = sp.index('cg69')
                    Freq.append(col_cg)
                else:
                    print ('You have not annotated cg69 database yet.')
                    sys.exit()
                if 'esp6500siv2_all' in sp:
                    col_esp = sp.index('esp6500siv2_all')
                    Freq.append(col_esp)
                else:
                    print ('You have not annotated esp6500 database yet.')
                    sys.exit()
                # if 'grandfreq' in sp:
                    # col_grandfreq = sp.index('grandfreq')
                    # Freq.append(col_grandfreq)
                # else:
                    # print ('You have not annotated grandfreq database yet.')
                    # sys.exit()
                if 'hgmd_variantType' in sp and 'hgmd_pmid' in sp:
                    col_hgmd_v = sp.index('hgmd_variantType')
                else:
                    print ('You have not annotated hgmd database yet.')
                    sys.exit()
                if 'Qual' in sp:
                    col_qual = sp.index('Qual')
                else:
                    print ('missing Qual column.')
            elif line_id == 2:
                col_info = 0
                for col in sp:
                    if 'AC=' in col:
                        return ( Func, Freq, col_clinsig,
                            col_hgmd_v, col_qual, col_info )
                    col_info += 1
                print ('bad vcf: missing INFO column.')
                sys.exit()

def Filter_Vcf( input, output, Func, Freq, col_clinsig,
               col_hgmd_v, col_qual, col_info ):
    with open( input, 'r' ) as open_vcf, open( output, 'w' ) as open_filter:
        count = 0
        # header = open_vcf.readlines()[0]
        # open_filter.write(header)
        for eachline in open_vcf.readlines():
            count += 1
            if count == 1:
                open_filter.write(eachline)
                continue
            sp = eachline.strip().split('\t')
            INFO = sp[col_info].split(';')
            print_out = True
            # for info in INFO:
                # if 'QD=' in INFO:
                    # if float(info.split('=')[1]) < 2:
                        # print_out = False
                        # break
                # elif 'FS=' in info:
                    # if float(info.split('=')[1]) > 60:
                        # print_out = False
                        # break
                # elif 'MQ=' in info:
                    # if float(info.split('=')[1]) < 40:
                        # print_out = False
                        # break
                # elif 'MQRankSum=' in info:
                    # if float(info.split('=')[1]) < -12.5:
                        # print_out = False
                        # break
                # elif 'ReadPosRankSum=' in info:
                    # if float(info.split('=')[1]) < -8:
                        # print_out = False
                        # break
            # if 'Pathogenic' in sp[col_clinsig] or 'Likely pathogenic' in sp[col_clinsig] or 'VUS' in sp[col_clinsig]:
                # open_filter.write( eachline )
                # continue
            # elif sp[col_hgmd_v] in ['DP', 'DFP', 'FTV', 'DM']:
                # open_filter.write( eachline )
                # continue
            # qual = sp[col_qual]
            # if qual == 'Qual':
                # open_filter.write( eachline )
                # continue
            # else:
                # if float(qual) < 20:
                    # print_out = False
                    # continue
            # for func in Func:
                # func_type = sp[func].split(';')
                # if ( ('exonic' not in func_type) and
                    # ('splicing' not in func_type) ):
                    # print_out = False
                    # continue
            for freq in Freq:
                # pdb.set_trace()
                f = sp[freq]
                if f != '.' and f != 'null':
                    if float(f) > 0.05:
                        print_out = False
                        continue
            if print_out:
                open_filter.write( eachline )

def Tran_excel( output ):
    df = pd.read_csv( output, sep='\t', header=0, index_col=False )
    writer = pd.ExcelWriter( output+'.xlsx' )
    df.to_excel( writer, 'Sheet1' )
    writer.save()

if __name__ == '__main__':
    opt = Get_Opt()
    ( Func, Freq, col_clinsig,
        col_hgmd_v, col_qual, col_info ) = Read_Vcf( opt.input )
    Filter_Vcf( opt.input, opt.output, Func, Freq, col_clinsig,
        col_hgmd_v, col_qual, col_info )
    # Tran_excel( opt.output )
