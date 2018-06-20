#!/usr/bin/env python3
import sys
import time

def mummer_bed2vcf(mummerBed,sample_id,out_file):
    with open(mummerBed,"r") as io:
        #header
        io.readline()
        records = io.readlines()
        vcf_lines = []
        counter = 1
        for record in records:
            (reference,ref_start,ref_stop,ID,size,strand,
                    sv_type,ref_gap_size,query_gap_size,
                    query_coordinates,method) = record.strip().split()
            #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO FORMAT
            CHROM = reference
            CHR2 = reference
            POS = ref_start
            END = ref_stop
            ID = str(counter)
            REF = "N"
            #sv_type
            sv_dict = {"Deletion":"DEL",
                    "Insertion":"INS",
                    "Repeat_contraction":"DEL",
                    "Repeat_expansion":"DUP",
                    "Tandem_contraction":"DEL",
                    "Tandem_expansion":"DUP"
                    }
            SVTYPE = sv_dict[sv_type]
            #SUB_SVTRPE = sv_type
            ALT = "<{}>".format(SVTYPE)
            QUAL = "."
            FILTER = "PASS"
            SVLEN = size
            SVMETHOD = "mummer"
            INFO = "SVMETHOD={};CHR2={};END={};SVTYPE={};SVLEN={}".format(SVMETHOD,CHR2,END,SVTYPE,SVLEN)
            FORMAT = "GT"
            SAMPLE = "./."
            vcf_line = "\t".join([CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE]) + "\n"
            vcf_lines.append(vcf_line)
            counter += 1
    vcf_header = """##fileformat=VCFv4.2
##source=mummer
##fileDate={}
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr4,length=190214555>
##contig=<ID=chr5,length=181538259>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr12,length=133275309>
##contig=<ID=chr13,length=114364328>
##contig=<ID=chr14,length=107043718>
##contig=<ID=chr15,length=101991189>
##contig=<ID=chr16,length=90338345>
##contig=<ID=chr17,length=83257441>
##contig=<ID=chr18,length=80373285>
##contig=<ID=chr19,length=58617616>
##contig=<ID=chr20,length=64444167>
##contig=<ID=chr21,length=46709983>
##contig=<ID=chr22,length=50818468>
##contig=<ID=chrX,length=156040895>
##contig=<ID=chrY,length=57227415>
##contig=<ID=chrM,length=16569>
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=TRA,Description="Translocation">
##ALT=<ID=INS,Description="Insertion">
##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	{}""".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),sample_id)
    with open(out_file,"w") as io:
        io.write(vcf_header+"\n")
        io.writelines(vcf_lines)


def main():
    mummer_bed2vcf(sys.argv[1],sys.argv[2],sys.argv[3])


if __name__ == "__main__":
    main()

