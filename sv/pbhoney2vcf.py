#!/usr/bin/env python3
import sys
import time

def pbhoney2vcf(spot,tail,sample_id,out_file):
    main_chr = ["chr"+str(i) for i in range(1,22)] + ["chrX","chrY"]
    #read spot
    vcf_lines = []
    counter = 1
    with open(spot,"r") as io:
        ##CHROM  START       END   TYPE    SIZE    INFO
        #header
        io.readline()
        records = io.readlines()
        for record in records:
            (CHROM,START,END,TYPE,SIZE,SPOT_INFO) = record.strip().split("\t")

            CHR2 = CHROM
            POS = START
            ID = str(counter)
            REF = "N"
            SVTYPE = TYPE
            ALT = "<{}>".format(SVTYPE)
            QUAL = "."
            FILTER = "PASS"
            SVLEN = SIZE
            SVMETHOD = "pbhoney_spots"

            INFO_fields = SPOT_INFO.split(";")
            INFO_dict = {}
            for info in INFO_fields:
                if "=" in info:
                    info_id,info_value = info.split("=")
                    INFO_dict[info_id] = info_value
                else:
                    INFO_dict[info] = info

            #support reads number
            RE = float(INFO_dict["szCount"])

            #filter
            if int(SVLEN) < 30 or int(SVLEN) > 1000000 or RE < 30:
                continue

            if CHROM not in main_chr or CHR2 not in main_chr:
                continue

            INFO = "SVMETHOD={};CHR2={};END={};SVTYPE={};SVLEN={}".format(SVMETHOD,CHR2,END,SVTYPE,SVLEN)
            FORMAT = "GT"
            SAMPLE = "./."

            vcf_line = "\t".join([CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE])+"\n"
            vcf_lines.append(vcf_line)
            counter += 1

    with open(tail,"r") as io:
        #id     chrKey  uRef    uBreak  uMapq   dRef    dBreak  dMapq   remainSeq       annot   numReads        numZMWs evidence
        #header
        io.readline()
        io.readline()
        records = io.readlines()
        for record in records:
            (sv_id,chrKey,uRef,uBreak,uMapq,dRef,dBreak,dMapq,remainSeq,
                    annot,numReads,numZMWs,evidence) = record.strip().split("\t")

            CHROM = uRef
            CHR2 = dRef
            POS = uBreak
            ID = str(counter)
            REF = "N"
            SVTYPE = annot
            if SVTYPE == "TLOC":
                SVTYPE = "TRA"
            ALT = "<{}>".format(SVTYPE)
            QUAL = "."
            FILTER = "PASS"
            if SVTYPE == "DEL" or SVTYPE == "INV":
                if CHROM != CHR2:
                    continue
                SVLEN = int(dBreak) - int(uBreak)
                if SVLEN < 0:
                    raise RuntimeError("file format error")
            else:
                SVLEN = "NA"


            if SVLEN != "NA":
                if SVLEN < 30 or SVLEN > 1000000:
                    continue

            SVLEN = str(SVLEN)

            #support reads num
            RE = int(numReads)
            if RE < 20:
                continue

            if CHROM not in main_chr or CHR2 not in main_chr:
                continue

            SVMETHOD = "pbhoney_tails"
            INFO = "SVMETHOD={};CHR2={};END={};SVTYPE={};SVLEN={}".format(SVMETHOD,CHR2,END,SVTYPE,SVLEN)
            FORMAT = "GT"
            SAMPLE = "./."
            vcf_line = "\t".join([CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE]) + "\n"
            vcf_lines.append(vcf_line)
            counter += 1


    vcf_header = """##fileformat=VCFv4.2
##source=pbhoney
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
    pbhoney2vcf(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])


if __name__ == "__main__":
    main()

