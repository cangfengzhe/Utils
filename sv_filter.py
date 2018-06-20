#!/usr/bin/env python3
"""
a simple sv filter
filter by support reads number, chromosome
"""
import sys
import sv_vcf

def sv_filter(vcf,min_re,out_name):
    with open(vcf,"r") as io:
        lines = io.readlines()

        # main chrom
        main_chrom = [str(i) for i in range(1,23)] + ["X","Y"]

        # output bedlines
        bed_lines = []

        # check if id is unique
        id_dict = {}
        #chrom1,pos1,chrom2,pos2
        previous_sv_breakpoint = ["NA","NA","NA","NA"]
        for line in lines:
            #skip comment lines
            if line.strip()[0] == "#":
                continue
            sv = sv_vcf.SV(line)
            if sv.id not in id_dict:
                id_dict[sv.id] = 1
            else:
                raise RuntimeError("Duplicated SV ID in you VCF "
                "file {}".format(sv.id))

            # filter
            if int(sv.re) < min_re:
                continue
            if sv.chrom1 not in main_chrom or sv.chrom2 not in main_chrom:
                continue

            sv_breakpoint = [sv.chrom1,sv.pos1,sv.chrom2,sv.pos2]
            # remove duplicate adjacency BND record in picky vcf
            # Exactly the same
            if ((sv.svtype == "TRA" or sv.svtype == "INV") and
                    sv_breakpoint[:4] == previous_sv_breakpoint[:4]):
                continue
            # just swap breakpoint1 and breakpoint2, still the same
            if ((sv.svtype == "TRA" or sv.svtype == "INV") and
                    sv_breakpoint[:2] == previous_sv_breakpoint[2:] and
                    sv_breakpoint[2:] == previous_sv_breakpoint[:2]):
                previous_sv_breakpoint = sv_breakpoint
                continue

            previous_sv_breakpoint = sv_breakpoint
            # convert to bed format
            # chrom,start,end,svtype,id,svlen,re,info
            if sv.chrom1 == sv.chrom2:
                if int(sv.pos1) > 1:
                    sv.pos1 = str(int(sv.pos1)-1)
                bed_line = "\t".join([sv.chrom1,sv.pos1,sv.pos2,sv.svtype,
                    sv.id,sv.svlen,sv.re,sv.info])+"\n"
            else: #TRA
                if int(sv.pos1) > 1:
                    pos1_1 = str(int(sv.pos1)-1)
                    pos1_2 = sv.pos1
                if int(sv.pos2) > 1:
                    pos2_1 = str(int(sv.pos2)-1)
                    pos2_2 = sv.pos2
                bed_line1 = "\t".join([sv.chrom1,pos1_1,pos1_2,sv.svtype,
                    sv.id+"_1",sv.svlen,sv.re,sv.info])+"\n"
                bed_line2 = "\t".join([sv.chrom2,pos2_1,pos2_2,sv.svtype,
                    sv.id+"_2",sv.svlen,sv.re,sv.info])+"\n"
                bed_line = bed_line1+bed_line2
            bed_lines.append(bed_line)
    with open(out_name,"w") as io:
        io.writelines(bed_lines)

def main():
    #test
    sv_filter(sys.argv[1],int(sys.argv[2]),sys.argv[3])

if __name__ == "__main__":
    main()

