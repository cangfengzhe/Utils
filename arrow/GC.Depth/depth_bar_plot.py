#!/usr/bin/env python

from __future__ import division
import sys
import argparse
import gzip
import zipfile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def xySoap(xyFile, check):

    if zipfile.is_zipfile(xyFile):
        z = zipfile.ZipFile(xyFile)
        for zip_file in z.namelist():
            f = z.open(zip_file, 'r')
    elif xyFile.endswith('.gz'):
        f = gzip.open(xyFile, 'r')
    else:
        f = open(xyFile)
    rawData = {}
    totalBase = 0
    #with open(xyFile, 'r') as f:
    for line in f:
        x, y = line.strip().split()
        rawData.setdefault(int(x), 0)
        rawData[int(x)] += int(y)
        #if int(x) != 0:
        if check or (not check and int(x)!=0): # Modified by TY 2015-03-23
            totalBase += int(y)
    f.close()
    return rawData, totalBase


def xyGet(xyFile, check):

    if zipfile.is_zipfile(xyFile):
        z = zipfile.ZipFile(xyFile)
        for zip_file in z.namelist():
            f = z.open(zip_file, 'r')
    elif xyFile.endswith('.gz'):
        f = gzip.open(xyFile, 'r')
    else:
        f = open(xyFile)
    rawData = {}
    totalBase = 0
    numLine = 0
    for line in f:
        numLine += 1
        if numLine%2 == 0:
            strings = line.strip().split()
            for string in strings:
                number = float(string)
                rawData.setdefault(number, 0)
                rawData[number] += 1
                #if number != 0:
                if check or (not check and number!=0): # Modified by TY 2015-03-23
                    totalBase += 1
    f.close()
    return rawData, totalBase

def xySave(xyFile, column, check):
    if zipfile.is_zipfile(xyFile):
        z = zipfile.ZipFile(xyFile)
        for zip_file in z.namelist():
            f = z.open(zip_file, 'r')
    elif xyFile.endswith('.gz'):
        f = gzip.open(xyFile, 'r')
    else:
        f = open(xyFile)
    rawData = {}
    totalBase = 0
    for line in f:
        number = float(line.strip().split()[column - 1])
        rawData.setdefault(number, 0)
        rawData[number] += 1
        #if number != 0:
        if check or (not check and number!=0): # Modified by TY 2015-03-23
            totalBase += 1
    f.close()
    return rawData, totalBase


def xyExtract(rawData, totalBase, maxStart, step, check):

    xMerge = {}
    xMergeMean, xMergeSite, xMergeLable, yMergeSimple = [], [], [], []
    out = '#Total: ' + str(totalBase) + '\nBlock\tNumber\tPercent(%)\n'
    stepAllNum = int((maxStart + step)/step)
    #print >> sys.stderr, 'stepAllNum', stepAllNum
    start, end = 1, step
    for x in sorted(rawData):
        #if x != 0:
        if check or (not check and x!=0):
            number = int(x/step)
            if x%step == 0 and x!= 0:
                number -= 1
            #start, end = (number) * step + 1, (number + 1)* step 
            start, end = (number) * step, (number + 1)* step # Modified by TY 2015-03-23
            #print >> sys.stderr, 'judge:\t', x, '\t', start, '\t', end
            xMerge.setdefault((start, end), 0)
            xMerge[(start, end)] += rawData[x]
    blockNum, numBlock = 0, 0
    #print >> sys.stderr, xMerge
    for block in sorted(xMerge):
        out += '(%d-%d]\t%d\t%f\n' % (block[0], block[1], xMerge[block], xMerge[block]/totalBase*100)
        #print >> sys.stderr, 'merge', block, xMerge[block]
        if numBlock < stepAllNum - 1:
            #xMergeMean.append(int((block[0] + block[1])/2)) 
            xMergeMean.append(int((block[0] + 1 + block[1])/2)) # Modified by TY 2015-03-23
            if block[0] < maxStart:
                #xMergeLable.append(str(block[0]) + '-' + str(block[1])) 
                xMergeSite.append(block[0]) # Modified by TY 2015-03-23
                xMergeLable.append(str(block[0]))
            yMergeSimple.append(xMerge[block]/totalBase * 100)
        else:
            blockNum += xMerge[block]
        numBlock += 1
    xMergeMean.append(int((maxStart * 2 + step)/2))
    #xMergeLable.append(str(maxStart + 1) + '~') # Modified by TY 2015-03-23
    xMergeSite.append(maxStart)
    xMergeLable.append(str(maxStart))
    xMergeSite.append(maxStart + step) # Modified by TY 2015-03-23
    xMergeLable.append('~') # Modified by TY 2015-03-23
    yMergeSimple.append(blockNum/totalBase * 100)
    return xMergeMean, yMergeSimple, xMergeSite, xMergeLable, out


def depthPlot(xMergeMean, yMergeSimple, xMergeSite, xMergeLable, maxStart, step, picFile):
    
    yBiggest = round( sorted(yMergeSimple)[-1] )
    plt.figure(figsize = (16, 12))
    v = [0, maxStart + step, 0, yBiggest + 1]
    #plt.xlim(0, maxStart + step)
    #plt.ylim(0, yBiggest)
    plt.axis(v, linewidth = 20)
    plt.xlabel('Sequencing depth (x)', fontweight = 'bold', fontsize = 20)
    plt.ylabel('Percent of bases (%)', fontweight = 'bold', fontsize = 20)
    #plt.xticks(xMergeMean, xMergeLable, rotation = 45, fontweight = 'bold') 
    plt.xticks(xMergeSite, xMergeLable, fontweight = 'bold') # Modified by TY 2015-03-23
    plt.yticks(fontweight = 'bold')
    plt.bar(xMergeMean, yMergeSimple, int(step/2), alpha = 0.7, color = 'g', align = 'center')
    plt.savefig(picFile + '.png', format = 'png')
    plt.savefig(picFile + '.pdf', format = 'pdf')
    plt.savefig(picFile + '.svg', format = 'svg')


def depthNum(rawData, resultFile):

    out = 'depth\tnumber\n'
    for depth in sorted(rawData):
        out += str(depth) + '\t' + str(rawData[depth]) + '\n'
    f = open(resultFile, 'w')
    f.write(out)
    f.close()


def main(args):

    # Calculate number of all bases in genome.
    #totalBase = lenRead(args.length)
    if args.soap:
        # Extract depth types and numbers of matched bases from coverge_depth file
        # a type of depth and matched number per line
        # 0    25666811
        # 1      460593
        # 2      358162
        # .       .
        # .       .
        # .       .
        rawData, totalBase = xySoap(args.soap, args.add)
    elif args.depth:
        # Extract depth per base from coverage_depth file
        # This file is like fasta
        # Sequence id in first line
        # depth per base seperated by white space in one line or several lines
        # Contig0
        # 0 2 31 2 3 1 4 24 3 5 ...
        rawData, totalBase = xyGet(args.depth, args.add)
        depthNum(rawData, args.result)
    else:
        rawData, totalBase = xySave(args.list, args.col, args.add)
        depthNum(rawData, args.result)
    # Merge depth type 0 1 2 3 ... 100 ... n
    # to type block 1-10 11-20 21-30 ... 151-160 161-170 ... (n-9)-n by step
    # Trans number of block into percent
    # Store the medium value of type block
    # Decrease number of type block to 1-10 11-20 21-30 ... 141-150 150~
    # Decrease number of type block to [0,10], (10,20], (20,30], (30, 40], ..., (150~] Modified by TY 2015-03-20
    #print >> sys.stderr, 'raw', totalBase, '\n', rawData
    xMergeMean, yMergeSimple, xMergeSite, xMergeLable, out = xyExtract(rawData, totalBase, args.maximum, args.step, args.add)
    f = open(args.result + '.block', 'w')
    print >> f, out
    #print >> sys.stderr, 'Last', xMergeMean, '\n', yMergeSimple, '\n', xMergeLable
    # plot the bar diagram
    depthPlot(xMergeMean, yMergeSimple, xMergeSite, xMergeLable, args.maximum, args.step, args.picture)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter, description = \
'Bar diagram of sequencing depth distribution result After soap.coverage or get_coverage_depth.pl. Must use only one result for plotting! Edit: TY 2014-05-16')
    parser.add_argument('-o', '--soap', metavar = 'file', help = 'depths vs base number in two column generated by soap.coverage -plot')
    parser.add_argument('-d', '--depth', metavar = 'file', help = 'depths distribution per base generated by get_coverage_depth.pl. file can be *.gz or *.zip')
    parser.add_argument('-l', '--list', metavar = 'file', help = 'depths per line.')
    parser.add_argument('-c', '--col', metavar = 'int', type = int, default = 1, help = 'column number of depth.Use only when -l/--list.')
    parser.add_argument('-m', '--maximum', metavar = 'int', type = int, default = 150, help = 'Maximum start in x axis')
    parser.add_argument('-s', '--step', metavar = 'int', type = int, default = 10, help = 'step length for x axis')
    parser.add_argument('-a', '--add', action = 'store_true', default = False, help = 'Count 0 or not. Default: not') # Modified by TY 2015-03-23
    parser.add_argument('-p', '--picture', metavar = 'file', default = 'percent_depth', help = 'prefix of picture file name')
    parser.add_argument('-r', '--result', metavar = 'file', default = 'depth_number.txt', help = 'depth and its number when -d choosed')
    args = parser.parse_args()
    main(args)
