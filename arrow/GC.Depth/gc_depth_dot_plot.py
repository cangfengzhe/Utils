#!usr/bin/env python

from __future__ import division
import sys
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def listLoop(dataPart, dataAll):
    for data in dataPart:
        dataAll.append(float(data))

def numberStore(dataFile):

    data = {}
    idData = []
    seqId = ''
    with open(dataFile, 'r') as f:
        for line in f:
            start = line[0]
            if 48 <= ord(start) <= 57:
            # 0-9 ASCII : 48-57
                strings = line.strip().split()
                numbers = []
                listLoop(strings, numbers)
                data[seqId] = numbers
            else:
                seqId = line.strip()
                idData.append(seqId)
    return data, idData


def xyExtract(gcData, depthData, ids):

    x, y = [], []
    for idCommon in ids:
        xPart, yPart = gcData[idCommon], depthData[idCommon]
        if len(xPart) != len(yPart):
            print >> sys.stderr, "Error: gc cann't correspond to depth in " + idCommon
            print >> sys.stderr, "gc number is " + str(len(xPart)) + "but depth number is " + str(len(yPart))
            exit(0)
        else:
            listLoop(xPart, x)
            listLoop(yPart, y)
    return x, y


def gcVSdepthBar(x, y, picFile, xmin, xmax, ymin, ymax, xname, yname):

    plt.figure(figsize = (12, 9))
    #plt.xlim(0.2, 0.7)
    #plt.ylim(0, 150)
    #v = [0.2, 0.7, 0, 150]
    v = [xmin, xmax, ymin, ymax]
    plt.axis(v, linewidth = 10)
    plt.plot(x, y, color = 'red', marker = '.', linestyle = 'None')
    # Use when xlable not complete !
    xSite, xTick = [], []
    x = xmin
    while x <= xmax:
        xSite.append(x)
        xTick.append(str(x))
        x += 0.1
    plt.xlabel(xname, fontweight = 'bold', fontsize = 18)
    plt.ylabel(yname, fontweight = 'bold', fontsize = 18)
    plt.xticks(xSite, xTick, fontweight = 'bold')
    plt.yticks(fontweight = 'bold')
    plt.savefig(picFile + '.png', format = 'png')
    plt.savefig(picFile + '.pdf', format = 'pdf')
    plt.savefig(picFile + '.svg', format = 'svg')


def gcDepthOut(x, y, outFile):

    out = 'GC \t depth\n'
    for i in xrange(len(x)):
        out += str(x[i]) + '\t' + str(y[i]) + '\n'
    f = open(outFile, 'w')
    f.write(out)
    f.close()


def main(args):

    gcData, gcId = numberStore(args.gc)
    depthData, depthId= numberStore(args.depth)
    if sorted(gcId) != sorted(depthId):
        gcSet, depthSet = set(gcId), set(depthId)
        gcMore = list(gcSet - depthSet)
        depthMore = list(depthSet - gcSet)
        print >> sys.stderr, 'Error: Ids in gc_distribution are not same as ids in depth_distribution!'
        print >> sys.stderr, 'Number of ids in gc_distribution are ', str(len(gcId))#, 'They are:\n', gcId
        print >> sys.stderr, 'Number indepth_distribution are ', str(len(depthId))#, 'They are:\n', depthId
        if len(gcMore) > 0: print >> sys.stderr, 'Ids in gc but not in depth are:\n', gcMore
        if len(depthMore) > 0: print >> sys.stderr, 'Ids in depth but not in gc are:\n', depthMore
        exit(0)
    else:
        ids = gcId
    x, y = xyExtract(gcData, depthData, ids)
    gcDepthOut(x, y, args.result)
    gcVSdepthBar(x, y, args.picture, args.xmin, args.xmax, args.ymin, args.ymax, args.xname, args.yname)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter, description = '\
Dot diagram of sequencing depth distribution. Png and svg pictures are generated.\
Edit: TY 2014-05-19')
    parser.add_argument('-g', '--gc', metavar = 'file', required = True, help = 'gc distribution')
    parser.add_argument('-d', '--depth', metavar = 'file', required = True, help = 'depth distribution')
    parser.add_argument('-x', '--xmin', metavar = 'float', type = float, default = 0.20, help = 'Minimum value in x axis.')
    parser.add_argument('-X', '--xmax', metavar = 'float', type = float, default = 0.70, help = 'Maximum value in x axis.')
    parser.add_argument('-y', '--ymin', metavar = 'float', type = float, default = 0, help = 'Minimum value in y axis.')
    parser.add_argument('-Y', '--ymax', metavar = 'float', type = float, default = 150, help = 'Maximum value in y axis.')
    parser.add_argument('-n', '--xname', metavar = 'string', default = 'GC', help = 'Title of x axis.')
    parser.add_argument('-N', '--yname', metavar = 'string', default = 'Average depth (X)', help = 'Title of y axis.')
    parser.add_argument('-p', '--picture', metavar = 'picture', default = 'gc_depth', help = 'perfix of picture file name')
    parser.add_argument('-r', '--result', metavar = 'file', default = 'gc_depth.txt', help = 'GC vs depth in file')
    args = parser.parse_args()
    main(args)
