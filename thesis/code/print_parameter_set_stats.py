#!/usr/bin/env python
'''
Created on Aug 2, 2016
@author: henrikahl
'''

import numpy as np
import argparse
import matplotlib.pyplot as plt
from scipy import stats


def parse_arg():
    'Use argparse to parse the arguments from the command line'
    descrip = 'Make a histogram of data from a column (default: first)             \
    in a file. Note that this script cannot handle varying column sizes.         \
    This can be circumvented in the case of recurring time-step or counter         \
    fields (on separate fields) by using awk, for example by sorting out         \
    all rows in infile with more than N fields:\n\n\t\tawk "NF>[N]" >             \
    tmpfile && ./extract_histogram.py [ARGUMENTS] infile > outfile && rm         \
    tmpfile\n\n(with the script in the relevant folder)'

    parser = argparse.ArgumentParser(
        description=descrip, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('files', nargs='+')
    parser.add_argument(
        '-b', '--bins', help='number of bins to use', default=10)
    parser.add_argument('-c', '--col', help='use col (0 indexed)', default=0)
    parser.add_argument(
        '-s', '--skip', help='skip number of initial rows', default=0)
    parser.add_argument('-l', '--log', help='use log spacing on data',
                        required=False, dest='log', action='store_true')
    parser.add_argument('-save', '--save', help='save figure',
                        required=False, dest='save', action='store_true')
    parser.add_argument(
        '-r', '--histrange', help='set range of histogram', nargs=2)
    parser.add_argument(
        '-f', '--filter', help='Read rows with more than FILTER fields (default=0)', default=0)
    parser.add_argument(
        '--size', help='Set dimensions of plot', nargs=2, default=[12, 9])
    parser.add_argument('--single', help='', nargs=2, default=[2, 5])
    parser.add_argument('--xlabel', default='')
    parser.add_argument('--ylabel', default='')
    parser.add_argument('--rows')
    parser.add_argument('--columns')
    parser.set_defaults(relative=False)
    parser.set_defaults(save=False)
    parser.set_defaults(log=False)
    args = parser.parse_args()
    return args


def read_blocks(input_file):
    empty_lines = 0
    blocks = []
    blocks.append([])
    with open(input_file, 'r') as fin:
        for line in fin.readlines()[1:]:
            line = line.rstrip('\n')
            # Check for empty/commented lines
            if not line or line.startswith('#'):
                # If 1st one: new block
                if empty_lines == 0:
                    blocks.append([])
                    empty_lines += 1
                    # Non empty line: add line in current(last) block
            else:
                empty_lines = 0
                blocks[-1].append(line)
    return blocks


def get_avg_stddev(input_file):
    stddevs = []
    for block in read_blocks(input_file):
        if len(block) == 0:
            continue
        templines = []

        ''' Read in our data a bit more properly '''
        for line in block[1:]:
            templines = np.append(templines, line.split("\n"))
        lines = [line.split() for line in templines]

        ''' Take out central zone '''
        lines = sorted(lines, key=lambda l: float(l[4]), reverse=True)
        for line in lines:
            for ii in range(len(line)):
                line[ii] = float(line[ii])
#         lines = [line for line in lines if line[4] > .5]
        # Take all over the mean
        mean = np.mean([line[6] for line in lines])
        lines = [line for line in lines if float(line[6]) > mean]
        for ii, line in enumerate(lines):
            lines[ii] = np.array(lines[ii], dtype='float64')
        lines = np.array(lines)
        stddevs.append(np.std(lines))

    return np.mean(stddevs), np.std(stddevs)


def main(args):

    files = args.files
    files = sorted(files, key=lambda l: int(l.split("_")[-2]), reverse=False)
    counter = 0
#     for file in files:
#         print file

    ''' GO! '''
#     for infile in args.files:
    while (counter + 12 < len(files)):
        printData = [0 for ii in xrange(12)]
        stddevData = [0 for ii in xrange(12)]

        fileBatch = [files[counter + ii] for ii in xrange(12)]
        for ii in xrange(12):
            printData[ii], stddevData[ii] = get_avg_stddev(fileBatch[ii])

        print printData
        with open("../noise_output_LOF", 'a') as chan:
            chan.write('%s\t' % (str(files[int(counter)].split("_")[-2])))
            for ii in xrange(12):
                chan.write('%s\t%s\t' % (printData[ii], stddevData[ii]))
            chan.write("\n")
        counter += 12



############################################################
if __name__ == "__main__":
    args = parse_arg()
    main(args)
