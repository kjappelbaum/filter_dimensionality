# -*- coding: utf-8 -*-

from __future__ import absolute_import

import argparse

from .filter_metal_channels.screening import ScreenMetalChannel

parser = argparse.ArgumentParser(description='Run the screening.')
parser.add_argument('-f', '--folder', help='path to folder with structures', type=str, required=True, dest='folder')
parser.add_argument('-ext',
                    '--extension',
                    help='extension of structure files, default: cif',
                    default='cif',
                    dest='extension',
                    type=str)
parser.add_argument('-o',
                    '--outname',
                    help='name for results file, default: metal_channels_results.csv',
                    default='metal_channels_results.csv',
                    dest='outname',
                    type=str)

parser.add_argument('--njobs', help='maximum number of workers  (i.e. processes)', default='2', dest='njobs', type=int)

parser.add_argument('--featuremode',
                    help='mode for feature computation',
                    default='cheap',
                    choices=['all', 'cheap'],
                    dest='feature_mode',
                    type=str)

parser.add_argument(
    '--other_metals',
    help='list of elements that are considered to be metals in addition to transition metals and lanthanoides',
    default='Al, Ga, In',
    dest='other_metals',
    type=list)


def main():
    arguments = parser.parse_args()
    screener = ScreenMetalChannel.from_folders(arguments.folder, extension=arguments.extension, njobs=arguments.njobs)

    df = screener.run()
    df.to_csv(arguments.outname, index=False)


if __name__ == '__main__':
    main()
