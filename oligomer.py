
#sudo apt-get install python3-rdkit librdkit1 rdkit-data
#pip3 install matplotlib

import argparse

import os
import sys

external_lib_path = os.path.dirname(os.path.realpath(__file__)) + '/external/build/lib/python/site-packages'

if os.path.exists(external_lib_path):
    sys.path.insert(0, external_lib_path)

from MakeOligomers import MakeOligomers
from argparse import RawTextHelpFormatter

description = """
/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
- - - - - - - - - - - - - - - - - - - - - - - - -
               o l i g o m e r
- - - - - - - - - - - - - - - - - - - - - - - - -
\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ 
-------------------------------------------------
oligomer version 1.0
-------------------------------------------------

Generate short fragments of polymers.
"""

parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)

parser.add_argument('-i','--input', type=str, metavar='', required=True,
                    help='input file with list of substrates in format: Name;SMILES')
parser.add_argument('-s','--size', type=int, metavar='', required=True,
                    help='number of substrates units used for reaction: 2 (dimers), 3 (trimers), 4 (tetramers)')
parser.add_argument('-o','--output', type=str, metavar='',
                    help='name of output file with list of products in SMILES format')
parser.add_argument('-m2D','--molecules_2D', type=bool, metavar='',
                    help='creating folder with .mol files of products 2D structures  T/F (default: ''./oligomers_2D'')')
parser.add_argument('-m3D','--molecules_3D', type=bool, metavar='',
                    help='creating folder with .mol2 files of products T/F (default: ''./oligomers_3D'')')
parser.add_argument('-c','--conformers', type=bool, metavar='',
                    help='creating folder with conformers .mol2 files of products T/F (default: ''./oligomers_conformers'')')
parser.add_argument('-I','--images', type=bool, metavar='',
                    help='creating folder with png images of products T/F (default: ''./oligomers_images'')')


args = parser.parse_args()

if __name__ == '__main__':
    MakeOligomers(args.input, args.size, args.output, args.molecules_2D, args.molecules_3D, args.conformers, args.images)
