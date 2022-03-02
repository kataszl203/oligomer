#sudo apt-get install python3-rdkit librdkit1 rdkit-data
#export PYTHONPATH=${PYTHONPATH}:$usr/lib/python3/dist-packages
#pip3 install matplotlib

import argparse
from MakeOligomers import MakeOligomers

parser = argparse.ArgumentParser(description='Generate oligomers from specified substrates')

parser.add_argument('-i','--input', type=str, metavar='', required=True,
                    help='input file with list of substrates in format: Name;SMILES')
parser.add_argument('-s','--size', type=int, metavar='', required=True,
                    help='specify size of the oligomers: 2 (dimers), 3 (trimers), 4 (tetramers)')
parser.add_argument('-o','--output', type=str, metavar='',
                    help='name of output file with list of products in SMILES format')
parser.add_argument('-O','--output_molecules', type=bool, metavar='',
                    help='creating folder with .mol files of products (True/False)')
parser.add_argument('-I','--images', type=bool, metavar='',
                    help='creating folder with png images of products (True/False)')

args = parser.parse_args()

if __name__ == '__main__':
    MakeOligomers(args.input, args.size, args.output, args.output_molecules, args.images)