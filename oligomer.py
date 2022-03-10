#sudo apt-get install python3-rdkit librdkit1 rdkit-data
#sudo apt install libopenbabel-dev
#sudo apt install python-openbabel
#export PYTHONPATH=${PYTHONPATH}:$usr/lib/python3/dist-packages
#export PYTHONPATH=${PYTHONPATH}:$usr/lib/openbabel
#pip3 install matplotlib

import argparse
from MakeOligomers import MakeOligomers


parser = argparse.ArgumentParser(description='Generate oligomers from specified substrates')

parser.add_argument('-i','--input', type=str, metavar='', required=True,
                    help='input file with list of substrates in format: Name;SMILES')
parser.add_argument('-s','--size', type=int, metavar='', required=True,
                    help='specify size of the oligomers: 2 (dimers), 3 (trimers), 4 (tetramers)')
parser.add_argument('-o','--output', type=str, metavar='',
                    help='name of output file with list of products in SMILES format (default: oligomers.txt)')
parser.add_argument('-m2D','--molecules_2D', type=bool, metavar='',
                    help='creating folder with .mol files of products 2D structures  T/F (default: ''./oligomers_2D'')')
parser.add_argument('-m3D','--molecules_3D', type=bool, metavar='',
                    help='creating folder with .mol files of products  T/F (default: ''./oligomers_3D'')')
parser.add_argument('-I','--images', type=bool, metavar='',
                    help='creating folder with png images of products T/F (default: ''./oligomers_images'')')

args = parser.parse_args()

if __name__ == '__main__':
    MakeOligomers(args.input, args.size, args.output, args.molecules_2D, args.molecules_3D, args.images)