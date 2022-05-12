# oligomer
*version 1.0*

Generator of short oligomers of synthetic polymers (version for polyurethanes).

This tool using SMILES entry generates specified oligomers (dimers, trimers, tetramers) encoded in SMILES.
Products can be saved as 2D structures in .mol files or .png images, 3D structures in .mol2 files with option to generate conformers. 

Used packages:

RDKit https://www.rdkit.org/docs/index.html

Pybel https://openbabel.org/docs/dev/UseTheLibrary/Python_Pybel.html

## Dependencies for Debian/Ubuntu
### Libraries needed to run the script
```
$apt update
$apt -y install --no-install-recommends\
                build-essential\
                ca-certificates\
                cmake\
                git\
                zlib1g-dev\
                libcairo2-dev\
                libboost-dev\
                libboost-program-options-dev\
                libboost-iostreams-dev\
                libboost-regex-dev\
                rapidjson-dev\
                python3-dev\
                libbz2-dev\
                libeigen3-dev\
                libxml2-dev\
                swig\
                lzma\
                python3-rdkit\
                python3-pip\
                librdkit1\
                rdkit-data\
                wget
$pip3 install matplotlib
```
### Build patched openbabel
```
$chmod +x ./build_openbabel.sh
$./build_openbabel.sh
```
## Running the project locally
1. Clone this project locally.
2. Prepare your input with substrates in SMILES format (ex. substrates.txt)
3. Run ```$python oligomer.py -h``` to get the help message.
