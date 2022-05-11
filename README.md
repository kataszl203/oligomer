# oligomer
Generator of short oligomers of synthetic polymers (version for polyurethanes).

This tool using SMILES entry generates specified oligomers (dimers, trimers, tetramers) encoded in SMILES.
Products can be saved as 2D structures in .mol files or .png images, 3D structures in .mol2 files with option to generate conformers. 

Used packages:
RDKit https://www.rdkit.org/docs/index.html
Pybel https://openbabel.org/docs/dev/UseTheLibrary/Python_Pybel.html 

## build patched openbabel
```
chmod +x ./build_openbabel.sh
./build_openbabel.sh
```

