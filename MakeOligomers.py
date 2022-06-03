from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdmolfiles
from openbabel import openbabel as ob
from openbabel import pybel
import os
import os.path

#General settings
conformer = ob.OBConformerSearch()
conformer.SetLogStream(None)
x = ob.OBMinimizingEnergyConformerScore()
conformer.SetScore(x)
y = ob.OBStericConformerFilter()
conformer.SetFilter(y)
pff = ob.OBForceField_FindType("gaff")


#Substrates properties 
iso = Chem.MolFromSmiles('N=C=O')
diiso = Chem.MolFromSmiles('O=C=N.N=C=O')
ol = Chem.MolFromSmiles('CO')
ol_2 = Chem.MolFromSmiles('COCCO')
diol = Chem.MolFromSmiles('OC.CO')

#Substrates and products
a_list = [] #isocyanates and diisocyanates
b_list = [] #mono and poliols
products_list = [] #output

#Reactions in SMARTS
#Dimerization
ab = '[N:1]=[C:2]=[O:3].[C:4][O;H1:5]>>[N:1][C:2]([O:5][C:4])=[O:3]'
#Trimerization
aba = '[N:1]=[C:2]=[O:3].([O;H1:4][C:5].[C:6][O;H1:7]).[N:8]=[C:9]=[O:10]>>([N:1][C:2](=[O:3])[O:4][C:5].[C:6][O:7][C:9](=[O:10])[N:8])'
bab = '[C:1][O;H1:2].([O:3]=[C:4]=[N:5].[N:6]=[C:7]=[O:8]).[C:9][O;H1:10]>>([C:1][O:2][C:4](=[O:3])[N:5].[N:6][C:7](=[O:8])[O:10][C:9])'
#Tetramerization
abab = ('[N:1]=[C:2]=[O:3].([O;H1:4][C:5].[C:6][O;H1:7]).([N:8]=[C:9]=[O:10].[N:11]=[C:12]=[O:13]).[O;H1:14][C:15]>>'
        '([N:1][C:2](=[O:3])[O:4][C:5].[C:6][O:7][C:9](=[O:10])[N:8].[N:11][C:12](=[O:13])[O:14][C:15])')

def MakeOligomers(substrates, size, output=None, molecules_2D=False, molecules_3D=False, conformers=False, images=False):
    
    print("\n /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ \n",
          "- - - - - - - - - - - - - - - - - - - - - - - - -\n",
          "               o l i g o m e r\n",
          "- - - - - - - - - - - - - - - - - - - - - - - - -\n",
          "\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ \n\n",
          "-------------------------------------------------\n",
          "oligomer version 1.0\n",
          "-------------------------------------------------\n"
         )
  
    print("Loading input file...", end=" ")
    if os.path.exists(substrates):
        print("OK.\n")
        
        valid = prepare_reaction(substrates)
        #Perfom reaction
        print('Reaction type: ', end="")
        if valid:
            print('polyurethane fragments synthesis.\n')
            print('Number of units for reaction:', size)
            print('\n------------------------------------------------')
            if size == 2:
                mer='dimer'
                i = perform_dimerization()
            elif size == 3:
                mer='trimer'
                i = perform_trimerization()
            elif size == 4:
                mer='tetramer'
                i = perform_tetramerization()
            else:
                print('Error: Number of units out of range! (2=dimer, 3=trimer, 4=tetramer)')
                exit(1)
        else:
            print('invalid reaction properties')
            exit(1)
        
    else:
        print('\nError: No such file or directory. Wrong input file name!')
        exit(1)
    
    print('\nSynthesis done!')
    print('%i products were generated.\n'%i) 
    print('------------------------------------------------')
   
    process_output(output)
    process_molecules_2D(molecules_2D)
    process_molecules_3D(molecules_3D, mer) 
    process_conformers(conformers, mer) 
    process_images(images) 
    
    print('                     END')
    print('------------------------------------------------\n')

    return

def prepare_reaction(substrates):

    substrates = Chem.SmilesMolSupplier(substrates, delimiter=';', smilesColumn=1, nameColumn=0)

    #Assign properties to the substrates according to the functional group
    n_diiso=0 
    n_iso=0 
    n_ol=0 
    n_diol=0 
    n=0   
    for comp in substrates:
        if comp.HasSubstructMatch(diiso):
            comp.SetProp('func_group', 'diiso')
            a_list.append(comp)
            n_diiso+=1

        elif comp.HasSubstructMatch(iso):
            comp.SetProp('func_group', 'iso')
            a_list.append(comp)
            n_iso+=1

        elif comp.HasSubstructMatch(diol):
            if comp.HasSubstructMatch(ol_2):
                comp.SetProp('func_group', 'ol')
                b_list.append(comp)
                n_ol+=1
            else:
                comp.SetProp('func_group', 'diol')
                b_list.append(comp)
                n_diol+=1

        elif comp.HasSubstructMatch(ol):
            comp.SetProp('func_group', 'ol')
            b_list.append(comp)
            n_ol+=1
        n+=1
               
    print('Number of uploaded substrates: %i\n\n Substrates types:\n - %i isocyanates\n - %i diisocyanates\n - %i alcohols/phenols\n - %i diols\n' 
          % (n, n_iso, n_diiso, n_ol, n_diol))
    return (n_iso != 0 or n_diiso != 0) and (n_ol != 0 or n_diol != 0)

def perform_dimerization():
    print('-_-_-_-_-_-_-_-_ dimerization -_-_-_-_-_-_-_-_-_')
    i = 0
    for a in a_list:
        for b in b_list:
            reaction = ab
            reacts = (a,b)
            rxn = rdChemReactions.ReactionFromSmarts(reaction)
            products = rxn.RunReactants(reacts)
            products_list.append(products[0][0])
            i+=1
            print('\nReaction %i: %s + %s -> %s' %(i,Chem.MolToSmiles(reacts[0]),
                Chem.MolToSmiles(reacts[1]),Chem.MolToSmiles(products[0][0])))
    return i

def perform_trimerization():
    print('-_-_-_-_-_-_-_-_ trimerization -_-_-_-_-_-_-_-_-')
    i = 0
    for a in a_list:
        for b in b_list:
            if a.GetProp('func_group') == 'iso' and b.GetProp('func_group') == 'diol':
                reaction = aba
                reacts = (a,b,a)
                rxn = rdChemReactions.ReactionFromSmarts(reaction)
                products = rxn.RunReactants(reacts)
                products_list.append(products[0][0])
                i+=1
                print('\nReaction %i: %s + %s + %s -> %s' %(i,Chem.MolToSmiles(reacts[0]),
                    Chem.MolToSmiles(reacts[1]),Chem.MolToSmiles(reacts[2]),Chem.MolToSmiles(products[0][0])))

            elif a.GetProp('func_group') == 'diiso' and b.GetProp('func_group') == 'ol':
                reaction = bab
                reacts = (b,a,b)
                rxn = rdChemReactions.ReactionFromSmarts(reaction)
                products = rxn.RunReactants(reacts)
                products_list.append(products[0][0])
                i+=1
                print('\nReaction %i: %s + %s + %s -> %s' %(i,Chem.MolToSmiles(reacts[0]),
                    Chem.MolToSmiles(reacts[1]),Chem.MolToSmiles(reacts[2]),Chem.MolToSmiles(products[0][0])))
    return i

def perform_tetramerization():
    print('-_-_-_-_-_-_-_- tetramerization -_-_-_-_-_-_-_-_')
    i = 0
    for a in a_list:
        if a.GetProp('func_group') == 'diiso':
            for b in b_list:
                if b.GetProp('func_group') == 'diol':
                    reaction = abab
                    reacts = (a,b,a,b)
                    rxn = rdChemReactions.ReactionFromSmarts(reaction)
                    products = rxn.RunReactants(reacts)
                    products_list.append(products[0][0])
                    i+=1
                    print('\nReaction %i: %s + %s + %s + %s -> %s' %(i,Chem.MolToSmiles(reacts[0]),Chem.MolToSmiles(reacts[1]),
                        Chem.MolToSmiles(reacts[2]),Chem.MolToSmiles(reacts[3]),Chem.MolToSmiles(products[0][0])))
    return i


def process_output(output):
    if output:
        print('Products will be saved in SMILES format.')
        with open(output, 'w') as o:
            i = 0
            o.write('Nr\tSMILES\n')
            for product in products_list:
                i+=1
                o.write(str(i)+'\t'+Chem.MolToSmiles(product)+'\n')
        print('Saving file: %s\n'%(output))
        print('------------------------------------------------')

def process_molecules_2D(molecules_2D):
    if molecules_2D:
        print('Generating 2D structures.\n')
        if not os.path.exists(molecules_2D):
            os.makedirs(molecules_2D)
            i=0
            for product in products_list:
                i+=1
                mol_file = molecules_2D+'/'+str(i)+'.mol'
                rdmolfiles.MolToMolFile(product, mol_file)
            print('%i mol files saved in %s\n'%(i,(molecules_2D)))
            print('------------------------------------------------')
        else:
            print('Error: directory %s already exist!\n' %(molecules_2D))

def process_molecules_3D(molecules_3D, mer):
    if molecules_3D:
        print('Generating 3D structures.\n')
        if not os.path.exists(molecules_3D):
            os.makedirs(molecules_3D)
            i=0
            conversion = ob.OBConversion()
            mol = ob.OBMol()
            conversion.SetInAndOutFormats("smi", "mol2")
            gen3d = ob.OBOp.FindType("gen3D")
            main_file = molecules_3D+'/ligands.mol2' #file with all ligands together
            
            with open(main_file,'w') as outfile:
                for product in products_list:
                    i+=1
                    mol_file = molecules_3D+'/'+str(i)+'.mol2'
                    conversion.ReadString(mol, Chem.MolToSmiles(product))
                    gen3d.Do(mol, "--best")
                    print('Converting molecule %i.\n' % (i))
                    pybel.Molecule(mol).title = 'PU_'+mer+'_'+str(i) 
                    conversion.WriteFile(mol, mol_file)
                    with open(mol_file,'r') as infile:
                        outfile.write(infile.read())

                print('%i mol2 files saved in %s\n'%(i,(molecules_3D)))
                print('------------------------------------------------')
        else:
            print('Error: directory %s already exist!\n' % (molecules_3D))

def process_conformers(conformers, mer):
    if conformers:
        print('Generating conformers for 3D structures.')
        
        conversion = ob.OBConversion()
        conversion.SetOutFormat("mol2")
    
        if not os.path.exists(conformers):
            os.makedirs(conformers)
            i=0
            conf_number = 20
            main_file = conformers+'/ligands.mol2'
            with open(main_file,'w') as outfile:
                for product in products_list:
                    i+=1
                    print('\nMolecule %i\n'% i)
                    mol = pybel.readstring("smi", Chem.MolToSmiles(product))
                    mol.addh()
                    mol.make3D()
                    print('Number of rotatable bonds: %d\n' % (mol.OBMol.NumRotors()))
                    conformer.Setup(mol.OBMol, conf_number, 5, 5, 25)
                    conformer.GetConformers(mol.OBMol)
                    # checking if conformers are different
                    dif_conf = 1
                    mol.OBMol.SetConformer(0)
                    init_conf=conversion.WriteString(mol.OBMol)
                    
                    for n in range(1, mol.OBMol.NumConformers()):
                        mol.OBMol.SetConformer(n)
                        new_conf=conversion.WriteString(mol.OBMol)


                        if init_conf != new_conf:
                            dif_conf += 1

                    print('Number of conformers using confab: %d (%i different structures)\n' % (mol.OBMol.NumConformers(), dif_conf))


                    if dif_conf == 1:
                        pff.Setup(mol.OBMol)
                        pff.DiverseConfGen(0.5, 1000000, 50.0, False)
                        pff.GetConformers(mol.OBMol)
                        if mol.OBMol.NumConformers() > conf_number:
                            print('Number of conformers using genetic algorithm: %d (%i will be saved)\n' % (mol.OBMol.NumConformers(),conf_number))
                        else:
                            print('Number of conformers using genetic algorithm: %d\n' % (mol.OBMol.NumConformers()))
                    
                    
                    if mol.OBMol.NumConformers() <= conf_number:
                        for j in range(mol.OBMol.NumConformers()):
                            mol_file = conformers+'/'+str(i)+'_'+str(j)+'.mol2'
                            mol.OBMol.SetConformer(j)
                            pybel.Molecule(mol.OBMol).title = 'PU_'+mer+'_'+str(i)+'_'+str(j)
                            conversion.WriteFile(mol.OBMol, mol_file)
                            with open(mol_file,'r') as infile:
                                outfile.write(infile.read())
                    else:
                        for j in range(conf_number):
                            mol_file = conformers+'/'+str(i)+'_'+str(j)+'.mol2'
                            mol.OBMol.SetConformer(j)
                            pybel.Molecule(mol.OBMol).title = 'PU_'+mer+'_'+str(i)+'_'+str(j)
                            conversion.WriteFile(mol.OBMol, mol_file)
                            with open(mol_file,'r') as infile:
                                outfile.write(infile.read())

            print('\nConformers saved in %s\n'%(conformers))
            print('------------------------------------------------')
        else:
            print('Error: directory %s already exist!\n'%(conformers))

def process_images(images):
    if images:
        print('Generating 2D images.\n')
        conversion = ob.OBConversion()
        mol = ob.OBMol()

        if not os.path.exists(images):
            os.makedirs(images)
            i=0
            conversion.SetInAndOutFormats("smi", "_png2")
            for product in products_list:
                i+=1
                file = images+'/'+str(i)+'.png'
                conversion.ReadString(mol, Chem.MolToSmiles(product))
                conversion.WriteFile(mol, file)

            print('%i images saved in %s\n'%(i,(images)))
            print('------------------------------------------------')

        else:
            print('Error: directory %s already exist!\n'%(images))
