from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import Draw
from rdkit.Chem import rdmolfiles
from matplotlib.colors import ColorConverter
import os.path

def MakeOligomers(substrates, size, output=None, output_molecules=None, images=False, features=None):
    
    #substrates - list of substrates as Mol objects 
    #size - int (2=dimer, 3=trimer, 4=tetramer)
        
    #Properties 
    iso = Chem.MolFromSmiles('N=C=O')
    diiso = Chem.MolFromSmiles('O=C=N.N=C=O')
    ol = Chem.MolFromSmiles('CO')
    ol_2 = Chem.MolFromSmiles('COCCO')
    diol = Chem.MolFromSmiles('OC.CO')
    
    #Reactions
    #Dimerisation in SMARTS
    ab = '[N:1]=[C:2]=[O:3].[C:4][O;H1:5]>>[N:1][C:2]([O:5][C:4])=[O:3]'
    #Trimerisation in SMARTS
    aba = '[N:1]=[C:2]=[O:3].([O;H1:4][C:5].[C:6][O;H1:7]).[N:8]=[C:9]=[O:10]>>([N:1][C:2](=[O:3])[O:4][C:5].[C:6][O:7][C:9](=[O:10])[N:8])'
    bab = '[C:1][O;H1:2].([O:3]=[C:4]=[N:5].[N:6]=[C:7]=[O:8]).[C:9][O;H1:10]>>([C:1][O:2][C:4](=[O:3])[N:5].[N:6][C:7](=[O:8])[O:10][C:9])'
    #Tetramerisation in SMARTS
    abab = '[N:1]=[C:2]=[O:3].([O;H1:4][C:5].[C:6][O;H1:7]).([N:8]=[C:9]=[O:10].[N:11]=[C:12]=[O:13]).[O;H1:14][C:15]>>([N:1][C:2](=[O:3])[O:4][C:5].[C:6][O:7][C:9](=[O:10])[N:8].[N:11][C:12](=[O:13])[O:14][C:15])'
    
    a_list = [] #isocyanates and diisocyanates
    b_list = [] #mono and poliols
    
    products_list = [] #output
    
    if os.path.exists(substrates):
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
                   
        print('\n%i substrates uploaded.\n\n --- SUBSTRATES ---\n\n %i isocyanates\n %i diisocyanates\n %i alcohols/phenols\n %i diols\n' 
              % (n, n_iso, n_diiso, n_ol, n_diol))
            
        #Perfom reaction
        if size == 2:
            print('\n--- DIMERIZATION ---')
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

        elif size == 3:
            print('\n--- TRIMERIZATION ---')
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

        elif size == 4:
            print('\n--- TETRAMERIZATION ---')
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

        else:
            print('Error: oligomer size out of range! (2=dimer, 3=trimer, 4=tetramer)')
        
    else:
        print('Error: No such file or directory. Wrong input with substrates!')
    
    if output:
        o = open(output,'w')
        i = 0
        o.write('Nr\tSMILES\n')
        for product in products_list:
            i+=1
            o.write(str(i)+'\t'+Chem.MolToSmiles(product)+'\n')
        o.close()
        print('\nOutput file: %s'%output)
    
    if output_molecules:
        path = os.getcwd()
        if not os.path.exists(path+'/oligomers'):
            os.makedirs(path+'/oligomers')
            i=0
            for product in products_list:
                i+=1
                mol_file = './oligomers/'+str(i)+'.mol'
                rdmolfiles.MolToMolFile(product, mol_file)
            print('\n%i mol files saved in %s'%(i,path+'/oligomers'))
        else:
            print('Error: directory ''./oligomers'' already exist!')

    if images:
        path = os.getcwd()
        if not os.path.exists(path+'/oligomers_images'):
            os.makedirs(path+'/oligomers_images')
            i=0
            for product in products_list:
                i+=1
                file = './oligomers_images/'+str(i)+'.png'
                Draw.MolToFile(product, file, size=(600,600))
            print('\n%i images saved in %s'%(i,path+'/oligomers_images'))
        else:
            print('Error: directory ''./oligomers_images'' already exist!')
    
    print('--- DONE ---')

    return products_list