from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdChemReactions

def MakeOligomers(substrates, size):
    
    #substrates - list of Mol
    #size - int (2=dimer, 3=trimer, 4=tetramer - default 2)
        
    #Assigning properties to the substrates according to the functional group
    
    #Properties
    iso = Chem.MolFromSmiles('N=C=O')
    diiso = Chem.MolFromSmiles('O=C=N.N=C=O')
    ol = Chem.MolFromSmiles('CO')
    diol = Chem.MolFromSmiles('OC.CO')
    
    a_list = [] #isocyanates and diisocyanates
    b_list = [] #mono and poliols
    
    for comp in substrates:
        if comp.HasSubstructMatch(diiso):
            comp.SetProp('func_group', 'diiso')
            a_list.append(comp)
            
        elif comp.HasSubstructMatch(iso):
            comp.SetProp('func_group', 'iso')
            a_list.append(comp)

        elif comp.HasSubstructMatch(diol):
            comp.SetProp('func_group', 'diol')
            b_list.append(comp)

        elif comp.HasSubstructMatch(ol):
            comp.SetProp('func_group', 'ol')
            b_list.append(comp)
    
    
    #dimerisation in SMARTS
    ab = '[N:1]=[C:2]=[O:3].[C:4][O;H1:5]>>[N:1][C:2]([O:5][C:4])=[O:3]'

    #trimerisation in SMARTS
    aba = '[N:1]=[C:2]=[O:3].([O;H1:4][C:5].[C:6][O;H1:7]).[N:8]=[C:9]=[O:10]>>([N:1][C:2](=[O:3])[O:4][C:5].[C:6][O:7][C:9](=[O:10])[N:8])'
    bab = '[C:1][O;H1:2].([O:3]=[C:4]=[N:5].[N:6]=[C:7]=[O:8]).[C:9][O;H1:10]>>([C:1][O:2][C:4](=[O:3])[N:5].[N:6][C:7](=[O:8])[O:10][C:9])'

    #tetramerisation in SMARTS
    abab = '[N:1]=[C:2]=[O:3].([O;H1:4][C:5].[C:6][O;H1:7]).([N:8]=[C:9]=[O:10].[N:11]=[C:12]=[O:13]).[O;H1:14][C:15]>>([N:1][C:2](=[O:3])[O:4][C:5].[C:6][O:7][C:9](=[O:10])[N:8].[N:11][C:12](=[O:13])[O:14][C:15])'

    products_list = []
    
    if size == 2:
        print('--- DIMERIZATION ---')
        i = 0
        for a in a_list:
            for b in b_list:
                reaction = ab
                reacts = (a,b)
                rxn = rdChemReactions.ReactionFromSmarts(reaction)
                products = rxn.RunReactants(reacts)
                products_list.append(products[0][0])
                i+=1
                print('Reaction',i ,' :', Chem.MolToSmiles(reacts[0]), '+', Chem.MolToSmiles(reacts[1]), '->', Chem.MolToSmiles(products[0][0]))
    
    elif size == 3:
        print('--- TRIMERIZATION ---')
        i = 0
        for a in a_list:
            for b in b_list:
                if a.GetProp('func_group') == 'iso' and b.GetProp('func_group') == 'diol':
                    print(Chem.MolToSmiles(a), Chem.MolToSmiles(b))
                    reaction = aba
                    reacts = (a,b,a)
                    rxn = rdChemReactions.ReactionFromSmarts(reaction)
                    products = rxn.RunReactants(reacts)
                    products_list.append(products[0][0])
                    i+=1
                    print('Reaction',i ,' :', Chem.MolToSmiles(reacts[0]), '+', Chem.MolToSmiles(reacts[1]), '+', Chem.MolToSmiles(reacts[2]), '->', Chem.MolToSmiles(products[0][0]))
                
                elif a.GetProp('func_group') == 'diiso' and b.GetProp('func_group') == 'ol':
                    reaction = bab
                    reacts = (b,a,b)
                    rxn = rdChemReactions.ReactionFromSmarts(reaction)
                    products = rxn.RunReactants(reacts)
                    products_list.append(products[0][0])
                    i+=1
                    print('Reaction',i ,' :', Chem.MolToSmiles(reacts[0]), '+', Chem.MolToSmiles(reacts[1]), '+', Chem.MolToSmiles(reacts[2]), '->', Chem.MolToSmiles(products[0][0]))
        

    elif size == 4: #and a and b has 2 functional groups each
        print('--- TETRAMERIZATION ---')
        reaction = abab
        reacts = (a,b,a,b)
    else:
        print('Error - oligomer size out of range!')
                                     
    return products_list