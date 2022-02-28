from rdkit import Chem
from rdkit.Chem import rdChemReactions

substrates = Chem.SmilesMolSupplier('substrates.txt', delimiter=';', smilesColumn=1, nameColumn=0)

#giving a flag to each compound
# isocyanate in SMILES has 'N=C=O' or 'O=C=N'
# alcohol in SMILES has 'CO' or 'OC' or '(O)' *except COC

#for comp in substrates:
#    if comp is None: continue
#    print(comp.GetNumAtoms())


## Reactions
# a - isocyanate
# b - alcohol
# c - dialkyl glycol (eter bond)

# order of atoms in structure doesn't matter but order of substrates in reaction matter

#dimerisation
ab = '[N:1]=[C:2]=[O:3].[C:4][O:5]>>[N:1][C:2]([O:5][C:4])=[O:3]'
ba = '[C:1][O:2].[N:3]=[C:4]=[O:5]>>[N:3][C:4]([O:2][C:1])=[O:5]'

#trimerisation
aba = '[N:1]=[C:2]=[O:3].[O:4][C:5][C:6][O:7].[N:8]=[C:9]=[O:10]>>[N:1][C:2](=[O:3])[O:4][C:5][C:6][O:7][C:9](=[O:10])[N:8]'
bab = '[C:1][O:2].([O:3]=[C:4]=[N:5].[N:6]=[C:7]=[O:8]).[C:9][O:10]>>[C:1][O:2][C:4](=[O:3])[N:5][N:6][C:7](=[O:8])[O:10][C:9]' #wrong

#test aba
#rxn = rdChemReactions.ReactionFromSmarts(aba)
#reacts = (substrates[0],substrates[15],substrates[0])

#test bab
rxn = rdChemReactions.ReactionFromSmarts(bab)
reacts = (substrates[12],substrates[4],substrates[12])

products = rxn.RunReactants(reacts)

print(Chem.MolToSmiles(reacts[0]))

print(Chem.MolToSmiles(reacts[0]), '+', Chem.MolToSmiles(reacts[1]), '+', Chem.MolToSmiles(reacts[2]), '->', Chem.MolToSmiles(products[0][0]))

#print(substrates[0].GetNumAtoms())




