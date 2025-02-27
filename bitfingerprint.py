# Libraries

import rdkit
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem as AChem

# Smiles reading

smiles_list = pd.read_csv('molecule_data.csv')
smiles_list = smiles_list[['SMILES']]

# Fingerprint parameter

fpgen = AChem.GetMorganGenerator(radius=3, fpSize=4096) # Larger fpSize means less chance of colliding bits

# Generates the fingerprints and array

mol = [Chem.MolFromSmiles(smiles) for smiles in smiles_list['SMILES']]	
fpx = [fpgen.GetFingerprint(mols) for mols in mol]
arr = np.array(fpx)

# Prints fingerprints and the array

for i, fingerprint in enumerate(fpx):
	print(f"Molecule {i+1} Fingerprint`: {fingerprint}")

print(arr)
print(arr.shape)

# Exports the array to a csv

df = pd.DataFrame(arr)
df.to_csv("bit_array_fingerprint.csv", header=False, index=False)