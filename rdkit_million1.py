# pip install rdkit-pypi
# pip install mordredpy
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
import pandas as pd
from tqdm import tqdm
import os

# Load the data
# url = 'https://raw.githubusercontent.com/gashawmg/molecular-descriptors/main/Orbital_Energies_input_data.csv'
dataset = pd.read_csv('chembl.csv')

# Drop rows with missing SMILES and remove duplicates
df2 = dataset.dropna(subset=['Canonical_smiles'])
df2 = df2.drop_duplicates(subset=['Canonical_smiles'])

# Define a function to calculate RDKit descriptors for a list of SMILES
def RDkit_descriptors(smiles):
    mols = [Chem.MolFromSmiles(i) for i in smiles] 
    calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
    desc_names = calc.GetDescriptorNames()
    
    Mol_descriptors = []
    for mol in mols:
        # add hydrogens to molecules
        mol = Chem.AddHs(mol)
        # Calculate all 200 descriptors for each molecule
        descriptors = calc.CalcDescriptors(mol)
        Mol_descriptors.append(descriptors)
    return Mol_descriptors, desc_names 

# Split the SMILES into chunks of 100,000 for faster processing
chunk_size = 100000
chunks = [df2[i:i+chunk_size] for i in range(0, len(df2), chunk_size)]

total_chunks=len(chunks)
total_time=0

# Check if there is an existing output file
if os.path.isfile('chembl_set.csv'):
    existing_data = pd.read_csv('chembl_set.csv', index_col=0)
else:
    existing_data = pd.DataFrame()

# Calculate descriptors for each chunk and concatenate the results
for i, chunk in enumerate(tqdm(chunks, desc='Processing', total=len(chunks))):
    # Check if this chunk has already been processed
    if len(existing_data) >= len(chunk):
        continue
    # Calculate descriptors for this chunk
    descriptors, desc_names = RDkit_descriptors(chunk['Canonical_smiles'])
    # Convert the descriptors to a dataframe
    df_with_200_descriptors = pd.DataFrame(descriptors, columns=desc_names)
    # Add the chunk index as a new column
    df_with_200_descriptors['chunk_index'] = i
    # Append the data to the existing data
    existing_data = pd.concat([existing_data, df_with_200_descriptors], axis=0)
    # Save the data after each chunk
    existing_data.to_csv('chembl_set.csv')

# Save the final data
existing_data.to_csv('chembl_set.csv')