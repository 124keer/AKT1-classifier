import pandas as pd
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import DataStructs

# Read the CSV file containing known compounds
known_compounds_df = pd.read_csv('Known_compounds.csv')  # Replace 'Known_compounds.csv' with your CSV file name
known_compounds_smiles = known_compounds_df['Canonical_smiles']  # Assuming the column name for SMILES is 'Canonical_smiles'

# Generate MACCS fingerprints for known compounds
known_compounds = [Chem.MolFromSmiles(smiles) for smiles in known_compounds_smiles]
known_maccs_fingerprints = [MACCSkeys.GenMACCSKeys(mol) for mol in known_compounds]

# Create a DataFrame to store the known compound data
known_compounds_data = pd.DataFrame({'Known_Compound_SMILES': known_compounds_smiles,
                                     'Known_MACCS_Fingerprint': known_maccs_fingerprints})

# Save the known compound data as CSV
known_compounds_data.to_csv('maccs_fingerprints_known_kstar11.csv', index=False)  # Replace 'maccs_fingerprints_known_RF.csv' with the desired output file name

# Read the CSV file containing unknown compounds
unknown_compounds_df = pd.read_csv('kstar_similar11.csv')  # Replace 'Unknown_RF.csv' with your CSV file name
unknown_compounds_smiles = unknown_compounds_df['Canonical_smiles']  # Assuming the column name for SMILES is 'Canonical_smiles'
unknown_compounds_chembl = unknown_compounds_df['ChEMBL ID']

# Generate MACCS fingerprints for unknown compounds
unknown_compounds = [Chem.MolFromSmiles(smiles) for smiles in unknown_compounds_smiles]
unknown_maccs_fingerprints = [MACCSkeys.GenMACCSKeys(mol) for mol in unknown_compounds]

# Create a DataFrame to store the unknown compound data
unknown_compounds_data = pd.DataFrame({'Unknown_Compound_SMILES': unknown_compounds_smiles,
                                       'Unknown_MACCS_Fingerprint': unknown_maccs_fingerprints})

# Save the unknown compound data as CSV
unknown_compounds_data.to_csv('maccs_fingerprints_unknown_kstar11.csv', index=False)  # Replace 'maccs_fingerprints_unknown_RF.csv' with the desired output file name

# Calculate the Tanimoto coefficient similarity for each pair of known and unknown compounds
similarities = []
for unknown_fingerprint in unknown_maccs_fingerprints:
    similarity_scores = [DataStructs.TanimotoSimilarity(unknown_fingerprint, known_fingerprint) for known_fingerprint in known_maccs_fingerprints]
    similarities.append(similarity_scores)

# Create a new DataFrame to store the similarity scores
output_df2_similarity = pd.DataFrame(similarities, columns=known_compounds_smiles)
output_df2_similarity.insert(0, 'Unknown_Compound_SMILES', unknown_compounds_chembl)

# Find the closest known compound for each unknown compound
output_df2_similarity['Closest_Known_Compound'] = output_df2_similarity.iloc[:, 1:].idxmax(axis=1)

# Save the output DataFrame as CSV
output_df2_similarity.to_csv('Similar_compounds_maccskstar11.csv', index=False)
