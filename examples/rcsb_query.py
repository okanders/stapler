from rcsbsearchapi.search import AttributeQuery, Attr
from Bio import PDB
import os
import random

# Query components
resolution_query = Attr("rcsb_entry_info.resolution_combined") <= 2.0
disulfide_bond_query = Attr("rcsb_entry_info.disulfide_bond_count") > 0

# Combine queries
combined_query = resolution_query & disulfide_bond_query

# Execute and shuffle results
results = list(combined_query.exec())  # Convert iterator to list
random.shuffle(results)  # Shuffle the list of results

# Setup directory
pdb_dir = './native_disulfide_pdb_NONSQLITE_BIG'
if not os.path.exists(pdb_dir):
    os.makedirs(pdb_dir)

# Initialize and download PDBs
pdbl = PDB.PDBList()
for result in results[:30000]:  # Limit to the first 30,000 shuffled results
    print('result', result)
    # Retrieve the PDB file
    try:
        filename = pdbl.retrieve_pdb_file(result, pdir=pdb_dir, file_format='pdb')
        if filename:  # Checks if the file was downloaded
            # Construct new filename by directly using the result variable
            new_filename = os.path.join(pdb_dir, result.upper() + '.pdb')  # result is converted to uppercase
            # Check if the original file exists before renaming
            if os.path.exists(filename):
                os.rename(filename, new_filename)
            else:
                print(f"File does not exist after download: {filename}")
        else:
            print(f"Download failed for {result}")
    except Exception as e:
        print(f"An error occurred with {result}: {e}")
