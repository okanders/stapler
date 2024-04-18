import numpy as np
import itertools
import glob
import tqdm
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem
import xbin  # Assuming xbin is a hypothetical module you have access to
import getpy as gp  # For the multidimensional hash table
from homog import hstub
import torch
from torch_geometric.data import Data
import pkg_resources




# Assumed definition for atom type one-hot encoding
def one_hot_encode(atom_type):
    # Define atom types present in cysteine residues
    types = ['N', 'CA', 'C', 'O', 'CB', 'SG']  # Common atoms in cysteine
    return [1 if atom_type == t else 0 for t in types]

def retrieve_chi_angles_from_pdb(mol, hash_function, hash_table):
    # Verify molecule validity
    if mol is None or mol.GetNumAtoms() == 0:
        return None, "Invalid or empty molecule"
    if mol.GetNumConformers() == 0:
        return None, "No conformer found"

    # Use the first conformer
    conf = mol.GetConformer(0)
    chi_angles = []

    # Define atom indices for chi angle calculations based on your molecule
    # These indices should match those used in your hash table construction
    atom_indices = [0, 1, 2, 4, 5, 11, 10, 7, 6, 8]
    xyz = np.array([list(conf.GetAtomPosition(idx)) for idx in atom_indices])

    # Check for NaN values in coordinates, which indicate missing data
    if np.isnan(xyz).any():
        return None, "NaN values found in atom positions"

    try:
        # Generate transformation stubs based on specific atom positions
        # Generate transformation stubs based on atom positions
        stubs_a = hstub(xyz[:3], xyz[1:4], xyz[2:5])  # Stub for the first cysteine
        stubs_b = hstub(xyz[6:9], xyz[7:], xyz[8:])  # Stub for the second cysteine

        # Calculate transformations
        xforms_ab = np.linalg.inv(stubs_a) @ stubs_b
        xforms_ba = np.linalg.inv(stubs_b) @ stubs_a


        print("Xforms AB:", xforms_ab)
        print("Xforms BA:", xforms_ba)



        # Retrieve corresponding bin indices
        keys_ab = hash_function.get_bin_index(xforms_ab)
        keys_ba = hash_function.get_bin_index(xforms_ba)

        print("Hash keys AB:", keys_ab)
        print("Hash keys BA:", keys_ba)


        # Fetch chi angles for the given keys from the hash table
        chi_abc_values = [hash_table.get(key, None) for key in keys_ab]
        chi_bac_values = [hash_table.get(key, None) for key in keys_ba]



        # Print statements to verify the extracted data
        print("XYZ coordinates:", xyz)
        print("Stubs A:", stubs_a)
        print("Stubs B:", stubs_b)
        print("Chi angles AB:", chi_abc_values)
        print("Chi angles BA:", chi_bac_values)

        # Verify that chi angles were found
        if chi_abc_values is None or chi_bac_values is None:
            return None, "Chi angles not found for given keys"

        # Append retrieved chi angles to the list
        chi_angles.append((chi_abc_values, chi_bac_values))
    except Exception as e:
        # Handle any errors during processing
        return None, f"Error in processing: {str(e)}"

    return chi_angles, None



def preprocess_and_retrieve_chi_angles(directory):
    pdb_files = glob.glob(f'{directory}/*.pdb')
    dataset = []

    # Initialize hash function and hash table
    hash_function_kwargs = {'cart_resl': 1.0, 'ori_resl': 15.0, 'cart_bound': 512.0}
    hash_table_kwargs = {
        'key_type': np.dtype('i8'), 
        'value_type': np.dtype('i4'), 
        'filename': pkg_resources.resource_filename('stapler', 'hash_tables/default_native_disulfide_helicon_base_case.bin')
    }
    hash_function = xbin.XformBinner(**hash_function_kwargs)
    hash_table = gp.MultiDict(**hash_table_kwargs)
    hash_table.load(pkg_resources.resource_filename('stapler', 'hash_tables/default_native_disulfide_helicon_base_case.bin'))

    for pdb_file in tqdm.tqdm(pdb_files):
        try:
            with open(pdb_file, 'r') as file:
                pdb_content = file.read()
            mol = Chem.MolFromPDBBlock(pdb_content, removeHs=False)
            if mol is None or mol.GetNumAtoms() == 0 or mol.GetNumConformers() == 0:
                raise ValueError(f"Invalid or empty molecule for {pdb_file}")
            AllChem.AssignAtomChiralTagsFromStructure(mol)
            AllChem.Compute2DCoords(mol)
        except Exception as e:
            print(f"Error processing {pdb_file}: {str(e)}")
            continue  # Skip this file due to conversion error

        nodes, pos, edge_index = [], [], []
        for atom in mol.GetAtoms():
            nodes.append(one_hot_encode(atom.GetSymbol()))
            position = mol.GetConformer().GetAtomPosition(atom.GetIdx())
            pos.append([position.x, position.y, position.z])

        for bond in mol.GetBonds():
            edge_index.append([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
            edge_index.append([bond.GetEndAtomIdx(), bond.GetBeginAtomIdx()])  # Ensure graph is undirected

        graph_data = Data(x=torch.tensor(nodes, dtype=torch.float),
                          edge_index=torch.tensor(edge_index, dtype=torch.long).t().contiguous(),
                          pos=torch.tensor(pos, dtype=torch.float))
                          
        chi_angles, error = retrieve_chi_angles_from_pdb(mol, hash_function, hash_table)
        if error:
            print(f"Error processing {pdb_file}: {error}")
            continue  # Skip this file due to chi angle retrieval error
        else:
            graph_data.chi_angles = chi_angles
            dataset.append(graph_data)

        # Optional: Print information about each molecule after processing
        print(f"Finished processing {pdb_file}. Node count: {len(nodes)}, Edge count: {len(edge_index)//2}")

    return dataset


# Example usage
directory = '../examples/helicon_base_case_bin'
dataset = preprocess_and_retrieve_chi_angles(directory)
print(dataset) 

##VERSION 2

# # Assumed definition for atom type one-hot encoding
# def one_hot_encode(atom_type):
#     # Define atom types present in cysteine residues
#     types = ['N', 'CA', 'C', 'O', 'CB', 'SG']  # Common atoms in cysteine
#     return [1 if atom_type == t else 0 for t in types]

# # Main processing function
# def preprocess_and_retrieve_chi_angles(directory):
#     """
#     Processes PDB files to create graph structures and retrieves corresponding chi angles.
    
#     Args:
#         directory: Path to the directory containing PDB files.

#     Returns:
#         A list of Data objects, each augmented with chi angle information.
#     """
#     pdb_files = glob.glob(f'{directory}/*.pdb')
#     dataset = []

#     # Hash function and hash table settings
#     hash_function_kwargs = {'cart_resl': 1.0, 'ori_resl': 15.0, 'cart_bound': 512.0}
#     hash_table_kwargs = {
#         'key_type': np.dtype('i8'), 
#         'value_type': np.dtype('i4'), 
#         'filename': pkg_resources.resource_filename('stapler', 'hash_tables/default_native_disulfide_helicon_base_case.bin')
#     }
    
#     hash_function = xbin.XformBinner(**hash_function_kwargs)
#     hash_table = gp.MultiDict(**hash_table_kwargs)

#     for pdb_file in tqdm.tqdm(pdb_files):
#         with open(pdb_file, 'r') as file:
#             pdb_content = file.read()
#         mol = Chem.MolFromPDBBlock(pdb_content, removeHs=False)
#         if mol is None or mol.GetNumAtoms() == 0:
#             continue  # Skip invalid or empty molecules

#         AllChem.AssignAtomChiralTagsFromStructure(mol)
#         AllChem.Compute2DCoords(mol)

#         nodes, pos, edge_index = [], [], []
#         for atom in mol.GetAtoms():
#             nodes.append(one_hot_encode(atom.GetSymbol()))
#             position = mol.GetConformer().GetAtomPosition(atom.GetIdx())
#             pos.append([position.x, position.y, position.z])
        
#         for bond in mol.GetBonds():
#             edge_index.append([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
#             edge_index.append([bond.GetEndAtomIdx(), bond.GetBeginAtomIdx()])  # Undirected graph

#         graph_data = Data(x=torch.tensor(nodes, dtype=torch.float),
#                           edge_index=torch.tensor(edge_index, dtype=torch.long).t().contiguous(),
#                           pos=torch.tensor(pos, dtype=torch.float))
#         dataset.append(graph_data)

#         # Retrieve chi angles using previously defined logic
#         # This might need to be adapted to your specific data format and extraction logic
        
        
    
        
#         # Process XYZs similarly as in the Stapler class; this might require specific methods from your pose object
#         # Here, I assume XYZ processing similar to the Stapler's apply method; adapt as needed
#         xyzs_a = np.array([[residue.atom(atom).xyz() for atom in self.atom_selectors[0]] for residue in pose.residues if not residue.is_virtual_residue()])
#         xyzs_b = np.array([[residue.atom(atom).xyz() for atom in self.atom_selectors[0]] for residue in pose.residues if not residue.is_virtual_residue()])

        
#         # You might need to adapt the selection logic based on actual methods available in your pose object
#         sele = np.array([
#             [i, j] for i, j in itertools.product(range(len(xyzs_a)), range(len(xyzs_b)))
#             if np.abs(i - j) >= stapler.minimum_sequence_distance and np.linalg.norm(xyzs_a[i] - xyzs_b[j]) <= stapler.maximum_neighborhood_distance
#         ])

#         sele_a = np.squeeze(sele[:,0])
#         sele_b = np.squeeze(sele[:,1])        
        
#         chi_angles = []
#         for i, j in sele:
#             # Compute transformation matrices (stubs) based on XYZs, then apply them
#             stubs_a = hstub(xyzs_a[sele_a,0,:], xyzs_a[sele_a,1,:], xyzs_a[sele_a,2,:])
#             stubs_b = hstub(xyzs_b[sele_b,0,:], xyzs_b[sele_b,1,:], xyzs_b[sele_b,2,:])

#             xforms_ab = np.linalg.inv(stubs_a) @ stubs_b
#             xforms_ba = np.linalg.inv(stubs_b) @ stubs_a

#             # Hash transformation to get bin index
#             bin_index_ab = hash_function.get_bin_index(xforms_ab)

#             bin_index_ba = hash_function.get_bin_index(xforms_ba)

#             # Retrieve chi angles from hash table
#             chi_ab = hash_table.get(bin_index_ab, None)
#             chi_ba = hash_table.get(bin_index_ba, None)
            
#             # Convert binary data to chi angles if needed
#             chi_angle_touple = (chi_ab, chi_ba)
#             chi_angles.append(((i, j), chi_angle_touple))


#         graph_data.chi_angles = chi_angles  # Attach chi angles to graph data
#         dataset.append(graph_data)

#     return dataset

# # Example usage
# directory = '../examples/helicon_base_case_bin'
# dataset = preprocess_and_retrieve_chi_angles(directory)
# print(dataset) 







#Version 1





# def preprocess_structure(pdb_file):
#     mol = Chem.MolFromPDBFile(pdb_file)
#     if mol is None:
#         return None

#     nodes = []
#     edge_index = []
#     pos = []  # Positions for 3D visualization or geometric features

#     # Assuming you have a way to iterate over atoms and bonds in your molecule
#     for atom in mol.GetAtoms():
#         # Create node features, for example, one-hot encoding of atom types
#         nodes.append(one_hot_encode(atom.GetSymbol()))
#         pos.append(atom.GetPosition())  # Assuming 3D coordinates available

#     for bond in mol.GetBonds():
#         # Create edges based on bonds, could also include spatial distances
#         edge_index.append([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])

#     # Convert to tensors or suitable format for your GNN framework
#     x = torch.tensor(nodes, dtype=torch.float)
#     edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
#     pos = torch.tensor(pos, dtype=torch.float)

#     # No chi angle information included here; this is for structural data only
#     return Data(x=x, edge_index=edge_index, pos=pos)



# def retrieve_chi_angles(pose):
#     """
#     Retrieves chi angles for valid residue pairs in a given pose based on the stapler's configuration.
    
#     Args:
#         pose: The molecular pose (structure).

#     Returns:
#         A list of tuples containing residue pair indices and their corresponding chi angles, which are in tuples.
#     """
#     # Define your hash function and hash table parameters
#     hash_function_kwargs = {'cart_resl': 1.0, 'ori_resl': 15.0, 'cart_bound': 512.0}
#     hash_table_kwargs = {
#         'key_type': np.dtype('i8'), 
#         'value_type': np.dtype('i4'), 
#         'filename': pkg_resources.resource_filename('stapler', 'hash_tables/default_native_disulfide_helicon_base_case.bin')
#     }

#     # Initialize your hash function and hash table
#     hash_function = xbin.XformBinner(**hash_function_kwargs)
#     hash_table = gp.MultiDict(**hash_table_kwargs)
    
#     # Process XYZs similarly as in the Stapler class; this might require specific methods from your pose object
#     # Here, I assume XYZ processing similar to the Stapler's apply method; adapt as needed
#     xyzs_a = np.array([[residue.atom(atom).xyz() for atom in self.atom_selectors[0]] for residue in pose.residues if not residue.is_virtual_residue()])
#     xyzs_b = np.array([[residue.atom(atom).xyz() for atom in self.atom_selectors[0]] for residue in pose.residues if not residue.is_virtual_residue()])

    
#     # You might need to adapt the selection logic based on actual methods available in your pose object
#     sele = np.array([
#         [i, j] for i, j in itertools.product(range(len(xyzs_a)), range(len(xyzs_b)))
#         if np.abs(i - j) >= stapler.minimum_sequence_distance and np.linalg.norm(xyzs_a[i] - xyzs_b[j]) <= stapler.maximum_neighborhood_distance
#     ])

#     sele_a = np.squeeze(sele[:,0])
#     sele_b = np.squeeze(sele[:,1])


#     # List to store the chi angles
#     chi_angles = []

#     # Iterate over selected pairs and retrieve chi angles
#     for i, j in sele:
#         # Compute transformation matrices (stubs) based on XYZs, then apply them
#         stubs_a = hstub(xyzs_a[sele_a,0,:], xyzs_a[sele_a,1,:], xyzs_a[sele_a,2,:])
#         stubs_b = hstub(xyzs_b[sele_b,0,:], xyzs_b[sele_b,1,:], xyzs_b[sele_b,2,:])

#         xforms_ab = np.linalg.inv(stubs_a) @ stubs_b
#         xforms_ba = np.linalg.inv(stubs_b) @ stubs_a

#         # Hash transformation to get bin index
#         bin_index_ab = hash_function.get_bin_index(xforms_ab)

#         bin_index_ba = hash_function.get_bin_index(xforms_ba)

#         # Retrieve chi angles from hash table
#         chi_ab = hash_table.get(bin_index_ab, None)
#         chi_ba = hash_table.get(bin_index_ba, None)
        
#         if chis is not None:
#             # Convert binary data to chi angles if needed
#             chi_angle = ...  # Convert `chis` to actual angles based on your logic
#             chi_angles.append(((i, j), (chi_ab, chi_ba)))

#     return chi_angles





