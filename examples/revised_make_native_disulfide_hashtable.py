import glob
from multiprocessing import Pool

import tqdm
import numpy as np
import getpy as gp
from rdkit import Chem

import nerf
import xbin
import homog

# Constants and configurations for hashing
CHI_RESOLUTION = 30.0
assert 360.0 / CHI_RESOLUTION < 256.0, "CHI_RESOLUTION must be less than 12."

hash_function_kwargs = {'cart_resl': 1.0, 'ori_resl': 15.0, 'cart_bound': 512.0}
hash_function = xbin.XformBinner(**hash_function_kwargs)

hash_table_kwargs = {'key_type': np.dtype('i8'), 'value_type': np.dtype('i4')}
hash_table = gp.MultiDict(**hash_table_kwargs)

def process_pdb(pdb_file):
    mol = Chem.MolFromPDBFile(pdb_file)
    if not mol or mol.GetNumAtoms() == 0:
        return f"Invalid or empty molecule found in {pdb_file}, skipping."
    
    if mol.GetNumConformers() == 0:
        return f"No conformer found for {pdb_file}, skipping."

    # Get the first conformer for processing
    conf = mol.GetConformer(0)
    # Extract required atom positions from the conformer
    atom_indices = [0, 1, 2, 4, 5, 11, 10, 7, 6, 8]  # Specific atom indices used in processing
    atom_positions = [conf.GetAtomPosition(idx) for idx in atom_indices]

    # Prepare data for iNeRF and NeRF calculations
    xyz = np.array([(pos.x, pos.y, pos.z) for pos in atom_positions])
    deps = np.array([
        [0, 0, 0],  # Root atom, no dependencies
        [0, 0, 0],  # Another root atom, no dependencies
        [1, 0, 0],  # Dependent on first atom
        [1, 0, 2],  # Dependent on first and third atoms
        [3, 1, 0],  # Perturb SG1-CB1-CA1-N1
        [4, 3, 1],  # Perturb SG2-SG1-CB1-CA1
        [5, 4, 3],  # Perturb CB2-SG2-SG1-CB1
        [6, 5, 4],  # Perturb CA2-CB2-SG2-SG1
        [7, 6, 5],  # Perturb N2-CA2-CB2-SG2
        [7, 6, 8]   # Dependent on seventh, eighth, and tenth atoms
    ])
    # Internal coordinates manipulation using iNeRF and generating new conformations with NeRF

    dof = nerf.iNeRF(xyz, deps=deps)

    repeats = 100
    dofs = np.repeat(dof[np.newaxis], repeats, axis=0)
    for i in [4, 5, 6, 7, 8]: # [CA1-CB1, CB1-SG1, SG1-SG2, SG2-CB2, CB2-CA2]
        dofs[:,i,2] += np.random.uniform(-5, 5, repeats)
        dofs[:,i,2] = (dofs[:,i,2] + (180.0 + CHI_RESOLUTION/2.0)) % 360.0 - (180.0 + CHI_RESOLUTION/2.0)

    # Compute transformations and update hash table
    update_hash_table(dofs, deps, hash_function, hash_table)

    return f"Processed {pdb_file}, updated hash table."

def update_hash_table(dofs, deps, hash_function, hash_table):
    # Generate new conformations using NeRF based on the DOFs
    xyzs = nerf.NeRF(dofs, deps=deps)

    # Compute stubs for the transformation between conformations
    stubs_a = homog.hstub(xyzs[:, 0, :], xyzs[:, 1, :], xyzs[:, 2, :])  # N1 CA1 C1
    stubs_b = homog.hstub(xyzs[:, -2, :], xyzs[:, -3, :], xyzs[:, -1, :])  # N2 CA2 C2

    # Compute transformations from A to B and from B to A
    xforms_ab = np.linalg.inv(stubs_a) @ stubs_b
    xforms_ba = np.linalg.inv(stubs_b) @ stubs_a

    # Get bin indices for both transformations
    keys_ab = hash_function.get_bin_index(xforms_ab)
    keys_ba = hash_function.get_bin_index(xforms_ba)

    # Compute the chi angles for insertion into the hash table, adjusting for resolution
    chis_ab = np.around(dofs[:, [4, 8], 2] / CHI_RESOLUTION).astype('i2').flatten().view(np.dtype('i4'))
    chis_ba = np.around(dofs[:, [8, 4], 2] / CHI_RESOLUTION).astype('i2').flatten().view(np.dtype('i4'))

    # Update the hash table with the new keys and chi angles
    # This assumes your hash_table can handle bulk updates; if not, you may need a loop
    for k, v in zip(keys_ab, chis_ab):
        hash_table[k] = v
    for k, v in zip(keys_ba, chis_ba):
        hash_table[k] = v


def main():
    pdb_files = glob.glob('./disulfide_structures_BIG5/*.pdb')
    with Pool(processes=4) as pool:  # Adjust number of processes based on your machine's capability
        for result in tqdm.tqdm(pool.imap(process_pdb, pdb_files), total=len(pdb_files)):
            tqdm.tqdm.write(result)

    # Once all files are processed, dump the hash table
    hash_table.dump('../stapler/hash_tables/default_native_disulfide_BIG5_revised.bin')

if __name__ == "__main__":
    main()
