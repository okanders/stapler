import glob
import os

import pyrosetta
pyrosetta.init()

from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from pyrosetta.rosetta.core.select.residue_selector import SecondaryStructureSelector
from pyrosetta.rosetta.core.select.residue_selector import AndResidueSelector
from pyrosetta.rosetta.core.select.residue_selector import TrueResidueSelector

from stapler import NativeDisulfideStapler

# Initialize residue selectors
chain_selector = ChainSelector('B')  # Select residues from chain B
helix_selector = SecondaryStructureSelector('H')  # Select residues in helices within chain B

# Combine chain and helix selectors to select residues within the helix of chain B
intra_chain_helix_selector = AndResidueSelector(chain_selector, helix_selector)

default_residue_selectors = [TrueResidueSelector(), TrueResidueSelector()]

only_binder_residue_selectors = [ChainSelector('B'), ChainSelector('B')]

intra_chain_param = [intra_chain_helix_selector, intra_chain_helix_selector]

# Initialize the native disulfide stapler with the intra-chain helix selector
native_disulfide_stapler = NativeDisulfideStapler(
    residue_selectors=default_residue_selectors,
    minimum_sequence_distance=1
)

# Process each input PDB file
for pdb in glob.glob('_inputs/*.pdb'):
    pdb_filename = os.path.splitext(os.path.basename(pdb))[0]
    pose = pyrosetta.pose_from_file(pdb)

    # Counter for numbering crosslinked poses
    crosslink_counter = 0

    # Attempt to apply stapling, handle empty selections gracefully
    try:
        for crosslinked_pose in native_disulfide_stapler.apply(pose):
            crosslinked_pose.dump_pdb(f'_practice_outputs_helix/{pdb_filename}_staple_{crosslink_counter}.pdb')
            crosslink_counter += 1

        if crosslink_counter == 0:
            print(f"No disulfide staples were applied for {pdb_filename}.")
    except IndexError as e:
        print(f"An error occurred while processing {pdb_filename}: {str(e)}")
        print("This may be due to no valid residue pairs found for stapling.")






# import glob
# import os

# import pyrosetta
# pyrosetta.init()

# from pyrosetta.rosetta.core.select.residue_selector import TrueResidueSelector
# from pyrosetta.rosetta.core.select.residue_selector import NotResidueSelector
# from pyrosetta.rosetta.core.select.residue_selector import AndResidueSelector
# from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
# from pyrosetta.rosetta.core.select.residue_selector import SecondaryStructureSelector

# from stapler import NativeDisulfideStapler

# # Initialize selectors
# chain_selector = ChainSelector('B')
# ss_selector = SecondaryStructureSelector('HE')
# not_ss_selector = NotResidueSelector(ss_selector)

# intra_chain_non_loop_selector = AndResidueSelector(chain_selector, not_ss_selector)

# native_disulfide_stapler = NativeDisulfideStapler(
#     residue_selectors=[intra_chain_non_loop_selector, intra_chain_non_loop_selector],
#     minimum_sequence_distance=1
# )

# for pdb in glob.glob('_inputs/*.pdb'):
#     pdb_filename = os.path.splitext(os.path.basename(pdb))[0]
#     pose = pyrosetta.pose_from_file(pdb)

#     # Counter for numbering crosslinked poses
#     crosslink_counter = 0

#     # Attempt to apply stapling, handle empty selections gracefully
#     try:
#         for crosslinked_pose in native_disulfide_stapler.apply(pose):
#             crosslinked_pose.dump_pdb(f'_outputs/{pdb_filename}_staple_{crosslink_counter}.pdb')
#             crosslink_counter += 1

#         if crosslink_counter == 0:
#             print(f"No disulfide staples were applied for {pdb_filename}.")
#     except IndexError as e:
#         print(f"An error occurred while processing {pdb_filename}: {str(e)}")
#         print("This may be due to no valid residue pairs found for stapling.")


##OG


# import glob
# import os

# import pyrosetta
# pyrosetta.init()

# from pyrosetta.rosetta.core.select.residue_selector import TrueResidueSelector
# from pyrosetta.rosetta.core.select.residue_selector import NotResidueSelector
# from pyrosetta.rosetta.core.select.residue_selector import AndResidueSelector
# from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
# from pyrosetta.rosetta.core.select.residue_selector import SecondaryStructureSelector

# from stapler import NativeDisulfideStapler


# # Initialize selectors
# chain_selector = ChainSelector('B')  # Select only chain B
# ss_selector = SecondaryStructureSelector('HE')  # Select alpha helices and beta sheets, if desired
# not_ss_selector = NotResidueSelector(ss_selector)  # Select not in alpha helices and beta sheets, if desired

# # Combining selectors for specific criteria (adjust according to your needs)
# intra_chain_non_loop_selector = AndResidueSelector(chain_selector, not_ss_selector)

# # Initialize the native disulfide stapler with intra-chain and structural preferences.
# native_disulfide_stapler = NativeDisulfideStapler(
#     residue_selectors=[intra_chain_non_loop_selector, intra_chain_non_loop_selector],
#     minimum_sequence_distance=4  # Adjust this as necessary
# )

# for pdb in glob.glob('_inputs/*.pdb'):
#     pdb_filename = os.path.splitext(os.path.basename(pdb))[0]

#     pose = pyrosetta.pose_from_file(pdb)

#     # Apply the stapling and save each resulting pose
#     for i, crosslinked_pose in enumerate(native_disulfide_stapler.apply(pose)):
#         crosslinked_pose.dump_pdb(f'_outputs/{pdb_filename}_staple_{i}.pdb')
