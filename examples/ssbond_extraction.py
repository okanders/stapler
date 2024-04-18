from Bio.PDB import PDBParser, PDBIO, Select
import os

##Best attempt at hierachy of conformations within residues
# class DisulfideSelect(Select):
#     def __init__(self, disulfide_residues):
#         self.disulfide_residues = disulfide_residues
#         self.required_atoms = {'CYS': ['N', 'CA', 'C', 'O', 'CB', 'SG']}
#         # Store best seen conformation for each atom in each residue
#         self.best_seen_conformations = {}

#     def accept_atom(self, atom):
#         parent_residue = atom.get_parent()
#         residue_id = (parent_residue.get_parent().id, parent_residue.id[1])
#         atom_key = (residue_id, atom.name)

#         # Define preference for conformations
#         conformation_preference = [' ', 'A', 'B', 'C', 'D']  # Add more if needed

#         print(f"Checking atom: {atom.name}, Residue: {residue_id}, Conformation: '{atom.altloc}'")  # Added for debugging

#         # Check if the atom is part of the disulfide residue and is one of the required atoms
#         if residue_id in self.disulfide_residues and atom.name in self.required_atoms['CYS']:
#             current_best = self.best_seen_conformations.get(atom_key, ('D',))[0]  # Set 'D' as default which is the least preferred

#             # Update the best conformation seen for the atom if the current atom's conformation is preferred
#             if atom.altloc in conformation_preference:
#                 if conformation_preference.index(atom.altloc) < conformation_preference.index(current_best):
#                     self.best_seen_conformations[atom_key] = (atom.altloc, atom)
#                     print(f"Accepting better conformation for atom: {atom.name}, Residue: {residue_id}, Conformation: '{atom.altloc}'")  # Added for debugging
#                     return True
#                 else:
#                     print(f"Rejecting due to better or same conformation already seen: {atom.name}, Residue: {residue_id}, Conformation: '{atom.altloc}'")  # Added for debugging
#                     return False
#             elif current_best == 'D':  # Accept if no better conformation has been recorded yet
#                 self.best_seen_conformations[atom_key] = (atom.altloc, atom)
#                 print(f"Accepting first seen conformation for atom: {atom.name}, Residue: {residue_id}, Conformation: '{atom.altloc}'")  # Added for debugging
#                 return True
#             else:
#                 print(f"Rejecting due to better conformation already present: {atom.name}, Residue: {residue_id}, Conformation: '{atom.altloc}'")  # Added for debugging
#                 return False
#         else:
#             print(f"Rejecting non-required or non-disulfide atom: {atom.name}, Residue: {residue_id}, Conformation: '{atom.altloc}'")  # Added for debugging
#             return False

#     def accept_residue(self, residue):
#         residue_id = (residue.get_parent().id, residue.id[1])
#         print(f"Checking residue: {residue_id}")  # Already exists but confirm it's still there

#         if residue_id in self.disulfide_residues and residue.resname == 'CYS':
#             required_atoms_check = {atom for atom in self.required_atoms['CYS']}
#             actual_atoms = {atom.name: atom for atom in residue.get_atoms()}
#             missing_atoms = required_atoms_check - actual_atoms.keys()
#             if missing_atoms:
#                 print(f"Missing atoms in residue {residue_id}: {missing_atoms}")
#                 return False
#             print(f"All required atoms present in residue {residue_id}")
#             return True
#         return False

#     def validate_disulfides(self, structure):
#         self.best_seen_conformations.clear()  # Reset at the start of validation
#         valid_residues = True

#         for model in structure:
#             for chain in model:
#                 for residue in chain:
#                     residue_id = (chain.id, residue.id[1])
#                     if residue_id in self.disulfide_residues:
#                         # Explicitly populating best_seen_conformations for the residue
#                         for atom in residue:
#                             self.accept_atom(atom)  # This should fill in the best_seen_conformations
#                         # After populating, we check for any missing required atoms
#                         missing_atoms = [atom for atom in self.required_atoms['CYS'] if (residue_id, atom) not in self.best_seen_conformations]
#                         if missing_atoms:
#                             print(f"Missing atoms in final check for residue {residue_id}: {missing_atoms}")
#                             valid_residues = False
#                             break
#                         else:
#                             print(f"Residue {residue_id} passed final check.")
#                 if not valid_residues:
#                     break
#             if not valid_residues:
#                 break
#         return valid_residues







# class DisulfideSelect(Select):
#     """Class to select disulfide-bonded cysteines with all required atoms for writing to a new PDB file."""
#     def __init__(self, disulfide_residues):
#         self.disulfide_residues = disulfide_residues
#         self.required_atoms = {'CYS': ['N', 'CA', 'C', 'O', 'CB', 'SG']}
#         self.seen_atoms = set()  # Track seen atoms to handle alternate conformations
#         self.selected_atoms = {}  # Track selected atoms to ensure all required atoms are chosen

#     def accept_atom(self, atom):
#         parent_residue = atom.get_parent()
#         residue_id = (parent_residue.get_parent().id, parent_residue.id[1])
#         atom_id = (residue_id, atom.name)

#         # Check if the residue is part of a disulfide bond and the atom is required
#         if residue_id in self.disulfide_residues and atom.name in self.required_atoms['CYS']:
#             # Select only the first conformation of the atom, but ensure all required atoms are selected
#             if atom_id not in self.seen_atoms:
#                 self.seen_atoms.add(atom_id)
#                 if atom.altloc in [' ', 'A'] or not self.selected_atoms.get(atom_id, False):
#                     self.selected_atoms[atom_id] = True
#                     return True
#             return False  # Ignore additional conformations
#         return False

#     def validate_disulfides(self, structure):
#         """Validate if selected disulfide residues have all the required atoms."""
#         # Reset selected atoms for new validation
#         self.selected_atoms = {}
#         valid_residues = True  # Assume all residues are valid initially
#         for residue in structure.get_residues():
#             residue_id = (residue.get_parent().id, residue.id[1])
#             if residue_id in self.disulfide_residues:
#                 missing_atoms = [atom for atom in self.required_atoms['CYS'] if not self.selected_atoms.get((residue_id, atom), False)]
#                 if missing_atoms:
#                     print(f"Missing atoms in residue {residue_id}: {missing_atoms}")
#                     valid_residues = False
#                     break  # Stop checking if any invalid residue is found
#         return valid_residues  # Return the final validity status





class DisulfideSelect(Select):
    """Class to select disulfide-bonded cysteines with all required atoms for writing to a new PDB file."""
    def __init__(self, disulfide_residues):
        self.disulfide_residues = disulfide_residues
        self.required_atoms = {'CYS': ['N', 'CA', 'C', 'O', 'CB', 'SG']}
        self.seen_atoms = set()  # Track seen atoms to handle alternate conformations

    def accept_atom(self, atom):
        parent_residue = atom.get_parent()
        residue_id = (parent_residue.get_parent().id, parent_residue.id[1])
        atom_key = (residue_id, atom.name)

        # Print statement for debugging
        print(f"Checking atom: {atom.name}, Residue: {residue_id}, Conformation: '{atom.altloc}'")

        # Check if the residue is part of a disulfide bond and the atom is required
        if residue_id in self.disulfide_residues and atom.name in self.required_atoms['CYS']:
            # Allow atoms with ' ' or 'A' conformation only, ensure each atom is evaluated only once
            if atom.altloc in [' ', 'A']:
                if (residue_id, atom.name, atom.altloc) not in self.seen_atoms:
                    self.seen_atoms.add((residue_id, atom.name, atom.altloc))
                    print(f"Accepting atom: {atom.name}, Residue: {residue_id}, Conformation: '{atom.altloc}'")
                    return True
                else:
                    print(f"Rejecting already seen atom: {atom.name}, Residue: {residue_id}, Conformation: '{atom.altloc}'")
                    return False
            else:
                print(f"Rejecting atom due to unsuitable conformation: {atom.name}, Residue: {residue_id}, Conformation: '{atom.altloc}'")
                return False
        else:
            print(f"Rejecting non-required atom or not part of disulfide bond: {atom.name}, Residue: {residue_id}")
            return False


#newer is first
    # def accept_residue(self, residue):
    #     residue_id = (residue.get_parent().id, residue.id[1])
    #     if residue_id in self.disulfide_residues and residue.resname == 'CYS':
    #         atom_names = {atom.name for atom in residue if atom.altloc in [' ', 'A']}
    #         missing_atoms = [req_atom for req_atom in self.required_atoms['CYS'] if req_atom not in atom_names]
    #         if missing_atoms:
    #             print(f"Missing atoms in residue {residue_id}: {missing_atoms}")
    #         return all(req_atom in atom_names for req_atom in self.required_atoms['CYS'])
    #     return False


    def accept_residue(self, residue):
        residue_id = (residue.get_parent().id, residue.id[1])
        if residue_id in self.disulfide_residues and residue.resname == 'CYS':
            atom_names = {atom.name for atom in residue if atom.altloc in [' ', 'A']}
            return all(req_atom in atom_names for req_atom in self.required_atoms['CYS'])
        return False


#newer is first
    # def validate_disulfides(self, structure):
    #     """Validate if selected disulfide residues have all the required atoms."""
    #     valid_residues = True  # Assume all residues are valid initially
    #     for residue in structure.get_residues():
    #         if (residue.get_parent().id, residue.id[1]) in self.disulfide_residues:
    #             if not self.accept_residue(residue):
    #                 print(f"Invalid or incomplete residue detected: {residue.get_parent().id}, {residue.id}")
    #                 valid_residues = False
    #                 break  # Stop checking if any invalid residue is found
    #     return valid_residues  # Return the final validity status


    def validate_disulfides(self, structure):
        """Validate if selected disulfide residues have all the required atoms."""
        valid_residues = True  # Assume all residues are valid initially
        for residue in structure.get_residues():
            if (residue.get_parent().id, residue.id[1]) in self.disulfide_residues:
                if not self.accept_residue(residue):
                    valid_residues = False
                    break  # Stop checking if any invalid residue is found
        return valid_residues  # Return the final validity status

def parse_ssbonds(pdb_path):
    ssbonds = []
    with open(pdb_path, 'r') as file:
        for line in file:
            if line.startswith('SSBOND'):
                chain_id1, res_num1, chain_id2, res_num2 = line[15], int(line[17:21].strip()), line[29], int(line[31:35].strip())
                # Skip adding the bond if it is an intra-disulfide bond (same residue)
                if res_num1 != res_num2:
                    ssbonds.append((chain_id1, res_num1, chain_id2, res_num2))
    return ssbonds

def process_pdb_files(directory, output_directory):
    parser = PDBParser(PERMISSIVE=1)
    io = PDBIO()
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    for filename in os.listdir(directory):
        if filename.endswith(".pdb"):
            path = os.path.join(directory, filename)
            print(f"Processing file: {path}")

            ssbonds = parse_ssbonds(path)
            print(f"Found SSBOND annotations: {ssbonds}")

            structure = parser.get_structure('PDB_structure', path)
            
            for index, ssbond in enumerate(ssbonds):
                chain_id1, res_num1, chain_id2, res_num2 = ssbond
                disulfide_residues = {(chain_id1, res_num1), (chain_id2, res_num2)}
                print(f"Adding disulfide bond: {(chain_id1, res_num1)} - {(chain_id2, res_num2)}")

                selector = DisulfideSelect(disulfide_residues)
                # # Only for hierarchical conformation. Before writing, print out best conformations determined
                # print("Best seen conformations before writing:", selector.best_seen_conformations)

                if selector.validate_disulfides(structure):  # Validate disulfide residues
                    output_path = os.path.join(output_directory, f"disulfide_{index + 1}_{filename}")
                    print(f"Writing disulfide structure to: {output_path}")
                    io.set_structure(structure)
                    io.save(output_path, selector)
                else:
                    print(f"Disulfide bond validation failed for {filename}.")

# Specify the directories for input and output PDB files
# pdb_directory = './native_disulfide_pdbs_BIG'
# output_directory = './disulfide_structures_BIG5'
pdb_directory = './ssbond_practice'
output_directory = './ssbond_practice_output'
if not os.path.exists(output_directory):
    os.makedirs(output_directory)
process_pdb_files(pdb_directory, output_directory)





###BEST I HAVE NO ALTERNATIVE CONFORMATIONS

# from Bio.PDB import PDBParser, PDBIO, Select
# import os


# class DisulfideSelect(Select):
#     """Class to select disulfide-bonded cysteines with all required atoms for writing to a new PDB file."""
#     def __init__(self, disulfide_residues):
#         self.disulfide_residues = disulfide_residues
#         self.required_atoms = {'CYS': ['N', 'CA', 'C', 'O', 'CB', 'SG']}

#     def accept_atom(self, atom):
#         parent_residue = atom.get_parent()
#         residue_id = (parent_residue.get_parent().id, parent_residue.id[1])
#         return residue_id in self.disulfide_residues and atom.name in self.required_atoms['CYS']

#     def accept_residue(self, residue):
#         residue_id = (residue.get_parent().id, residue.id[1])
#         if residue_id in self.disulfide_residues and residue.resname == 'CYS':
#             atom_names = {atom.name for atom in residue}
#             return all(req_atom in atom_names for req_atom in self.required_atoms['CYS'])
#         return False

#     def validate_disulfides(self, structure):
#         """Validate if selected disulfide residues have all the required atoms."""
#         valid_residues = True  # Assume all residues are valid initially
#         for residue in structure.get_residues():
#             if (residue.get_parent().id, residue.id[1]) in self.disulfide_residues:
#                 if not self.accept_residue(residue):
#                     valid_residues = False
#                     break  # Stop checking if any invalid residue is found
#         return valid_residues  # Return the final validity status

# def parse_ssbonds(pdb_path):
#     ssbonds = []
#     with open(pdb_path, 'r') as file:
#         for line in file:
#             if line.startswith('SSBOND'):
#                 chain_id1, res_num1, chain_id2, res_num2 = line[15], int(line[17:21].strip()), line[29], int(line[31:35].strip())
#                 # Skip adding the bond if it is an intra-disulfide bond (same residue)
#                 if res_num1 != res_num2:
#                     ssbonds.append((chain_id1, res_num1, chain_id2, res_num2))
#     return ssbonds

# def process_pdb_files(directory, output_directory):
#     parser = PDBParser(PERMISSIVE=1)
#     io = PDBIO()
#     if not os.path.exists(output_directory):
#         os.makedirs(output_directory)

#     for filename in os.listdir(directory):
#         if filename.endswith(".pdb"):
#             path = os.path.join(directory, filename)
#             print(f"Processing file: {path}")

#             ssbonds = parse_ssbonds(path)
#             print(f"Found SSBOND annotations: {ssbonds}")

#             structure = parser.get_structure('PDB_structure', path)

#             for index, ssbond in enumerate(ssbonds):
#                 chain_id1, res_num1, chain_id2, res_num2 = ssbond
#                 disulfide_residues = {(chain_id1, res_num1), (chain_id2, res_num2)}
#                 print(f"Adding disulfide bond: {(chain_id1, res_num1)} - {(chain_id2, res_num2)}")

#                 selector = DisulfideSelect(disulfide_residues)
#                 if selector.validate_disulfides(structure):  # Validate disulfide residues
#                     output_path = os.path.join(output_directory, f"disulfide_{index + 1}_{filename}")
#                     print(f"Writing disulfide structure to: {output_path}")
#                     io.set_structure(structure)
#                     io.save(output_path, selector)
#                 else:
#                     print(f"Disulfide bond validation failed for {filename}.")

# # Specify the directories for input and output PDB files
# pdb_directory = './native_disulfide_pdbs_BIG'
# output_directory = './disulfide_structures_BIG3'
# process_pdb_files(pdb_directory, output_directory)

###BEST I HAVE


# from Bio.PDB import PDBParser, PDBIO, Select
# import os


# class DisulfideSelect(Select):
#     """Class to select disulfide-bonded cysteines for writing to a new PDB file."""
#     def __init__(self, disulfide_residues):
#         # disulfide_residues should be a set of tuples (chain_id, residue_number)
#         self.disulfide_residues = disulfide_residues

#     def accept_atom(self, atom):
#         # Override accept_atom to ensure only atoms from the cysteine residues are selected
#         parent_residue = atom.get_parent()
#         residue_id = (parent_residue.get_parent().id, parent_residue.id[1])
#         return residue_id in self.disulfide_residues and parent_residue.resname == 'CYS'

#     def accept_residue(self, residue):
#         # Only accept residues that are part of the disulfide bonds and are cysteines
#         residue_id = (residue.get_parent().id, residue.id[1])  # (chain ID, residue number)
#         return residue_id in self.disulfide_residues and residue.resname == 'CYS'


# def parse_ssbonds(pdb_path):
#     ssbonds = []
#     with open(pdb_path, 'r') as file:
#         for line in file:
#             if line.startswith('SSBOND'):
#                 chain_id1, res_num1, chain_id2, res_num2 = line[15], int(line[17:21].strip()), line[29], int(line[31:35].strip())  # Corrected for string slicing
#                 ssbonds.append((chain_id1, res_num1, chain_id2, res_num2))
#     return ssbonds

# def process_pdb_files(directory, output_directory):
#     parser = PDBParser(PERMISSIVE=1)
#     io = PDBIO()

#     # Create output directory if it doesn't exist yet
#     if not os.path.exists(output_directory):
#         os.makedirs(output_directory)

#     for filename in os.listdir(directory):
#         if filename.endswith(".pdb"):
#             path = os.path.join(directory, filename)
#             print(f"Processing file: {path}")
            
#             # Directly parse the SSBOND entries from the PDB file
#             ssbonds = parse_ssbonds(path)
#             print(f"Found SSBOND annotations: {ssbonds}")

#             # Use Bio.PDB to parse the structure
#             structure = parser.get_structure('PDB_structure', path)

#             for index, ssbond in enumerate(ssbonds):
#                 chain_id1, res_num1, chain_id2, res_num2 = ssbond
#                 disulfide_residues = {(chain_id1, res_num1), (chain_id2, res_num2)}
#                 print(f"Adding disulfide bond: {(chain_id1, res_num1)} - {(chain_id2, res_num2)}")

#                 # If disulfide bonds are present, write out new PDB files for each bond
#                 if disulfide_residues:
#                     output_path = os.path.join(output_directory, f"disulfide_{index + 1}_{filename}")
#                     print(f"Writing disulfide structure to: {output_path}")
#                     io.set_structure(structure)
#                     io.save(output_path, DisulfideSelect(disulfide_residues))
#                 else:
#                     print(f"No disulfide bonds found in {filename}.")

# # Specify the directories for input and output PDB files
# pdb_directory = './native_disulfide_pdbs_BIG'
# output_directory = './disulfide_structures_BIG'
# process_pdb_files(pdb_directory, output_directory)



##OLD



# from Bio.PDB import PDBParser, PDBIO, Select
# import os

# class DisulfideSelect(Select):
#     """Class to select disulfide-bonded cysteines for writing to a new PDB file."""
#     def __init__(self, disulfide_residues):
#         self.disulfide_residues = disulfide_residues

#     def accept_residue(self, residue):
#         # Only accept residues that are part of the disulfide bonds
#         residue_id = (residue.get_parent().id, residue.id[1])  # (chain ID, residue number)
#         return residue_id in self.disulfide_residues



# def parse_ssbonds(pdb_path):
#     ssbonds = []
#     with open(pdb_path, 'r') as file:
#         for line in file:
#             if line.startswith('SSBOND'):
#                 parts = line.split()
#                 chain_id1 = line[15].strip()  # Chain ID for the first cysteine
#                 res_num1 = int(line[17:21].strip())  # Residue number for the first cysteine
#                 chain_id2 = line[29].strip()  # Chain ID for the second cysteine
#                 res_num2 = int(line[31:35].strip())  # Residue number for the second cysteine
#                 ssbonds.append((chain_id1, res_num1, chain_id2, res_num2))
#     return ssbonds


# def process_pdb_files(directory, output_directory):
#     parser = PDBParser(PERMISSIVE=1)
#     io = PDBIO()

#     # Create output directory if it doesn't exist yet
#     if not os.path.exists(output_directory):
#         os.makedirs(output_directory)

#     for filename in os.listdir(directory):
#         if filename.endswith(".pdb"):
#             path = os.path.join(directory, filename)
#             print(f"Processing file: {path}")
            
#             # Directly parse the SSBOND entries from the PDB file
#             ssbonds = parse_ssbonds(path)
#             print(f"Found SSBOND annotations: {ssbonds}")

#             # Use Bio.PDB to parse the structure
#             structure = parser.get_structure('PDB_structure', path)

#             for index, ssbond in enumerate(ssbonds):
#                 chain_id1, res_num1, chain_id2, res_num2 = ssbond
#                 disulfide_residues = {(chain_id1, res_num1), (chain_id2, res_num2)}
#                 print(f"Adding disulfide bond: {(chain_id1, res_num1)} - {(chain_id2, res_num2)}")

#                 # If disulfide bonds are present, write out new PDB files for each bond
#                 if disulfide_residues:
#                     output_path = os.path.join(output_directory, f"disulfide_{index + 1}_{filename}")
#                     print(f"Writing disulfide structure to: {output_path}")
#                     io.set_structure(structure)
#                     io.save(output_path, DisulfideSelect(disulfide_residues))
#             if not ssbonds:
#                 print(f"No disulfide bonds found in {filename}.")

# # Specify the directories for input and output PDB files
# pdb_directory = './native_disulfide_pdb'
# output_directory = './disulfide_structures'
# process_pdb_files(pdb_directory, output_directory)
