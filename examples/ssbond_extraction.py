from Bio.PDB import PDBParser, PDBIO, Select
import os


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

    def accept_residue(self, residue):
        residue_id = (residue.get_parent().id, residue.id[1])
        if residue_id in self.disulfide_residues and residue.resname == 'CYS':
            atom_names = {atom.name for atom in residue if atom.altloc in [' ', 'A']}
            return all(req_atom in atom_names for req_atom in self.required_atoms['CYS'])
        return False


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
pdb_directory = './ssbond_practice'
output_directory = './ssbond_practice_output'
if not os.path.exists(output_directory):
    os.makedirs(output_directory)
process_pdb_files(pdb_directory, output_directory)