import os
import sqlite3

# Attempt to query Cys.sqlite to query native disulfide bonds in proteins

# Path to your Cys.sqlite database
db_path = './jcim_cys.sqlite'
if not os.path.exists(db_path):
    print(f"Database not found at {db_path}")
    exit(1)

# Connect to the SQLite database
conn = sqlite3.connect(db_path)
cursor = conn.cursor()

# Create a directory to store PDB files
output_dir = './native_disulfide_pdb'
os.makedirs(output_dir, exist_ok=True)

# Adjusted SQL Query to fetch necessary coordinates, angles, and distances for constructing disulfide bridges
query = """
SELECT cc.pdb_id, 
       conf1.N_x, conf1.N_y, conf1.N_z,
       conf1.CA_x, conf1.CA_y, conf1.CA_z,
       conf1.C_x, conf1.C_y, conf1.C_z,
       conf1.SG_x, conf1.SG_y, conf1.SG_z,
       conf2.SG_x as SG2_x, conf2.SG_y as SG2_y, conf2.SG_z as SG2_z,
       cc.SGi_SGj, cc.CAi_CAj, cc.CBi_SGi_SGj, cc.SGi_SGj_CBj,
       cc.CAi_CBi_SGi_SGj, cc.CBi_SGi_SGj_CBj, cc.SGi_SGj_CBj_CAj
FROM Cys_Cys cc
JOIN Cys_Conf conf1 ON cc.cys_conf_idi = conf1.id
JOIN Cys_Conf conf2 ON cc.cys_conf_idj = conf2.id
WHERE conf1.pdb_id = conf2.pdb_id
AND conf1.model_num = conf2.model_num
AND conf1.alt_id = conf2.alt_id
AND cc.SGi_SGj < 2.05  -- This condition ensures that the distance between the two sulfur atoms is within a reasonable range for a disulfide bridge
LIMIT 30000;
"""

try:
    cursor.execute(query)
    entries = cursor.fetchall()
    print(f"Found {len(entries)} entries to create PDB files for.")
    
    # Generate PDB files based on fetched coordinates and disulfide bridge information
    for entry in entries:
        pdb_id, *data = entry  # Unpack all fetched coordinates and disulfide bridge data
        pdb_file_path = os.path.join(output_dir, f"{pdb_id}.pdb")
        with open(pdb_file_path, 'w') as pdb_file:
            pdb_file.write("REMARK   Generated from Cys.sqlite\n")
            atom_names = ["N", "CA", "C", "CB", "SG1", "SG2"]
            for i, (x, y, z) in enumerate(zip(data[::3], data[1::3], data[2::3])):
                if i < 5:  # SG1 coordinates
                    atom_type = 'S' if i == 4 else 'C'
                    pdb_file.write(f"ATOM  {i+1:5d}  {atom_names[i]:<4} CYS A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           {atom_type}\n")
                else:  # SG2 coordinates
                    pdb_file.write(f"ATOM  {i+1:5d}  {atom_names[i]:<4} CYS B   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           S\n")
            pdb_file.write("SSBOND   1 CYS A    1    CYS B    1\n")
            pdb_file.write("TER\n")
        print(f"Created PDB file: {pdb_id}.pdb")
except Exception as e:
    print(f"An error occurred: {e}")
finally:
    cursor.close()
    conn.close()






### BEST So far


# import os
# import sqlite3

# # Path to your Cys.sqlite database
# db_path = './jcim_cys.sqlite'
# if not os.path.exists(db_path):
#     print(f"Database not found at {db_path}")
#     exit(1)

# # Connect to the SQLite database
# conn = sqlite3.connect(db_path)
# cursor = conn.cursor()

# # Create a directory to store PDB files
# output_dir = './native_disulfide_pdb'
# os.makedirs(output_dir, exist_ok=True)

# # Adjusted SQL Query to fetch necessary coordinates and angles
# # Adjusted SQL Query to fetch necessary coordinates and angles
# query = """
# SELECT conf.pdb_id, 
#        conf.N_x, conf.N_y, conf.N_z,
#        conf.CA_x, conf.CA_y, conf.CA_z,
#        conf.C_x, conf.C_y, conf.C_z,
#        conf.SG_x, conf.SG_y, conf.SG_z  
# FROM Cys_Conf conf
# JOIN PDB pdb ON conf.pdb_id = pdb.id
# WHERE pdb.exp_method = 'X-RAY DIFFRACTION' AND pdb.resolution < 2.0
# LIMIT 30000;
# """


# try:
#     cursor.execute(query)
#     entries = cursor.fetchall()
#     print(f"Found {len(entries)} entries to create PDB files for.")
    
#     # Generate PDB files based on fetched coordinates
#     for entry in entries:
#         pdb_id, *coords = entry  # Unpack all fetched coordinates
#         pdb_file_path = os.path.join(output_dir, f"{pdb_id}.pdb")
#         with open(pdb_file_path, 'w') as pdb_file:
#             pdb_file.write("REMARK   Generated from Cys.sqlite\n")
#             atom_names = ["N", "CA", "C", "CB", "SG"]
#             for i, (x, y, z) in enumerate(zip(coords[::3], coords[1::3], coords[2::3])):
#                 pdb_file.write(f"ATOM  {i+1:5d}  {atom_names[i]:<4} CYS A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n")
#             pdb_file.write("TER\n")
#         print(f"Created PDB file: {pdb_id}.pdb")
# except Exception as e:
#     print(f"An error occurred: {e}")
# finally:
#     cursor.close()
#     conn.close()









#### SIMPLE



# import os
# import sqlite3
# from urllib.request import urlretrieve

# # Path to your Cys.sqlite database
# db_path = './jcim_cys.sqlite'
# if not os.path.exists(db_path):
#     print(f"Database not found at {db_path}")
#     exit(1)

# # Connect to the SQLite database
# conn = sqlite3.connect(db_path)
# cursor = conn.cursor()

# # Create a directory to store PDB files
# output_dir = './native_disulfide_pdbs'
# os.makedirs(output_dir, exist_ok=True)

# # SQL query to fetch PDB IDs with resolution under 2.0 and disulfide bonds
# query = """
# SELECT DISTINCT pdb.id
# FROM PDB pdb
# JOIN Cys_Cys cc ON pdb.id = cc.pdb_id
# WHERE pdb.exp_method = 'X-RAY DIFFRACTION'
# AND pdb.resolution <= 2.0
# LIMIT 32000;
# """

# try:
#     cursor.execute(query)
#     pdb_ids = cursor.fetchall()
#     print(f"Found {len(pdb_ids)} entries.")
    
#     # Download PDB files and save them in the output directory
#     for pdb_id in pdb_ids:
#         pdb_id = pdb_id[0]  # Unpack the tuple
#         pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
#         pdb_file = os.path.join(output_dir, f"{pdb_id}.pdb")
#         urlretrieve(pdb_url, pdb_file)
#         print(f"Downloaded PDB file: {pdb_id}.pdb")
# except Exception as e:
#     print(f"An error occurred: {e}")
# finally:
#     conn.close()
