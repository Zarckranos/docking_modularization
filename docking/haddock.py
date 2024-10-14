import os
import subprocess

from Bio.PDB import PDBParser, PDBIO
from Bio import SeqIO

# HADDOCK TUTORIAL:
# https://www.bonvinlab.org/education/HADDOCK24/HADDOCK24-local-tutorial/
class Haddock:
    # lista de proteinas a serem modificadas
    PROTEINS = ['HIE', 'HID', 'HIP', 'CYX']

    def __init__(self, pdb_files_path, antibody, antigen, ab_type):
        self.pdb_files_path = pdb_files_path

        self.ab = antibody
        self.ag = antigen
        self.ab_type = ab_type

        self.ab_name = antibody.split('/')[-1].split('.')[0]
        self.ag_name = antigen.split('/')[-1].split('.')[0]

    """
    Metodo para criar tabela de restricao e criar o arquivo unambig. Para isso,
    e preciso renumerar a cadeia leve, dando um shift de 500 em sua numeracao 
    """
    def restrain_table(self):
        def renumber_pdb(input_pdb, output_pdb, shift):
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure('structure', input_pdb)

            for model in structure:
                for chain in model:
                    for residue in chain:
                        new_id = (' ', residue.id[1] + shift, ' ')
                        residue.id = new_id
            
            io = PDBIO()
            io.set_structure(structure)
            io.save(output_pdb)

        def combine_pdbs(chain1, chain2, output_pdb):
            with open(output_pdb, 'w') as outfile:
                with open(chain1) as infile1:
                    lines1 = infile1.readlines()
                    # Remover a linha "END" do primeiro PDB
                    if lines1[-1].startswith("END"):
                        lines1 = lines1[:-1]
                    outfile.writelines(lines1)
                
                with open(chain2) as infile2:
                    lines2 = infile2.readlines()
                    outfile.writelines(lines2)

        pdb_h = f'{self.ab_name}_H'
        pdb_h_path = os.path.join(self.pdb_files_path, f'{pdb_h}.pdb')
        pdb_l = f'{self.ab_name}_L'
        pdb_l_path = os.path.join(self.pdb_files_path, f'{pdb_l}.pdb')

        # Renumerar cadeias
        pdb_h_shift_path = f'{pdb_h_path[:-4]}_shift.pdb'
        pdb_l_shift_path = f'{pdb_l_path[:-4]}_shift.pdb'

        renumber_pdb(pdb_h_path, pdb_h_shift_path, 0)
        renumber_pdb(pdb_l_path, pdb_l_shift_path, 500)

        # Juntar arquivos PDB
        combined_pdb_path = os.path.join(self.pdb_files_path, f'antibody_{self.ab_name}.pdb')
        combine_pdbs(pdb_h_shift_path, pdb_l_shift_path, combined_pdb_path)

        # Executar script externo
        unambig_path = os.path.join(self.pdb_files_path, f'antibody_{self.pdb_name}-unambig.tbl')
        subprocess.run(["python", "restrain_bodies.py", combined_pdb_path], stdout=open(unambig_path, 'w'), check=True)

    def run_freesasa(self, pdb_file, rsa_file):
        command = f"freesasa {pdb_file} --format=rsa -o {rsa_file}"
        subprocess.run(command, shell=True, check=True)

    def parse_rsa(self, rsa_file, threshold=40):
        active_residues = []
        passive_residues = []
        with open(rsa_file, 'r') as file:
            for line in file:
                columns = line.split()
                if len(columns) in [13, 14]:  # Considera lines com 13 ou 14 colunas
                    sasa_all_atoms = float(columns[7]) if len(columns) == 13 else float(columns[8])
                    sasa_side_chain = float(columns[9]) if len(columns) == 13 else float(columns[10])
                    if sasa_all_atoms > threshold or sasa_side_chain > threshold:
                        residue = columns[3] if len(columns) == 13 else columns[4]
                        active_residues.append(residue)
                    else:
                        residue = columns[3] if len(columns) == 13 else columns[4]
                        passive_residues.append(residue)
        return active_residues, passive_residues

    def save_list(self, residues, filename):
        with open(filename, 'w') as file:
            file.write(" ".join(residues))
        

    """
    Metodo para identificar residuios ativos e passivos atraves do FREESASA.

    """
    def interface_map(self):
        pdb_file = os.path.join(self.pdb_files_path, f'{self.ab_name}.pdb')
        active_rsa_file = f"{self.pdb_files_path}/{self.ab_name}_active.rsa"
        passive_rsa_file = f"{self.pdb_files_path}/{self.ab_name}_passive.rsa"

        # Run FREESASA to calculate SASA
        self.run_freesasa(pdb_file, pdb_file, active_rsa_file)
        self.run_freesasa(pdb_file, pdb_file, passive_rsa_file)

        # Parse RSA files to get active and passive residues
        active_residues, passive_residues = self.parse_rsa(active_rsa_file)
        _, passive_residues = self.parse_rsa(passive_rsa_file)

        # Save active and passive residues to list files
        self.save_list(active_residues, f"{self.pdb_files_path}/{self.pdb_name}_active.list")
        self.save_list(passive_residues, f"{self.pdb_files_path}/{self.pdb_name}_passive.list")

        print(f"Active residues saved to {self.pdb_files_path}/{self.pdb_name}_active.list")
        print(f"Passive residues saved to {self.pdb_files_path}/{self.pdb_name}_passive.list")