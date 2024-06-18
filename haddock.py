import os
import subprocess

from Bio.PDB import PDBParser, PDBIO
from Bio import SeqIO

# HADDOCK TUTORIAL:
# https://www.bonvinlab.org/education/HADDOCK24/HADDOCK24-local-tutorial/
class Haddock:
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

        def combine_pdbs(file1, file2, output_pdb):
            with open(output_pdb, 'w') as output_file:
                with open(file1) as infile1:
                    lines1 = infile1.readlines()
                    # Remover a linha "END" do primeiro PDB
                    if lines1[-1].startswith("END"):
                        lines1 = lines1[:-1]
                    output_file.writelines(lines1)
                
                with open(file2) as infile2:
                    lines2 = infile2.readlines()
                    output_file.writelines(lines2)

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
        self.unambig_path = os.path.join(self.pdb_files_path, f'antibody_{self.pdb_name}-unambig.tbl')
        subprocess.run(["python", "restrain_bodies.py", combined_pdb_path], stdout=open(unambig_path, 'w'), check=True)

    def run_freesasa(self, pdb_file, rsa_file):
        command = f"freesasa {pdb_file} --format=rsa {rsa_file}"
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
    def interface_map(self, active_res, passive_res, dock_type):
        if self.ab_type == 'vhvl':
            # (antibody_{self.ab_name}.pdb) -> arquivo combinado com H-chian e L-chain
            pdb_file = os.path.join(self.pdb_files_path, f'antibody_{self.ab_name}.pdb')
        else:
            pdb_file = os.path.join(self.pdb_files_path, f'{self.ab_name}.pdb')
        active_rsa_file = f"{self.pdb_files_path}/{self.ab_name}_active.rsa"
        #passive_rsa_file = f"{self.pdb_files_path}/{self.ab_name}_passive.rsa"

        # Run FREESASA to calculate SASA
        self.run_freesasa(pdb_file, pdb_file, active_rsa_file)
        #self.run_freesasa(pdb_file, pdb_file, passive_rsa_file)

        # Parse RSA files to get active and passive residues
        active_residues, passive_residues = self.parse_rsa(active_rsa_file)
        #_, passive_residues = self.parse_rsa(passive_rsa_file)

        # Verificar residuos passados pelo usuário
        if active_res:
            ...
        if passive_res:
            ...

        # Save active and passive residues to list files
        self.save_list(active_residues, f"{self.pdb_files_path}/{self.pdb_name}_active.list")
        # self.save_list(passive_residues, f"{self.pdb_files_path}/{self.pdb_name}_passive.list")

        print(f"Active residues saved to {self.pdb_files_path}/{self.pdb_name}_active.list")

    """
    Método para geração de arquivos de restrição ambig entre ag e ab.
    """
    def ambig_table(self):
        command = f"active-passive-to-ambig.py antibody_{pdb_name}.list"
        subprocess.run(command, shell=True, check=True)

    """
    Método para criar arquivo run.param.
    """
    def create_rum_param(self):
        segid_1 = ''
        if ab.type == 'vhvl':
            segid_1 = 'H,L'
        elif ab.type == 'scfv':
            sigid_1 = 'S'
        elif ab.type == 'vhh':
            segid_1 = 'h'

        run_param_content = f"""\
AMBIG_TBL={self.ambig_path}
HADDOCK_DIR=/usr/local/softwares/haddock/2.4/
N_COMP=2
PDB_FILE1={self.ab}
PDB_FILE2={self.ag}
PROJECT_DIR={self.pdb_files_path}
PROT_SEGID_1={segid_1}
PROT_SEGID_2=A
RUN_NUMBER=1
UNAMBIG_TBL={self.unambig.path}
        """

        with open('run.param', 'w') as file:
            file.write(run_param_content)

    def create_haddock_condor_sub(self):
        content = f"""
Universe = vanilla
Environment = SHELL=/usr/bin/bash;PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
Executable = haddock-condor.sh

output = job.out.\$(Process)
error = job.error.\$(Process)
log = job.log.\$(Process)

requirements = (Machine != \"prainha\")

#Transferring Input Files
should_transfer_files = YES
when_to_transfer_output = ON_EXIT

transfer_input_files = ${folder}/${pdb_cdr}/antigeno_otimizado.pdb, ${folder}/${pdb_cdr}/haddock-condor.sh, ${folder}/${pdb_cdr}/anticorpo_otimizado.pdb, ${folder}/${pdb_cdr}/run.param, ${folder}/${pdb_cdr}/$list_ag, ${folder}/${pdb_cdr}/antibody_${reduced_name}-residuos-das-CDR.list, ${folder}/${pdb_cdr}/antibody-antigen-ambig.tbl, ${folder}/${pdb_cdr}/antibody_${reduced_name}-unambig.tbl, ${folder}/${pdb_cdr}/run-haddock.sh, ${folder}/${pdb_cdr}/antigeno_HISTIDINAS.txt, ${folder}/${pdb_cdr}/anticorpo_HISTIDINAS.txt

#Transferring Output Files
should_transfer_files = yes
when_to_transfer_output = ON_EXIT
transfer_output_files = run1

request_cpus = 12
Queue 1
        """

        with open('haddock-condor.sub', 'w') as file:
            file.write(content)

    """
    Método para criar arquivo hrun-addock.sh
    """
    def create_run_haddock_sh(self):
        content = """
#!/bin/bash
#sleep 15m

#source /usr/share/modules/init/bash
module load Haddock/2.4
echo $PATH
echo $LD_LIBRARY_PATH
echo $PYTHONPATH
echo $HADDOCK

python2 \$HADDOCK/Haddock/RunHaddock.py
        """

        with open('run-haddock.sh', 'w') as file:
            file.write(content)

