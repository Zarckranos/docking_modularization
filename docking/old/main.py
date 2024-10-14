import os
import glob

from argparse import ArgumentParser
from Bio.PDB import PDBParser, PDBIO

from antibody import Antibody
#from haddock import Condor_Docking
from haddock import Haddock

# CONSTANTES ==================
PDB_FILES_PATH = 'pdb_files'
# =============================

"""
Função responsável por identifcar os tipos de cadeias do anticorpo. VH/VL; VHH; scFv
  :input:
  :output:
"""
def identify_ab_chains():
    ...

if __name__ == "__main__":
    print('Script de HADDOCK!')
    # Criando lista com todas as estruturas
    all_pdbs = glob.glob(PDB_FILES_PATH + '/*.pdb')
    for pdb_file in all_pdbs:
        # extrair nome do arquivo sem a extensao
        pdb_name = pdb_file.split('/')[-1].split('.')[0]
        print(f'>> file: {pdb_name}')

        #cd_docking = Condor_Docking(PDB_FILES_PATH, pdb_name, pdb_file)
        cd_docking = Haddock(PDB_FILES_PATH, pdb_name, pdb_file)

        # Realizar limpeza do arquivo .pdb ==============================
        print('>> Cleaning pdb file...\n')
        cd_docking.cleaning()

        # Realizar identificação de cadeias =============================
        print('>> Identifieng antibody chain...')
        if cd_docking.identify_ab_chain(): continue
        print(f'  -> antibody of type {cd_docking.type.upper()}\n')
        
        # Reorganizar cadeias ===========================================
        print('>> Rearranging chains...\n')
        cd_docking.rearrange_chains()

        # Calcula residuos ati-pas ======================================
        print('>> Calculando residuos ativos e passivos...\n')
        cd_docking.interface_map()

        # Criar tabela de restricao =====================================
        if cd_docking.type == 'vhvl':
            print('>> Restrain table...\n')
            cd_docking.restrain_table()
        
        cd_docking.create_run_param()
        cd_docking.create_run_haddock_sh()
        cd_docking.create_haddock_condor_sh()

        print('>> FINISHED!\n\n\n')
        # ===============================================================