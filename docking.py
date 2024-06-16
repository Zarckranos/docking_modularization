import os
import subprocess
from argparse import ArgumentParser
from Bio.PDB import PDBParser, PDBIO, Select
from Bio import SeqIO

from haddock import Haddock

"""
Função para aplicar limpeza em um arquivo PDB de modo a preparar para o Haddock.
Mantendos apenas linhas com ATOM, TER e END.

Parameters:
    pdb_file: [path + pdb_name.pdb]
"""
def cleaning_pdb(pdb_file: str):
    parser = PDBParser(QUIET=True)
    io = PDBIO()

    # Carregar a estrutura PDB
    structure = parser.get_structure('structure', pdb_file)

    class SelectAtoms(Select):
        def accept_atom(self, atom):
            return True

        def accept_line(self, line):
            return line.startswith("ATOM") or line.startswith("TER") or line.startswith("END")
    
    # Gravar a estrutura filtrada
    io.set_structure(structure)
    io.save('teste.pdb', select=SelectAtoms())

"""
Função para extrair e criar arquivo com a sequencia fasta do pdb.

Parameters:
    pdb_file: [path + pdb_name.pdb]
    pdb_name: nome do arquivo pdb
    output_path: caminho da pasta para os arquivos fasta

Returns:
    fasta_path: caminho de pasta contendo os arquivos fasta
"""
def extract_fasta(pdb_file: str, pdb_name: str, output_path: str) -> str:
    # Retorna iterador da seq. do arquivo de entrada
    with open(pdb_file) as handle:
        sequences = list(SeqIO.parse(handle, "pdb-atom"))

    fasta_file = f'{pdb_name}.fasta'
    fasta_path = os.path.join(output_path, fasta_file)
    # Criar arquivo .fasta na subpasta './fasta_files'
    with open(fasta_path, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")

    return fasta_path

"""
Função para criar chamada do ANARCI e criar arquivo .anarci do pdb

Parameters:
    fasta_path: caminho da pasta contendo o arquivo fasta
    pdb_name: nome do arquivo pdb
    output_path: caminho da pasta para os arquivo anarci

Returns:
    anarci_path: caminho da pasta contendo os arquivo anarci
"""
def create_anarci_file(fasta_path: str, pdb_name: str, output_path: str) -> str:
    anarci_file = subprocess.run(f"ANARCI -i {fasta_path}", shell=True, capture_output=True).stdout.decode('utf-8')

    anarci_name = f'{pdb_name}.anarci'
    anarci_path = os.path.join(output_path, anarci_name)
    with open(anarci_path, 'w') as file:
        file.write(anarci_file)

    return anarci_path


if __name__ == '__main__':
    # Argumentos de passagem para inicialização =========================================
    parser = ArgumentParser(description='Processo de Docking através do Haddock 2.4.')
    parser.set_defaults(
        path = './'
    )

    parser.add_argument(
        '-p', '--path', metavar='STRING', type=str,
        help='Endereço da pasta dos arquivos de entrada.'
    )
    parser.add_argument(
        'input', metavar='STRING', type=str, nargs=2,
        help='Passe dois arquivos; 1 ab + 1 ag.'
    )

    args = parser.parse_args()
    # ===================================================================================

    print('Script de DOCKING!')
    # Realizar a limpeza dos arquivos PDBs =====================
    cleaning_pdb(args.path + args.input[0])
    # cleaning_pdb(args.path + args.input[1])

    # Identificar quais arquivos são anticorpo e antigeno ======
    # Criar caminho temp para os arquivos fasta e anarci
    temp_path = os.path.join('./', 'temp')
    if not os.path.exists(temp_path):
        os.mkdir(temp_path)

    extract_fasta(pdb_file=args.input[0], output_path=temp_path)
