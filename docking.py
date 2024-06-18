import os
import subprocess
from argparse import ArgumentParser
from Bio.PDB import PDBParser, PDBIO, Select
from Bio import SeqIO

from haddock import Haddock

def write_log(log, log_file):
    with open(log_file, 'w') as file:
        for line in log:
            file.write(line + "\n")

"""
Função para aplicar limpeza em um arquivo PDB de modo a preparar para o Haddock.
Mantendos apenas linhas com ATOM, TER e END.

Parameters:
    pdb_file: [path + pdb_name.pdb]
"""
def cleaning_pdb(pdb_file: str):
    lines = []
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("TER") or line.startswith("END"):
                lines.append(line)

    with open(pdb_file, 'w') as out_file:
        for line in lines:
            out_file.write(line)

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

"""
Função para identificar se um arquivo pdb é anticorpo ou não.

Parameters:
    pdb_file: [path + pdb_name.pdb]

Returns:
    bool: True se for anticorpo; False caso seja antigeno.
"""
def is_ab_ag(pdb_file: str) -> bool:
    pdb_name = pdb_file.split('/')[-1].split('.')[0]

    # Caminho temp para os arquivos fasta e anarci
    output_path = os.path.join('temp', pdb_name)
    if not os.path.exists(output_path):
        if not os.path.exists('temp'):
            os.mkdir("temp")
        os.mkdir(output_path)

    fasta_path = extract_fasta(pdb_file=pdb_file, pdb_name=pdb_name, output_path=output_path)
    anarci_path = create_anarci_file(fasta_path=fasta_path, pdb_name=pdb_name, output_path=output_path)

    with open(anarci_path, 'r') as file:
        lines = file.readlines()

    # Não preciso ler todo o arquivo
    for line in lines:
        if line[0] != "#":
            return True if line[0] != "/" else False

"""
Função para identificar se o anticorpo é VH/VL, scFv ou nanocorpo a
partir do arquivo .anarci

Parameters:
    anarci_file: caminho do arquivo anarci

Returns:
    ab_type: tipo de construção do anticorpo; ERROR caso entrada não prevista
"""
def identify_ab_construction(anarci_file: str) -> str:
    with open(anarci_file, 'r') as file:
        lines = file.readlines()
    
    count_chain = 0
    for line in lines:
        if line[0] == "#":
            if "# Domain 1 of 1" in line:
                count_chain += 1
            elif "# Domain 1 of 2" in line:
                ab_type = 'scfv'
                return 0
    
    if count_chain == 1:
        ab_type = 'vhh'
    elif count_chain == 2:
        ab_type = 'vhvl'
    else:
        print("Entrada nao prevista no algoritmo")
        print(f"Total de cadeias contadas igual a {count_chain}.")
        return 'ERROR'
    return ab_type
    
"""
Função para modificar o nomeclativo: conversão Amber p/ formato pdb

Parameters:
    pdb_file: [path + pdb_name.pdb]

Returns:
    res_log: lista contendo as linhas que foram modificadas
"""
def modify_pdb_residues(pdb_file: str) -> list[str]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)

    log = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.resname in ['HIE', 'HID', 'HIP']:
                    original_line = f"Chain {chain.id}, Residue {residue.id}, Original {residue.resname}"
                    log.append(original_line)
                    residue.resname = 'HIS'
                elif residue.resname in ['CYX']:
                    original_line = f"Chain {chain.id}, Residue {residue.id}, Original {residue.resname}"
                    log.append(original_line)
                    residue.resname = 'CYS'
    
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_file)

    return log

"""
Função para modificar a numeração do arquivo pdb, removendo as letras que vem
apos a numeração e continuando a sequencia.

Parameters:
    pdb_file: [path + pdb_name.pdb]

Returns:
    files_path: caminho dos arquivos pdbs separados por cadeias
"""
def change_numeration_pdb(pdb_file: str, pdb_name: str) -> list[str]:
    parser     = PDBParser(QUIET=True)
    io         = PDBIO()
    structure  = parser.get_structure(pdb_name, pdb_file)
    pdb_chains = structure.get_chains()

    files_path  = [] # Contem caminho dos arquivos das cadeias
    # Remover cabeçalho =============================================
    for chain in pdb_chains:
        io.set_structure(chain)
        chain_file = structure.get_id() + "_" + chain.get_id() + ".pdb"

        files_path.append(chain_file)

        io.save(chain_file)

    for file_path in files_path:
        with open(file_path, 'r') as file:
            # Ler e grava todas as linhas do pdb numa lista
            lines = file.readlines()

        with open(file_path, 'w') as file:
            res_code   = ' '   # code for insertions of residues
            conversao  = False # indicativo para converter numeração
            qtd_code   = 0     # contador de quantos code apareceram

            for line in lines:
                new_line = ""

                # Verificar se está no fim do arquivo
                if line[0:3] == 'TER':
                    res_number_strip = line[23:26].strip()
                    count_space = 3 - len(res_number_strip)
                    new_line += count_space * ' '

                    file.write(f"{line[0:23]} {count_space * ' '}{str(int(res_number_strip) + qtd_code)} {line[28:]}")
                    break

                # nova linha copiada até Chain identifier
                new_line += line[:23] + ' '

                # Conversão de numeração (Sequencial,   letras)
                if line[26] != ' ':
                    if line[26] != res_code:
                        qtd_code += 1

                    res_code = line[26]
                    conversao = True
                elif line[26] == ' ' and conversao:
                    res_code = ' '
                    conversao = False

                # Adicionar espaços pela diferença de tamanho das strings
                res_number_strip = line[23:26].strip()
                count_space = 3 - len(res_number_strip)
                new_line += count_space * ' '

                new_line += str(int(res_number_strip) + qtd_code)

                # Preenche até final da linha
                new_line += ' ' + line[28:]

                file.write(new_line)
    
    return files_path

"""
Renomeia o nome da cadeia do arquivo de estrutura .pdb para o formato
HADDOCK. Cadeia pesada = H, Cadeia leve = L e VHH = H. 
Caso seja antigeno muda para alguma letra exceto H ou L

Parameters:
    pdb_file: [path + pdb_name.pdb]
    is_ab: True se for anticorpo; False se for antigeno
    ab_type: tipo de construção do anticorpo
    chain_id: id da cadeia com base do arquivo anarci
"""
def rename_chain(pdb_file: str, is_ab: bool, ab_type: str, chain_id: str):
    with open(pdb_file, 'r') as file:
        lines = file.readlines()

    new_line = ''
    with open(pdb_file, 'w') as file:
        for line in lines:
            if is_ab:
                if ab_type == 'vhh':
                    new_line = line[:21] + 'H' + line[23:]
                elif ab_type == 'scfv':
                    new_line = line[:21] + 'S' + line[23:]
                else:
                    new_line = line[:21] + chain_id + line[23:]
            else:
                new_line = line[:21] + chain_id + line[23:]
            
            file.write(new_line)

"""
Método para Realizar a reorganização das caideas (Ab / Ag).
Renomear o identificador das Cadeias, Renumerar e Fundir.

Parameters:
    anarci_path: caminho para o arquivo anarci
    pdb_name: nome do arquivo pdb
    is_ab: True se for anticorpo; False se for antigeno
    ab_type: tipo de cosntrução do anticorpo
"""
def rearrange_chains(anarci_path: str, pdb_name: str, is_ab: bool, ab_type: str='ag'):
    with open(anarci_path, 'r') as file:
        anarci_lines = file.readlines()

    structure_name = ''
    ag_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'I', 'J', 'K']
    count = 0
    aux = False # Indicativo para pegar o caminho do arquivo
    for line in anarci_lines:
        if not aux:
            structure_id = line.split(':')[1].strip()
            structure_name = f'{pdb_name}_{structure_id}.pdb'
            file_path = os.path.join('./', structure_name)
            aux = True

        # Fim da leitura do residuo
        if line[:2] == '//':
            if is_ab:
                rename_chain(file_path, is_ab=True, ab_type=ab_type, chain_id=chain_id)
            else:
                rename_chain(file_path, is_ab=False, ab_type=ab_type, chain_id=ag_list[count])
                count += 1
            
            aux = False
        else:
            if(line[0] != "#"):
                chain_id = line[0]

if __name__ == '__main__':
    # Argumentos de passagem para inicialização =========================================
    parser = ArgumentParser(
        description='Processo de Docking através do Haddock 2.4.',
        usage='python %(prog)s [file_1.pdb] [file_2.pdb] [options]'
    )
    parser.set_defaults(
        path    = './',
        dc_type = 'cego' 
    )

    parser.add_argument(
        '-p', '--path', metavar='STRING', type=str,
        help='Endereço da pasta dos arquivos de entrada.'
    )
    parser.add_argument(
        'input', metavar='STRING', type=str, nargs=2,
        help='Passe dois arquivos; 1 ab + 1 ag.'
    )
    parser.add_argument(
        '-t', '--dc_type', metavar='STRING', type=str,
        help="Indique o tipo de docking a ser executado: 'cego' ou 'dirigido'."
    )
    parser.add_argument(
    	'-at', '--active_res', metavar='PATH', type=list[str,str],
    	help='Lista contendo os residuos ativos da estrutura.'
    )
    parser.add_argument(
    	'-ps', '--passive_res', metavar='PATH', type=list[str,str],
    	help='Lista contendo os residuos passivos da estrutura.'
    )

    args = parser.parse_args()
    print(args.dc_type)
    # ===================================================================================

    print('''
    #####################################################
    #                 Script de DOCKING!                #
    #####################################################
    ''')

    # Carregar modulos ======================================== #
    print('Carregando modulos . . .')                           #
    subprocess.run(['bash', 'utils/load_modules.sh'])           #
    # ========================================================= #

    pdb_file1 = args.path + args.input[0]
    pdb_file2 = args.path + args.input[1]

    # Realizar a limpeza dos arquivos PDBs =====================
    print('>> Limpeza e tratamento de arquivos . . .')
    cleaning_pdb(pdb_file1)
    cleaning_pdb(pdb_file2)

    # Identificar quais arquivos são anticorpo e antigeno ======
    # Criar caminho temp para os arquivos fasta e anarci
    temp_path = os.path.join('./', 'temp')
    if not os.path.exists(temp_path):
        os.mkdir(temp_path)

    if is_ab_ag(pdb_file1):
        antibody = pdb_file1
        antigen  = pdb_file2
    else: 
        antibody = pdb_file2
        antigen  = pdb_file1

    is_ab_ag(pdb_file2)

    ab_name = antibody.split('/')[-1].split('.')[0]
    ag_name = antigen.split('/')[-1].split('.')[0]

    # Some processament ========================================
    # Criar caminho para os arquivos de log
    log_path = os.path.join('./', 'log')
    if not os.path.exists(log_path):
        os.mkdir(log_path)

    # Modificar nomeclativo dos residuos
    log = modify_pdb_residues(antibody)
    log_file = os.path.join(log_path, f'log_res_{ab_name}.txt')
    write_log(log, log_file)
    log = modify_pdb_residues(antigen)
    log_file = os.path.join(log_path, f'log_res_{ag_name}.txt')
    write_log(log, log_file)

    # Modificar numeração dos pdbs
    ab_chains_path = change_numeration_pdb(antibody, ab_name)
    ag_chains_path = change_numeration_pdb(antigen, ag_name)

    # Preparamento para o Haddock ==============================
    anarci_path = os.path.join(temp_path, f'{ab_name}')
    anarci_path = os.path.join(anarci_path, f'{ab_name}.anarci')
    ab_type = identify_ab_construction(anarci_file=anarci_path)

    rearrange_chains(anarci_path, pdb_name=ab_name, is_ab=True, ab_type=ab_type)

    anarci_path = os.path.join(temp_path, f'{ag_name}')
    anarci_path = os.path.join(anarci_path, f'{ag_name}.anarci')

    rearrange_chains(anarci_path, pdb_name=ag_name, is_ab=False)

    # Chamada do Haddock =======================================
    haddock = Haddock(
        pdb_files_path=args.path,
        antibody=antibody, 
        antigen=antibody, 
        ab_type=ab_type
    )    

    # Criar tabela de restricao ==========
    if haddock.ab_type == 'vhvl':
        print('>> Restrain table...\n')
        haddock.restrain_table()
        
    # Calcula residuos ati-pas ===========
    print('>> Calculando residuos ativos e passivos...\n')
    active_res  = args.active_res
    passive_res = args.passive_res
    haddock.interface_map(active_res, passive_res, args.dc_type)
    
