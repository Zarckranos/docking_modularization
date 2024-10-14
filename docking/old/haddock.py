import os
import glob
import subprocess

from Bio.PDB import PDBParser, PDBIO
from Bio import SeqIO

# HADDOCK TUTORIAL
# https://www.bonvinlab.org/education/HADDOCK24/HADDOCK24-local-tutorial/
class Haddock:
    # lista de proteinas a serem modificadas
    PROTEINS = ['HIE', 'HID', 'HIP', 'CYX']

    def __init__(self, pdb_files_path, pdb_name, pdb_file):
        self.pdb_files_path = pdb_files_path

        self.pdb_name = pdb_name
        self.pdb_file = pdb_file

        # Indicadores para referir ao tipo do Ab: ['scfv' | 'vhvl' | 'vhh'].
        self.ab_type = ''

    """
    Método para extrair e criar arquivo com a sequencia fasta do arquivo pdb.
     :param path: caminho da pasta para os arquivos fasta
     :output: 
    """
    # @
    def extract_fasta(self, path):
        # Retorna iterador da seq. do arquivo de entrada
        with open(self.pdb_file) as handle:
            sequences = list(SeqIO.parse(handle, "pdb-atom"))

        fasta_file = f'{self.pdb_name}.fasta'
        self.fasta_path = os.path.join(path, fasta_file)
        # Criar arquivo .fasta na subpasta './fasta_files'
        with open(self.fasta_path, "w") as output_handle:
            SeqIO.write(sequences, output_handle, "fasta")

    """
    Método para criar chamada do ANARCI e criar arquivo .anarci do pdb
    """
    def create_anarci_file(self, path):
        anarci_file = subprocess.run(f"ANARCI -i {self.fasta_path}", shell=True, capture_output=True).stdout.decode('utf-8')

        anarci_name = f'{self.pdb_name}.anarci'
        self.anarci_path = os.path.join(path, anarci_name)
        with open(self.anarci_path, 'w') as file:
            file.write(anarci_file)

    """
    Método responsável por identificar se o anticorpo é VH/VL, scFv ou VHH
    a partir do arquivo .anarci.
    Retorna:
      int: Retorna 1 se a funcao encontrar alguma falha
    """
    def identify_ab_chain(self):
        # Caminho temp para os arquivos fasta e anarci
        path = os.path.join('temp', self.pdb_name)
        if not os.path.exists(path):
            if not os.path.exists('temp'):
                os.mkdir("temp")
            os.mkdir(path)

        self.extract_fasta(path)
        self.create_anarci_file(path)

        with open(self.anarci_path, 'r') as file:
            lines = file.readlines()

        # Identificar o tipo de anticorpo com base no 
        count_chain = 0
        for line in lines:
            if line[0] == "#":
                if "# Domain 1 of 1" in line:
                    count_chain += 1
                elif "# Domain 1 of 2" in line:
                    self.ab_type = 'scfv'
                    return 0
        
        if count_chain == 1:
            self.ab_type = 'vhh'
        elif count_chain == 2:
            self.ab_type = 'vhvl'
        else:
            print("Entrada nao prevista no algoritmo")
            print(f"Total de cadeias contadas igual a {count_chain}.")
            return 1
        return 0

    """
    Método que aplicará a limpeza do arquivo PDB
     :input:
     :output:
    """
    def cleaning(self):
        # Criar subpasta para colocar as estruturas do arquivo
        structures_fold = os.path.join(self.pdb_files_path, f'{self.pdb_name}_structures')
        if not os.path.exists(structures_fold):
            os.mkdir(structures_fold)

        # Remover Het_atm e Conexões ====================================
        # Extraindo as cadeias
        parser     = PDBParser(QUIET=True)
        io         = PDBIO()
        structure  = parser.get_structure(self.pdb_name, self.pdb_file)
        pdb_chains = structure.get_chains()

        files_path  = [] # Contem caminho dos arquivos das cadeias
        # Remover cabeçalho =============================================
        for chain in pdb_chains:
            io.set_structure(chain)
            chain_file = structure.get_id() + "_" + chain.get_id() + ".pdb"
            file_path = os.path.join(structures_fold, chain_file)

            files_path.append(file_path)

            io.save(file_path)

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

                        file.write(f"{line[0:23]}{count_space * ' '}{str(int(res_number_strip) + qtd_code)} {line[28:]}")
                        break

                    new_line += line[0:17]
                    # Nomeclativo: conversão Amber p/ formato PDB =======
                    res = line[17:20]
                    # @ POSSO USAR DICIONARIO
                    if res == 'HIE':
                        new_line += 'HIS'
                    elif res == 'HID':
                        new_line += 'HIS'
                    elif res == 'HIP':
                        new_line += 'HIS'
                    elif res == 'CYX':
                        new_line += 'CYS'
                    else:
                        new_line += res

                    # nova linha copiada até Chain identifier
                    new_line += ' ' + line[20:23]

                    # Conversão de numeração (Sequencial,   letras) =====
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

        # Função ========================================================

        # Remoção de cadeias desnecessárias ============================= (stand-by)
        return 0
    
    """
    Renomeia o nome da cadeia do arquivo de estrutura .pdb para o formato
    HADDOCK. Cadeia pesada = H, Cadeia leve = L e VHH = h. 
    Caso seja antigeno muda de forma sequencial a letra A,B,C,...
    """
    def rename_chain(self, file_path, is_ab, chain_id):
        with open(file_path, 'r') as file:
            lines = file.readlines()

        new_line = ''
        with open(file_path, 'w') as file:
            for line in lines:
                if is_ab:
                    if self.ab_type == 'vhh':
                        new_line = line[:21] + 'h' + line[23:]
                    else:
                        new_line = line[:21] + chain_id + line[23:]
                else:
                    new_line = line[:21] + chain_id + line[23:]
                
                file.write(new_line)
    """
    Método para Realizar a reorganização das caideas (Ab / Ag).
    Renomear o identificador das Cadeias, Renumerar e Fundir.
     :input:
     :output:
    """
    def rearrange_chains(self):
        # @ Renomeio o nome dos arquivos?
        with open(self.anarci_path, 'r') as file:
            anarci_lines = file.readlines()

        is_antibody = 0 # 0 - para não; 1 - para sim; 2 - possivel antigeno
        structure_name = ''
        ag_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'I', 'J', 'K']
        count = 0
        aux = False # Indicativo para pegar o caminho do arquivo
        for line in anarci_lines:
            if not aux:
                structure_id = line.split(':')[1].strip()
                structure_name = f'{self.pdb_name}_{structure_id}.pdb'
                file_path = os.path.join(self.pdb_files_path, f'{self.pdb_name}_structures/{structure_name}')
                aux = True

            # Fim da leitura do residuo
            if line[:2] == '//':
                if(is_antibody == 1):
                    self.rename_chain(file_path, is_ab=True, chain_id=chain_id)
                else:
                    self.rename_chain(file_path, is_ab=False, chain_id=ag_list[count])
                    count += 1
                
                aux = False
                is_antibody = 0
            else:
                # Se não for anticorpo, pega o nome da cadeia e atualizada isantibody para 2 - para checar se é um antigeno
                if(is_antibody != 1):
                    if(is_antibody == 0):
                        is_antibody = 2

                    # Caso contrario, não foi encontrado o fim do residuo e com certeza é um anticorpo
                    else:
                        is_antibody = 1

                else:
                    if(line[0] != "#"):
                        chain_id = line[0]
    

    """
    Metodo para criar tabela de restricao e criar o arquivo unambig. Para isso,
    e preciso renumerar a cadeia leve, dando um shift de 500 em sua numeracao 
    """
    def restrain_table(self):
        pdb_h = f'{self.pdb_name}_H'
        pdb_h_path = os.path.join(self.pdb_files_path, f'{self.pdb_name}_structures/{pdb_h}')
        pdb_l = f'{self.pdb_name}_L'
        pdb_l_path = os.path.join(self.pdb_files_path, f'{self.pdb_name}_structures/{pdb_l}')

        with open(f'{pdb_h_path}_shift.pdb', 'w') as outfile:
            subprocess.run(["pdb_reres", "-1", f'{pdb_h_path}.pdb'], stdout=outfile, check=True)
        
        # Deslocamento de numeracao da cadeia leve (-500)
        with open(f"{pdb_l_path}_shift.pdb", 'w') as outfile:
            subprocess.run(["pdb_reres", "-500", f"{pdb_l_path}.pdb"], stdout=outfile, check=True)
        # Juncao de arquivos para PDB unico
        out_path = os.path.join(self.pdb_files_path, f'antibody_{self.pdb_name}.pdb')
        with open(out_path, 'w') as outfile:
            subprocess.run(["cat", f"{pdb_h_path}_shift.pdb", f"{pdb_l_path}_shift.pdb"], stdout=outfile, check=True) 

        unambig_path = os.path.join(self.pdb_files_path, f'antibody_{self.pdb_name}-unambig.tbl')
        with open(unambig_path, 'w') as outfile:
            subprocess.run(["python", "restrain_bodies.py", out_path], stdout=outfile, check=True)

    """
    Metodo para identificar residuios ativos e passivos atraves do FREESASA.

    """
    def interface_map(self):
        active_rsa_file = f"{self.pdb_name}_active.rsa"
        active_list_file = f"{self.pdb_name}_active.list"
        
        subprocess.run(["freesasa", f"{self.pdb_name}.pdb", "--format=rsa", f">{active_rsa_file}"], shell=True)
        
        awk_active_command = f"awk '{{if (NF==13 && ($7>40 || $9>40)) printf \"%s \",$3; if (NF==14 && ($8>40 || $10>40)) printf \"%s \",$4}}' {active_rsa_file} > {active_list_file}"
        subprocess.run(awk_active_command, shell=True)
        
        # Calcular SASA e identificar resíduos passivos
        passive_rsa_file = f"{self.pdb_name}_passive.rsa"
        passive_list_file = f"{self.pdb_name}_passive.list"
        
        subprocess.run(["freesasa", f"{self.pdb_name}.pdb", "--format=rsa", f">{passive_rsa_file}"], shell=True)
        
        with open(passive_list_file, "w") as file:
            file.write(" ")
        
        awk_passive_command = f"awk '{{if (NF==13 && ($7>40 || $9>40)) printf \"%s \",$3; if (NF==14 && ($8>40 || $10>40)) printf \"%s \",$4}}' {passive_rsa_file} >> {passive_list_file}"
        subprocess.run(awk_passive_command, shell=True)

    """
    Metodo para criar o arquivo run.param
    """
    def create_run_param(self, haddock_dir, pdb_file1, pdb_file2, prot_segid_1, prot_segid_2, run_number,  n_comp=2):
        ambig_tbl=f'antibody_{self.pdb_name}-unambig.tbl'
        project_dir = './'
        
        run_param_content = f"""\
AMBIG_TBL={ambig_tbl}
HADDOCK_DIR={haddock_dir}
N_COMP={n_comp}
PDB_FILE1={pdb_file1}
PDB_FILE2={pdb_file2}
PROJECT_DIR={project_dir}
PROT_SEGID_1={prot_segid_1}
PROT_SEGID_2={prot_segid_2}
RUN_NUMBER={run_number}
    """
        
        with open("run.param", "w") as file:
            file.write(run_param_content)

    def create_haddock_condor_sub(self, folder, pdb_cdr, list_ag, reduced_name):
        content = f"""\
Universe = vanilla
Environment = SHELL=/usr/bin/bash;PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
Executable = haddock-condor.sh

output = job.out.$(Process)
error = job.error.$(Process)
log = job.log.$(Process)

requirements = (Machine != "prainha")

# Transferring Input Files
should_transfer_files = YES
when_to_transfer_output = ON_EXIT

transfer_input_files = {folder}/{pdb_cdr}/antigeno_otimizado.pdb, {folder}/{pdb_cdr}/haddock-condor.sh, {folder}/{pdb_cdr}/anticorpo_otimizado.pdb, {folder}/{pdb_cdr}/run.param, {folder}/{pdb_cdr}/{list_ag}, {folder}/{pdb_cdr}/antibody_{reduced_name}-residuos-das-CDR.list, {folder}/{pdb_cdr}/antibody-antigen-ambig.tbl, {folder}/{pdb_cdr}/antibody_{reduced_name}-unambig.tbl, {folder}/{pdb_cdr}/run-haddock.sh, {folder}/{pdb_cdr}/antigeno_HISTIDINAS.txt, {folder}/{pdb_cdr}/anticorpo_HISTIDINAS.txt

# Transferring Output Files
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_output_files = run1

request_cpus = 12
Queue 1
    """

        with open("haddock-condor.sub", "w") as file:
            file.write(content)

    def create_run_haddock_sh(self, filename="run-haddock.sh"):
        content = """#!/bin/bash
#source /usr/share/modules/init/bash
module load Haddock/2.4
echo $PATH
echo $LD_LIBRARY_PATH
echo $PYTHONPATH
echo $HADDOCK

python2 $HADDOCK/Haddock/RunHaddock.py
    """
        with open(filename, "w") as file:
            file.write(content)

    def create_haddock_condor_sh(self, filename="haddock-condor.sh"):
        content = """#!/bin/bash
if [ -r /usr/share/modules/init/bash ]; then
    source /usr/share/modules/init/bash
fi

echo $SHELL
echo $PATH
echo $SHELL
echo $PATH
chmod +x run-haddock.sh
./run-haddock.sh

cd run1

sed -i -e 's,cpunumber_1=8,cpunumber_1=12,' run.cns
sed -i -e 's,autohis=true,autohis=false,' run.cns
sed -i -e 's,structures_0=1000,structures_0=1000,' run.cns
sed -i -e 's,structures_1=200,structures_1=250,' run.cns
sed -i -e 's,anastruc_1=200,anastruc_1=250,' run.cns
sed -i -e 's,waterrefine=200,waterrefine=250,' run.cns

######### obter numeros das histidinas ######################
cat << EOF > mudar_histidinas.py
import os

# Define directories and initialize histidine lists
HISE_anticorpo = []
HISD_anticorpo = []
HISE_antigeno = []
HISD_antigeno = []

# Get the root directory
diretorio_raiz = os.path.abspath(os.path.join(os.getcwd(), os.pardir))

# Read histidines for antibody
with open("../anticorpo_HISTIDINAS.txt", "r") as anticorpo:
    for linha in anticorpo:
        if "HISE" in linha:
            HISE_anticorpo.append(linha.split()[2])
        elif "HISD" in linha:
            HISD_anticorpo.append(linha.split()[2])

# Sort histidine lists
HISE_anticorpo.sort()
HISD_anticorpo.sort()

# Read histidines for antigen
with open("../antigeno_HISTIDINAS.txt", "r") as antigeno:
    for linha in antigeno:
        if "HISE" in linha:
            HISE_antigeno.append(linha.split()[2])
        elif "HISD" in linha:
            HISD_antigeno.append(linha.split()[2])

# Sort histidine lists
HISE_antigeno.sort()
HISD_antigeno.sort()

# Modify run.cns for histidines in antibody
for i, his in enumerate(HISE_anticorpo, 1):
    replace_str = f"hise_1_{i}=0"
    new_str = f"hise_1_{i}={his}"
    with open("./run.cns", "r") as cns:
        content = cns.read().replace(replace_str, new_str)
    with open("./run.cns", "w") as cns:
        cns.write(content)

for i, his in enumerate(HISD_anticorpo, 1):
    replace_str = f"hisd_1_{i}=0"
    new_str = f"hisd_1_{i}={his}"
    with open("./run.cns", "r") as cns:
        content = cns.read().replace(replace_str, new_str)
    with open("./run.cns", "w") as cns:
        cns.write(content)

# Modify run.cns for histidines in antigen
for i, his in enumerate(HISE_antigeno, 1):
    replace_str = f"hise_2_{i}=0"
    new_str = f"hise_2_{i}={his}"
    with open("./run.cns", "r") as cns:
        content = cns.read().replace(replace_str, new_str)
    with open("./run.cns", "w") as cns:
        cns.write(content)

for i, his in enumerate(HISD_antigeno, 1):
    replace_str = f"hisd_2_{i}=0"
    new_str = f"hisd_2_{i}={his}"
    with open("./run.cns", "r") as cns:
        content = cns.read().replace(replace_str, new_str)
    with open("./run.cns", "w") as cns:
        cns.write(content)

# Update the number of histidines in run.cns
with open("./run.cns", "r") as cns:
    content = cns.read()
    content = content.replace("numhise_1=0", f"numhise_1={len(HISE_anticorpo)}")
    content = content.replace("numhisd_1=0", f"numhisd_1={len(HISD_anticorpo)}")
    content = content.replace("numhise_2=0", f"numhise_2={len(HISE_antigeno)}")
    content = content.replace("numhisd_2=0", f"numhisd_2={len(HISD_antigeno)}")

with open("./run.cns", "w") as cns:
    cns.write(content)

EOF

echo "###########################"
echo $PATH    
chmod +rw mudar_histidinas.py
chmod +x mudar_histidinas.py
python3 mudar_histidinas.py

find ./ -type f | xargs -n1 sed -i 's,/usr/bin/grep,/usr/bin/grep,'
../run-haddock.sh

chmod +rw mudar_histidinas.py
chmod +x mudar_histidinas.py

module load compilers/gcc/10.2.0
condor_submit haddock-condor.sub
"""
        with open(filename, "w") as file:
            file.write(content)

        os.chmod(filename, 0o755)
        os.system(f"./{filename}")