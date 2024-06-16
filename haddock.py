
# HADDOCK TUTORIAL:
# https://www.bonvinlab.org/education/HADDOCK24/HADDOCK24-local-tutorial/
class Haddock:
    # lista de proteinas a serem modificadas
    PROTEINS = ['HIE', 'HID', 'HIP', 'CYX']

    def __init__(self, pdb_files_path, antibody, antigen):
        self.pdb_files_path = pdb_files_path

        self.ab = antibody
        self.ag = antigen

        # Indicadores para referir ao tipo do Ab: ['scfv' | 'vhvl' | 'vhh'].
        self.ab_type = ''