"""
Take ASE Diamond structure, input into PySCF and run
"""

import pyscf.pbc.gto as pbcgto
import pyscf.pbc.dft as pbcdft
from pyscf import lib, gto, scf, dft, mp, cc, ci

from ase.calculators.PYSCF_calculator import ase_atoms_to_pyscf as pyscf_atom
from ase.calculators.PYSCF_calculator import PySCF
from ase.lattice.cubic import Diamond

class parameters():
    # holds the calculation mode and user-chosen attributes of post-HF objects
    def __init__(self):
        self.mode = 'hf'
    def show(self):
        print('------------------------')
        print('calculation-specific parameters set by the user')
        print('------------------------')
        for v in vars(self):
            print('{}:  {}'.format(v,vars(self)[v]))
        print('\n\n')

# Use cubic diamond structure from ASE and
# calculate its volume
ase_atom=Diamond(symbol='C', latticeconstant=3.5668)
print(ase_atom.get_volume())


### Build the structure ###
    # Reference: https://github.com/pyscf/pyscf/blob/master/examples/pbc/09-talk_to_ase.py 
molcell = pbcgto.Cell() # periodic boundary conditions, esp for crystals --> molcell
molcell.verbose = 5
molcell.a = ase_atom.cell # initialize crystal shape using lattice vectors  
molcell.basis = 'gth-szv' # initialize basis
molcell.pseudo = 'gth-pade' # initialize pseudopotential


### Set up calculation ###
mf_class = pbcdft.RKS # --> mf_class
mf_dict = {'xc':'lda,vwn'} # -->mf_dict

index = 4

mf_p = parameters()
mf_p.mode = ['hf','mp2','cisd','ccsd','ccsd(t)'][index]
mf_p.verbose = 5
mf_p.show()

molcell.verbose = mf_p.verbose

args = {'molcell': molcell,'mf_class': mf_class, 'mf_dict': mf_dict,'mf_p': mf_p}



### Call PySCF class for access to functions ###
# calc = PySCF(molcell=molcell, mf_class=mf_class, mf_dict=mf_dict)

calc = PySCF(**args)
calc.calculate(atoms=ase_atom)


