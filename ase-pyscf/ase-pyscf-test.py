"""
Take ASE Diamond structure, input into PySCF and run
"""

import pyscf.pbc.gto as pbcgto
import pyscf.pbc.dft as pbcdft
from pyscf.pbc import scf

from ase.calculators.PYSCF_calculator import ase_atoms_to_pyscf as pyscf_atom
from ase.calculators.PYSCF_calculator import PySCF, make_kpts
from ase.lattice.cubic import Diamond


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
args = {'molcell': molcell,'mf_class': mf_class, 'mf_dict': mf_dict}

#  Set up for k-point sampling
nk = [5,5,5]
kpts =molcell.make_kpts(nk)

calc = PySCF(**args)
calc.calculate(atoms=ase_atom)
print(kpts.shape)




