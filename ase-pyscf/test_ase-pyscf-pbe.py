"""
Take ASE Diamond structure, input into PySCF and run
"""

import pyscf.pbc.gto as pbcgto
import pyscf.pbc.dft as pbcdft
from ase.calculators.PYSCF_calculator import PySCF
from ase.lattice.cubic import Diamond


# Use cubic diamond structure from ASE
ase_atom=Diamond(symbol='C', latticeconstant=3.5668)

### Build the structure ###
    # Reference: https://github.com/pyscf/pyscf/blob/master/examples/pbc/09-talk_to_ase.py 
molcell = pbcgto.Cell() 
molcell.verbose = 5
molcell.a = ase_atom.cell # initialize crystal shape using lattice vectors  
molcell.basis = 'gth-szv' 
molcell.pseudo = 'gth-pade'

### Set up calculation to pass arguments
# to PYSCF_calculator ###
mf_class = pbcdft.RKS 
mf_dict = {'xc':'pbe'} 
args = {'molcell': molcell,'mf_class': mf_class, 'mf_dict': mf_dict}

# Pass args and calculate
calc = PySCF(**args)
calc.calculate(atoms=ase_atom)