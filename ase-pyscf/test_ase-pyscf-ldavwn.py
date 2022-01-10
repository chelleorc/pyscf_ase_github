"""
Take ASE Diamond structure, input into PySCF and run
"""

import numpy as np
import pyscf.pbc.gto as pbcgto
import pyscf.pbc.dft as pbcdft


from ase.calculators.pyscf import PySCF
from ase.units import kJ
from ase.utils.eos import EquationOfState
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

# Set up calculation 
mf_class = pbcdft.RKS
mf_dict = {'xc':'lda,vwn'}
args = {'molcell': molcell,'mf_class': mf_class, 'mf_dict': mf_dict}

# Pass args to pyscf calculator and calculate energy
calc = PySCF(**args)
calc.calculate(atoms=ase_atom)
