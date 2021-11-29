"""
Take ASE Diamond structure, input into PySCF and run
"""

import pyscf.pbc.gto as pbcgto
import pyscf.pbc.dft as pbcdft

from ase.calculators.PYSCF_calculator import ase_atoms_to_pyscf as pyscf_atom
from ase.calculators.PYSCF_calculator import PySCF
from ase.lattice.cubic import Diamond

# Use cubic diamond structure from ASE and
# calculate its volume
ase_atom=Diamond(symbol='C', latticeconstant=3.5668)
print(ase_atom.get_volume())

'''
# Convert ASE diamond into PySCF diamond
diamond_shape = pyscf_atom(ase_atom)
print(diamond_shape)
'''


# Calculator PySCF energy

# Note: This is what is suppose to occur in the PYSCF_calculator
cell = pbcgto.Cell() # Modify cell object with periodic boundary conditions, esp for crystals
cell.verbose = 5

### This is a step-by-step way of building a crystal ###
#cell.atom = diamond_shape # holds ase-turned-pyscf diamond
cell.a = ase_atom.cell # .a initializes crystal shape using lattice vectors --> molcell 
cell.basis = 'gth-szv' # a basis is a coordinate system (i.e. x-y coord.)
cell.pseudo = 'gth-pade' # this is a pseudopotential, which is like a set energy for the crystal (since there are no interactions with other molecules, i.e. deformations)
#cell.build() # constructs the entire crystal with all previous component parts


# Rerence: https://github.com/pyscf/pyscf/blob/master/examples/pbc/09-talk_to_ase.py 
mf_class = pbcdft.RKS

mf_dict = {'xc', 'lda,vwn'}

# Call PySCF class for access to functions
# Use calculate to build diamond and calculate diamond

initialize_atom = ase_atom.set_calculator(PySCF.initialize(cell,mf_class,mf_dict)) # TypeError for mf_dict

print(initialize_atom)



