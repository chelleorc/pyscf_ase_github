"""
Take ASE Diamond structure, input into PySCF and run
"""

import pyscf.gto as pbcgto
import pyscf.dft as dft

from ase.calculators.PYSCF_calculator import ase_atoms_to_pyscf as pyscf_atom
from ase.lattice.cubic import Diamond

ase_atom=Diamond(symbol='C', latticeconstant=3.5668)
print(ase_atom.get_volume())

diamond_shape = pyscf_atom(ase_atom)
print(diamond_shape)

mf_hf_diamond = dft.RKS(diamond_shape)
mf_hf_diamond.xc = 'lda,vwn' # default
mf_hf_diamond = mf_hf.newton() # second-order algortihm
mf_hf_diamond.kernel() 



