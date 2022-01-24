from distutils.command.build import build
import ase.calculators.pyscf as pyscf
from ase import Atoms
from ase.lattice.cubic import Diamond
import pytest 
import sys
from pyscf.pbc import scf
from pyscf.pbc import gto as pbcgto
import numpy as np
import math


# Test package imports correctly
def test_pyscf_import():
    """Will pass if import statement works"""
    assert "pyscf" in sys.modules

def test_Atoms_import():
    """Will pass if import statement works"""
    assert "Atoms" in sys.modules

# Get ASE diamond and changed to PySCF atom
ase_atoms = Diamond(symbol='C', latticeconstant=3.5668)
pyscf_ase_atom = pyscf.ase_atoms_to_pyscf(ase_atoms)
#print(pyscf_ase_atom)


# Build pyscf diamond
cell_diamond = pbcgto.Cell()
cell_diamond = pbcgto.M(atom = [
                    ['C', np.array((0.,      0.,      0.))],    
                    ['C', np.array((0.8917,  0.8917,  0.8917))],
                    ['C', np.array((1.7834,  1.7834,  0.))],    
                    ['C', np.array((2.6751,  2.6751,  0.8917))],
                    ['C', np.array((1.7834,  0.,      1.7834))],          
                    ['C', np.array((2.6751,  0.8917,  2.6751))],
                    ['C', np.array((0.,      1.7834,  1.7834))],
                    ['C', np.array((0.8917,  2.6751,  2.6751))]
                        ],
                    basis = {'C': 'gth-szv'},
                    pseudo = 'gth-pade',
                    a = np.eye(3) * 3.5668
                    )

cell_diamond.build()

#print("PYSCF atom:")
#print(cell_diamond.atom)
print(pyscf_ase_atom == cell_diamond.atom[0:])
# np.all(math.isclose(pyscf_ase_atom, cell_diamond.atom))


'''
def test_ase_atoms_to_pyscf():
    """ Test that calculated_ase_atoms_to_pyscf_vol function calculates expected pyscf_atom_vol.
    This means that the ase_atoms_to_pyscf function is also what we expect. """

    # Convert ase atoms to pyscf atoms
    pyscf_atom = pyscf.ase_atoms_to_pyscf(ase_atoms)

    expected_format = cell_diamond._atom
   
   # ASE diamond volume
    #calculated_ase_atoms_to_pyscf_vol = ase_atoms.get_volume() 
    
    assert expected_format == pyscf_atom
'''
