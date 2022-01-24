from distutils.command.build import build
from tabnanny import verbose
import ase.calculators.pyscf as pyscf
from ase import Atoms
from ase.lattice.cubic import Diamond
import pytest 
import sys
import pyscf.pbc.dft as pbcdft
from pyscf.pbc import scf
from pyscf.pbc import gto as pbcgto
import numpy as np


# Test package imports correctly
def test_pyscf_import():
    """Will pass if import statement works"""
    assert "pyscf" in sys.modules

# Get ASE diamond and changed to PySCF atom
ase_atoms = Diamond(symbol='C', latticeconstant=3.5668)
pyscf_ase_atom = pyscf.ase_atoms_to_pyscf(ase_atoms)


# Note: wrap this code in a function
def pyscf_atom():
    cell_pyscf_atom = pbcgto.Cell()
    cell_pyscf_atom = pbcgto.M(atom = [
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
                    verbose = 5,
                    a = np.eye(3) * 3.5668
                    )
    cell_pyscf_atom.build()
    return cell_pyscf_atom.atom


def test_assert_atoms_are_equal():
    """ Assert that calculated_ase_atoms_to_pyscf function calculates expected pyscf_atom.
    This means that the ase_atoms_to_pyscf function is also what we expect. 
    
    array1: [["C", np.array([x, y, z])]

    array2: [["C", np.array([x, y, z])]

    """
    array1 = pyscf_ase_atom
    array2 = pyscf_atom()
    for a, b in zip(array1, array2):
        for x, y in zip(a, b):
            assert np.all(x == y), f"Molcell {x} != {y}! Check pyscf atom"

#test_assert_atoms_are'_equal(pyscf_ase_atom, pyscf_atom())
if __name__=="__main__":
    test_assert_atoms_are_equal()
    print("All modules executed when ran directly")
    
'''
def test_assert_pyscf_energy_calculator(pyscf_ase_atom, cell_pyscf_atom):
    # Calculate energy of ase_to_pyscf atom using pyscf.py calculator
    molcell = cell_pyscf_atom
    mf_class = pbcdft.RKS
    mf_dict = {'xc':'lda,vwn'}
    args = {'molcell': molcell,'mf_class': mf_class, 'mf_dict': mf_dict}
    calc = pyscf.PySCF(**args)
    pyscf_ase_energy = calc.calculate(atoms=ase_atoms)

    # Calculate energy of pyscf atom using pyscf package
    cell_pyscf_atom_energy = pbcdft.RKS(cell_pyscf_atom).run()
    #assert cell_pyscf_atom_energy == pyscf_ase_energy
    print("pyscf_ase_energy: ", type(pyscf_ase_energy))
    print("cell_pyscf_atom_energy: ", type(cell_pyscf_atom_energy))

test_assert_pyscf_energy_calculator(pyscf_ase_atom, cell_pyscf_atom)
'''