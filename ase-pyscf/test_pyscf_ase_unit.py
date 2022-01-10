import ase.calculators.pyscf as pyscf
from ase import Atoms
from ase.lattice.cubic import Diamond
from numpy import array

# import ASE diamond
ase_atoms = Diamond(symbol='C', latticeconstant=3.5668)


# Test method that converts ase to pyscf atom

def test_ase_atoms_to_pyscf():
    pyscf_atom = pyscf.ase_atoms_to_pyscf(ase_atoms)
    
    expected = 
    
    
    assert observed == expected

'''
ase_atoms = Diamond(symbol='C', latticeconstant=3.5668)

diamond = pyase.ase_atoms_to_pyscf(ase_atoms)

print("ASE diamond", ase_atoms)
print("PYSCF diamond", diamond)
print(type(diamond))
'''