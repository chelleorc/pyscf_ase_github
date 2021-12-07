"""Python-based Simulations Chemistry Framework"""

'''
This will use the PySCF framework to perform 
electronic structure calculations.
'''
import numpy as np
from ase.calculators.calculator import (Calculator, CalculatorError, 
                    CalculatorSetupError, all_changes, all_properties)



def ase_atoms_to_pyscf(ase_atoms):
    '''Convert ASE atoms to PySCF atom.

    Note: ASE atoms always use A.
    '''
    return [[atom.symbol, atom.position] for atom in ase_atoms]

atoms_from_ase = ase_atoms_to_pyscf


class PySCF(Calculator):
    implemented_properties = ['energy','forces']

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='PySCF', atoms=None, scratch=None, **kwargs):
        """Construct PySCF-calculator object.

        Parameters
        ==========
        label: str
            Prefix to use for filenames (label.in, label.txt, ...).
            Default is 'PySCF'.

        mfclass: PySCF mean-field class
        molcell: PySCF :Mole: or :Cell:
        """
        Calculator.__init__(self, restart=None, ignore_bad_restart_file=False,
                            label='PySCF', atoms=None, scratch=None, **kwargs)

        # TODO
        # This explicitly refers to "cell". How to refer
        # to both cell and mol together?

        self.mf=None
        self.initialize(**kwargs)


    # Check if ASE molecule begins with A or a
    def initialize(self, molcell, mf_class, mf_dict,mf_p=None):
        if not molcell.unit.startswith(('A','a')):
            raise RuntimeError("PySCF unit must be A to work with ASE")

        self.molcell=molcell
        self.mf_class=mf_class
        self.mf_dict=mf_dict
        self.mf_p=mf_p


    def set(self, **kwargs):
        changed_parameters = Calculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    # Set up ASE base calculator and calculate energy
    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=['positions', 'numbers', 'cell',
                                  'pbc', 'charges','magmoms']):

        Calculator.calculate(self, atoms)

        calc_molcell = self.molcell.copy()
        calc_molcell.atom = ase_atoms_to_pyscf(atoms)
        calc_molcell.a = np.asarray(atoms.cell)
        calc_molcell.build(None,None)
        self.mf = self.mf_class(calc_molcell)
        
        for key in self.mf_dict:
            self.mf.__dict__[key] = self.mf_dict[key]

        self.results['energy']=self.mf.scf()
        self.results['mf']=self.mf
        
        if th == 'mp2': # Need to fix this
            self._ham = mp.MP2(self._mf)
            hamkwargs = p['mp']
        elif th == 'ccsd':
            self._ham = cc.CCSD(self._mf)
            hamkwargs = p['cc']
        elif th == 'cisd':
            self._ham = ci.CISD(self._mf)
            hamkwargs = p['ci']
        else:
            self._ham = self._mf
            hamkwargs = dict()

        # Set Hamiltonian arguments
        # Is this the best way to do this?
        for key, val in hamkwargs.items():
            setattr(self._ham, key, val)

        # Energy/gradient scanner for more efficient evaluations
        self._escanner = None
        self._gscanner = None

        if 'forces' in properties:
            if self._gscanner is None:
                self._gscanner = self._ham.nuc_grad_method().as_scanner()
            e, g = self._gscanner(self._mol)
            self.results['forces'] = -(g * Hartree / Bohr).reshape((-1, 3))
        else:
            if self._escanner is None:
                self._escanner = self._ham.as_scanner()
            e = self._escanner(self._mol)
        '''  
        if 'forces' in properties:
            gf = self.mf.nuc_grad_method()
            gf.verbose = self.mf.verbose
            if self.p.mode.lower() == 'dft':
                gf.grid_response = True
            forces = -1. * gf.kernel() * (Ha / Bohr)
            totalforces = []
            totalforces.extend(forces)
            totalforces = numpy.array(totalforces)
            self.results['forces'] = totalforces
        '''
        self.results['energy'] = e * Hartree

def make_kpts(cell, nks):
    raise DeprecationWarning('Use cell.make_kpts(nks) instead.')


