"""Python-based Simulations Chemistry Framework"""


from ase.calculators.calculator import Calculator, CalculatorError, CalculatorSetupError, all_changes, all_properties
from ase.lattice.cubic import Diamond
from pyscf import scf
from pyscf.pbc import gto as pbcgto
from pyscf.pbc.tools.pyscf_ase import ase_atoms_to_pyscf



# Keys to be passed to hf.SCF objects
_scf_keys = [
        'chkfile',
        'conv_tol',
        'conv_tol_grad',
        'max_cycle',
        'init_guess',
        'diis',
        'diis_space',
        'diis_start_cycle',
        'diis_file',
        'level_shift',
        'direct_scf',
        'direct_scf_tol',
        'callback',
        'conv_check',
        'check_convergence',
        ]

# This will be converted to general module that converts ase atoms to pyscf atoms
ase_atom=Diamond(symbol='C', latticeconstant=3.5668) # Placeholder for ase atom
cell = pbcgto.Cell()
cell.atom=ase_atoms_to_pyscf(ase_atom) #function that converts ase to pyscf atoms and is needed for pyscf calculator
cell.a=ase_atom.cell
cell.basis = 'gth-szv'
cell.pseudo = 'gth-pade'
cell.build()

# This might be the generalized version needed to convert ase atoms to pyscf atoms
def init_geo(mf, atoms):
    # convert ASE structural information to PySCF information
    if atoms.pbc.any():
        cell = mf.cell.copy()
        cell.atom = ase_atoms_to_pyscf(atoms)
        cell.a = atoms.cell.copy()
        cell.build()
        mf.reset(cell=cell.copy())
    else:
        mol = mf.mol.copy()
        mol.atom = ase_atoms_to_pyscf(atoms)
        mol.build()
        mf.reset(mol=mol.copy())


# Splits parameters dict into a dict of dicts, where each
# top-level key indicates the pyscf object to which the key pertains.
def _split_params(params):
    out = dict(gto=dict(),
               scf=dict())

    return out


class PySCF(Calculator):
    implemented_properties = ['energy']

    default_parameters = dict(theory=None,
                              verbose=False,
                              nthreads=1,
                              df=False)

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='pyscf', atoms=None, directory='.', **kwargs):
        
        """Construct PySCF-calculator object.
        Parameters
        ==========
        label: str
            Prefix to use for filenames (label.in, label.txt, ...).
            Default is 'PySCF'.
        mfclass: PySCF mean-field class
        molcell: PySCF :Mole: or :Cell:
        """
        Calculator.__init__(self, restart, ignore_bad_restart_file, label,
                            atoms, directory, **kwargs)
        
        th = self.parameters.theory
        if th is not None and th.lower() == 'ccsd(t)':
            raise ValueError("CCSD(T) calculations currently not supported.")
            self._mol = None

    # Check if there is an atom
    def _pyscf_init(self, atoms):
        if atoms is None:
            return

        npbc = atoms.pbc.sum()
        p = _split_params(self.parameters)

        if npbc == 2 and atoms.pbc[2]:
            raise CalculatorSetupError("Atoms object has 2d-periodicity, but "
                                       "the c-vector direction is periodic. "
                                       "Re-orient your system so that the "
                                       "aperiodic dimension corresponds to "
                                       "the c-vector.")
        elif npbc == 1 and not atoms.pbc[0]:
            raise CalculatorSetupError("Atoms object has 1d-periodicity, but "
                                       "the a-vector direction is not "
                                       "periodic. Re-orient your system so "
                                       "that the periodic dimension "
                                       "coresponds to the a-vector.")

        # Set up Hartree-Fock dictionary_
        _mfs = dict(hf=scf.HF)

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
       
