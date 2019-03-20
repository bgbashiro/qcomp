"""gate module exposes Gate classes and common instances and means to join them to create circuit using operators * and +

(*) :: joins the gates that acts on whole register. E.g. a this how to make gate that acts Hadamard on 1st and 3rd qbit

>>> HIH = H*I*H
>>> HIH.apply(mk_reg([ONE,ONE,PLUS])) 
0.00+0.00j       |000> 
0.00+0.00j       |001> 
0.71+0.00j       |010> 
0.00+0.00j       |011> 
0.00+0.00j       |100> 
0.00+0.00j       |101> 
-0.71+0.00j      |110> 
0.00+0.00j       |111> 

(+) :: the gates will be applied consecutively. It is important that they can act on registers of same size
"""
import numpy as np
from .qbit import *
from .qregister import *
from .utils import *
from functools import reduce

def construct_unitary_F(fdef, in_size=None):
    """Given fdef construct appropriate unitary transformation matrix, probably useful for grovers algorithm. This is how to create a function {0,1}^n->{0,1}. It assumes last bit is output while others are input
    Parameters
    ---
    fdef: list of inputs that will calculate to 1. E.g. ["00","10"] -> will prepare gate that operates on 3 qbits. If first 2 qbits are 00 or 10 will AND 3rd one is 0, it will turn last one to 1. For making it unitary, if last bit is 1 it will turn into 0 for same inputs.
    in_size: expected input size, only necessary when fdef is empty list, ie function always produces 0

    Returns
    ---
    unitary_gate
    """
    if not in_size:
        # when input size is not given it is inferred from fdef
        in_size = len(fdef[0])
    reg_size = in_size+1
    dummy_bases=QReg(in_size).bases
    cmatrix = np.zeros([2**reg_size,2**reg_size])
    for base, i_col in zip(dummy_bases, range(0,2**reg_size,2)):
        if base in fdef:
            target_state = int(base+"1",2)
            cmatrix[target_state,i_col] = 1
            target_state = int(base+"0",2)
            cmatrix[target_state,i_col+1] = 1
        else:
            target_state = int(base+"0",2)
            cmatrix[target_state,i_col] = 1
            target_state = int(base+"1",2)
            cmatrix[target_state,i_col+1] = 1
    return MGate(cmatrix,reg_size)

class Gate():
    """Generic gate interface
    
    instance variables
    -----------
    *qreg_size* : number of qubits suitable for gate (e.g. Hadamard -> 1, CNOT -> 2)

    methods
    -----------------
    apply

    pre-defined gate instances:
    --------------------------
    **ID** :  identity gate, does not do anything to register<br/>
    **Hadamard** : Hadamard gate
    """
    def __init__(self,qreg_size):
        self.qreg_size = qreg_size

    def __add__(self, other):
        ss = Sequence(self.qreg_size)
        ss.add(self)
        ss.add(other)
        return ss


    def apply(self, qreg):
        """Apply the gate to the given register 
        Parameters
        -----------
        **qreg** : QReg instance 

        Returns
        ----------
        new_qreg : QReg object created after transformation
        """
        assert qreg.nbits == self.qreg_size, "This gate cannot be applied to register of size {}. Expected size: {}".format(qreg.nbits, self.qreg_size)
        raise NotImplementedError("Apply method not implemented yet!")

class Sequence(Gate):
    """Sequence of gates that will be applied in order. Note that this class is not supposed to be used on its own. Use + operator (it uses Sequence under the hood)
    """

    def __init__(self, qreg_size):
        self.seq = []
        super(Sequence, self).__init__(qreg_size)
    
    def add(self, other):
        if isinstance(other, Sequence):
            self.seq += other.seq
        else:
            self.seq += [other]

    def apply(self, qreg):
        tqreg = qreg.copy()
        for g in self.seq:
            tqreg = g.apply(tqreg)
        return tqreg

class MGate(Gate):
    """Gate that is formed by providing matrix form
    """

    def __init__(self, matrix, qreg_size):
        """ Initalize a gate from matrix
        *matrix* : the matrix form of the gate, both np and sparse should work
        """
        assert (matrix.shape[0] == matrix.shape[1]), "Only Square matrices allowed"
        assert (2**qreg_size == matrix.shape[0]), "Matrix dimensions and intended qreg size do not match"
        super(MGate, self).__init__(qreg_size)
        self.matrix = matrix 

    def __mul__(self, other):
        return MGate(kron(self.matrix, other.matrix), self.qreg_size + other.qreg_size)

    def __pow__(self, i):
        if i==0:
            return MGate(np.array([[1]]),0)
        elif i==1:
            return self
        else:
            return self*self**(i-1)

    def apply(self, qreg):
        new_state = self.matrix.dot(qreg.state)
        return QReg(len(qreg), new_state)


########################
# Common matrix gates  #

class PShiftGate(MGate):
    """Phase Shift Gate (on single qbit) 
    """
    def __init__(self, phi):
        """initalize phase shift gate for given phi
        *phi* : angle of shift (in radians)
        """
        matrix = np.array([
            [1,0],
            [0,np.exp(phi*1j)]
        ])
        super(PShiftGate, self).__init__(matrix,1)

class CPSGate(MGate):
    def __init__(self, phi):
        matrix = np.array([
            [1,0,0,0],
            [0,1,0,0],
            [0,0,1,0],
            [0,0,0,np.exp(phi*1j)]
        ])
        super(CPSGate, self).__init__(matrix, 2)

I = MGate(np.eye(2), 1)
H = MGate( 1 / np.sqrt(2) * np.array( [ [1,1] , [1,-1]] ), 1)

NOT = MGate(np.array([
    [0,1],
    [1,0]
]), 1)
CNOT = MGate(np.array([
    [1,0,0,0],
    [0,1,0,0],
    [0,0,0,1],
    [0,0,1,0]
]), 2)

CNOT_r = MGate(np.array([
    [1,0,0,0],[0,0,0,1],[0,0,1,0],[0,1,0,0]
]), 2)

cc_array = np.eye(8)
cc_array[-2:,-2:] = np.array([[0,1],[1,0]])
CCNOT = MGate(np.array(cc_array),3)

SWAP = MGate(np.array([
    [1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]
]), 2)

###########################

##############
# Lazy gates #

class LazyGate(Gate):
    """A gate that will act lazily on register rather than whole. User specifies matrix form of gate that will act on subqubits, other bits will be regarded as spectators. It does so by shuffling indices appropriately and multiplying the matrix in chunks by sliding down. Here is representation of how it works
    initial state --> sort indices according to values of spectator bits(!) and group the same values -> sort each `sub`register according to values of bits to be acted on(!!) (so it becomes like 00, 01, 10, 11) -> multiply matrix with each subgroup -> put all new values to their original position

    !: to extract spectator bits: apply bitwise a | ((base & (1<<q)) >> (q-i)), where q is the position of spectator, and put them into ith position  (i=0;i++)
    !!: similar to spectator bits, a | ((base & (1<<q)) >> q) << len(self.qbpos)-i-1, last shift is done so that first bit goes to leftmost rather than rightmost position
    """

    def __init__(self, matrix, qbpos, qreg_size):
        """create a lazily applied gate
        *matrix* : matrix form of base gate
        *qbpos* : the position of qbits that gate will apply to
        *qreg_size* : size of register it will act on
        """
        super().__init__(qreg_size)
        
        self.matrix = matrix
        self.msize = 2**len(qbpos)
        self.qsize = 2**qreg_size
        self.bases = np.arange(self.qsize)

        self.qbpos = qbpos
        self.specpos = list(filter(lambda x : not x in qbpos, range(qreg_size)))

        # time to shuffle indices
        # first group the indices which has same spectators
        self.bases = self.bases[np.argsort(self.gather_spectators())]
        # next for each group sort the necessary bits in order
        gbits = self.gather_qbpos()
        for i in range(0, self.qsize,self.msize):
            self.bases[i:i+self.msize] = self.bases[i:i+self.msize][np.argsort(gbits[i:i+self.msize])]

    def apply(self, qreg):
        new_state=np.zeros(len(qreg.state), dtype=np.complex)
        for i in range(0, self.qsize,self.msize):
            new_state[self.bases[i:i+self.msize]] = self.matrix.dot(qreg.state[self.bases[i:i+self.msize]])

        return QReg(len(qreg), new_state)

    def gather_qbpos(self):
        return np.array([self.get_qbpos(b) for b in self.bases])

    def gather_spectators(self):
        return np.array([self.get_spectator(b) for b in self.bases])

    def get_qbpos(self, base):
        a = 0
        for i in range(len(self.qbpos)):
            q = self.qbpos[i]
            a = a | ((base & (1<<q)) >> q) << len(self.qbpos)-i-1
        return a
    
    def get_spectator(self, base):
        a = 0
        for i in range(len(self.specpos)):
            q = self.specpos[i]
            a = a | ((base & (1<<q)) >> (q-i))
        return a

#########################################################################################
# Lazy (functional rather than matrix, not related to LazyGate) versions of commongates #
# !NOTE these are not complete and should not be used

class LazyCNOT:
    """Lazy controlled gate, in theory it should be possible to use many bits as control
    !WARNING: This class is incomplete do not use!!!
    """

    def __init__(self, ncontrols):
        self.ncontrols = ncontrols
    
    def dot(self, qreg_state):
        new_state = qreg_state.copy()
        new_state[-1], new_state[-2] = new_state[-2], new_state[-1]
        return new_state

class LazyNOT:
    """Lazy NOT gate. Should be able to reverese all bits at once very efficiently
    !WARNING: This class is incomplete do not use!!!
    """

    def __init__(self, nbits):
        self.nbits = nbits

    def dot(self, qreg_state):
        return qreg_state[::-1]
