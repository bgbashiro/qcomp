"""gate module exposes generic Gate interface and various implementations. Main idea is the instances of these objects can be joined together in different forms using HChain and VChain (super)classes to form network and then single apply method can be called for the algorithm. 

HChain joins gates in the flow of algorithm (?), so each gate must use same number of qbits. E.g two 2qbit gates A and be chained together
<pre>
 \* - |   | - |   | - \*
     | A |   | B |
 \* - |   | - |   | - \*
</pre>
VChain joins gates in the direction of the register (?). Thus, resulsting gate will always need register of size equal to sum of each individual gate making the register. E.G 1qbit C and 2qbit D gates forming a complex E gate. 
<pre>

\*   |- | C | -|- \*
        ---
\*   |- |   | -|- \*\*
    |  | D |  |
\*   |- |   | -|- \*\*
        ---  
       **E**
</pre>
This can be used to represent spectator bits by using identity gate. By making corresponding gate ID, those bits will not influenced. For example, applying hadamard gate to 1st and 3rd qubits of |011> state
::

    #!python
    >>> from qcomp.qbit import ZERO,ONE
    >>> from qcomp.qregister import mk_reg
    >>> from qcomp.gate import VChain, HChain, H, I
    >>> qr = mk_reg([ZERO,ONE,ONE])
    >>> qr
    0.00+0.00j       |000> 
    0.00+0.00j       |001> 
    0.00+0.00j       |010> 
    1.00+0.00j       |011> 
    0.00+0.00j       |100> 
    0.00+0.00j       |101> 
    0.00+0.00j       |110> 
    0.00+0.00j       |111> 

    >>> HIH = VChain([H,I,H])
    >>> qr_t = HIH.apply(qr)
    >>> qr_t
    0.00+0.00j       |000> 
    0.00+0.00j       |001> 
    0.50+0.00j       |010> 
    -0.50+0.00j      |011> 
    0.00+0.00j       |100> 
    0.00+0.00j       |101> 
    0.50+0.00j       |110> 
    -0.50+0.00j      |111> 
"""
import numpy as np
from .qbit import *
from .qregister import *
from .utils import *
from functools import reduce

def construct_unitary_F(fdef, in_size=None):
    """Given fdef construct appropriate unitary transformation matrix, probably useful for grovers algorithm
    Parameters
    ---
    fdef: list of inputs that will calculate to 1. E.g. ["00","10"] -> will prepare gate that operates on 3 qbits. If first 2 qbits are 00 or 10 will AND 3rd one is 0, it will turn last one to 1. For making it unitary, if last bit is 1 it will turn into 0 for same inputs.
    in_size: expected input size, only necessary when fdef is empty list, ie function always produces 0
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

class MGate(Gate):
    """Gate that is formed by providing matrix form
    instance variables
    -----------
    *matrix* : the matrix form of the gate, both np and sparse should work

    methods
    -----------------
    apply
    """

    def __init__(self, matrix, qreg_size):
        assert (matrix.shape[0] == matrix.shape[1]), "Only Square matrices allowed"
        assert (2**qreg_size == matrix.shape[0]), "Matrix dimensions and intended qreg size do not match"
        super(MGate, self).__init__(qreg_size)
        self.matrix = matrix 

    def __add__(self, other):
        return MGate(np.matmul(self.matrix, other.matrix),self.qreg_size)

    def __mul__(self, other):
        return MGate(np.kron(self.matrix, other.matrix), self.qreg_size + other.qreg_size)

    def __pow__(self, i):
        if i==1:
            return self
        else:
            return self*self**(i-1)

    def apply(self, qreg):
        new_state = self.matrix.dot(qreg.state)
        return QReg(len(qreg), new_state)

class PShiftGate(MGate):
    """Phase Shift Gate (on single qbit)
    """
    def __init__(self, phi):
        matrix = np.array([
            [1,0],
            [0,np.exp(phi*1j)]
        ])
        super(PShiftGate, self).__init__(matrix,1)

    """Form a single gate (that will be applied to register in one step) 
    by joining different matrices. It is **different** from gate chain. Example use case:
    to apply hadamard gate to 3rd bit |000> we can form complex gate of form
    [ID, ID, Hadamard] 
    """

I = MGate(np.eye(2), 1)
H = MGate( 1 / np.sqrt(2) * np.array( [ [1,1] , [1,-1]] ), 1)

# not gates
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

def swap_bits(qreg,q0,q1):
    if q0 > q1:
        q0,q1=q1,q0
    t_qreg = qreg.copy()
    for c_pointer in range(q0,q1):
        glist = [I] * (c_pointer) + [SWAP] + [I] * (qreg.nbits-c_pointer-2)
        big_g = reduce(lambda g1,g2:g1*g2, glist)
        t_qreg = big_g.apply(t_qreg)
    for c_pointer in range(q0,q1)[:-1:-1]:
        glist = [I] * c_pointer + [SWAP] + [I] * (qreg.nbits-c_pointer-2)
        big_g = reduce(lambda g1,g2:g1*g2, glist)
        t_qreg = big_g.apply(t_qreg)
    return t_qreg

def apply_gate_atbits(gate, qreg, target_bits):
    
    spectator_bits_i = 0
    for b in target_bits:
        spectator_bits_i = spectator_bits_i | (1<<b)
    print(bin(spectator_bits_i))    
    spectator_bits_i = ((1<< qreg.nbits)-1) - spectator_bits_i
    print(bin(spectator_bits_i))
    # create base states
    bases = np.arange(0,2**qreg.nbits)
    
    # extract spectator bits
    specbits = np.argsort(bases & spectator_bits_i, kind='mergesort')
    print(specbits)
    new_qreg = qreg.copy()
    for i in range(0,2**new_qreg.nbits, 2**len(target_bits)):
        new_qreg.state[specbits[i:i+2**len(target_bits)]] = gate.matrix.dot(new_qreg.state[specbits[i:i+2**len(target_bits)]])
    return new_qreg

