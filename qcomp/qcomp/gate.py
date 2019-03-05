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

class HChain(Gate):
    """Chain of gates that will be applied one by one
    """

    def __init__(self, gates):
        qsizes = [g.qreg_size for g in gates]
        assert all([qsizes[0] == q for q in qsizes]), "All gates should have same qregsize, use ID gate for spectators"
        self.gates = gates
        super(HChain, self).__init__(qsizes[0])

    def apply(self, qreg):
        # copy the input register
        temp_reg = qreg.copy()

        # apply gates sequentally
        for gate in self.gates:
            temp_reg = gate.apply(temp_reg)
        return temp_reg

class VChain(MGate):
    """Form a single gate (that will be applied to register in one step) 
    by joining different matrices. It is **different** from gate chain. Example use case:
    to apply hadamard gate to 3rd bit |000> we can form complex gate of form
    [ID, ID, Hadamard] 
    """

    def __init__(self, gates):
        """Form the complex gate. Note that new qregister size will be sum of
        qregsize of all gates AND gates provided MUST have matrix form (instance of MGate)
        
        Parameters
        -------------
        **gates**: list of gates, they will be applied from top to bottom
        """
        matrix_form = np.array([[1]], dtype=np.complex64)
        qreg_size = 0
        for gate in gates:
            matrix_form = np.kron(matrix_form, gate.matrix)
            qreg_size += gate.qreg_size
        super(VChain, self).__init__(matrix_form, qreg_size)

class SGate(Gate):
    """Swapping Gate. Normally, matrix form of a gate which acts on 2qbits(usually using one of them as control), will apply to neighbouring bits. This one utilizes swapping to work on desired bits and swaps them back. Because SGates do not have matrix form they cannot be VChained
    """
    def __init__(self, matrix, cbit, abit, qreg_size):
        """Make a gate that has defined matrix form, which will entangle bits specified rather than neighbours.

        Parameters:
        matrix.: matrix form of gate (4x4 expected as it will act on 2 bits)
        control_bit: first bit it acts on (usually used as control)
        act_bit: second bit it acts on (usually bit that will change when control bit is 0)
        """
        self.inner_gate = VChain([MGate(matrix,2)] + [I]*(qreg_size-2))
        self.b0 = cbit
        self.b1 = abit
        super(SGate, self).__init__(qreg_size)

    def apply(self, qreg):
        temp_reg = qreg.copy()
        # move control bit to 0, act bit to 1
        temp_reg.swap(0,self.b0)
        temp_reg.swap(1,self.b1)
        temp_reg = self.inner_gate.apply(temp_reg)
        # swap moved bits back
        temp_reg.swap(self.b0, 0)
        temp_reg.swap(self.b1, 1)
        return temp_reg

class CGate(SGate):
    """Controlled gate. This one will act on the second bit (in some predefined way) when control bit is 0. Otherwise it stays intact. Controlled gates have matrix in the from
    [1 0 0  0]
    [0 1 0 0]
    [0 0 a b]
    [0 0 c d]
    where a,b,c,d are transformation matrix that will be applied (when cbit is 1). These are SGate, so cannot be V-Chained
    """
    def __init__(self, transformation_matrix, cbit=0, abit=1, qreg_size=2):
        M = np.eye(4, dtype=np.complex)
        M[-2:,-2:] = transformation_matrix
        super(CGate, self).__init__(M,cbit,abit,qreg_size)

class CPSGate(CGate):
    """Controlled Phase Shift Gate
    """
    def __init__(self, phi, cbit=0, abit=1, qreg_size=2):
        """Parameters
        phi: amount of shift
        cbit: control bit
        abit: bit that will it act on
        qreg_size: size of Quantum register
        """
        matrix = PShiftGate(phi).matrix
        super(CPSGate, self).__init__(matrix, cbit, abit, qreg_size)

class CCGate(Gate):
    """Similar to CGate except two bits are used for control. When both is one transformation is applied, otherwise nothing is done
    """
    def __init__(self, transformation_matrix, cbits=[0,1], abit=2, qreg_size=3):
        """by defualt prepares 3 qbit gate that uses first 2 as control, acts on 3rd
        """
        self.c0 = cbits[0]
        self.c1 = cbits[1]
        self.a = abit
        M = np.eye(2**3, dtype=np.complex)
        M[-2:,-2:] = transformation_matrix
        self.inner_gate = VChain([MGate(M,3)] + [I]*(qreg_size - 3))
        self.qreg_size = self.inner_gate.qreg_size

    def apply(self, qreg):
        temp_reg = qreg.copy()
        # swap control bits
        temp_reg.swap(0,self.c0)
        temp_reg.swap(1,self.c1)
        # move acting bit
        temp_reg.swap(2,self.a)
        temp_reg = self.inner_gate.apply(temp_reg)
        # put things back
        temp_reg.swap(0,self.c0)
        temp_reg.swap(1,self.c1)
        # move acting bit
        temp_reg.swap(2,self.a)
        return temp_reg

class MultiCGate(Gate):
    """Extension of CCGate to more than 2 qbits. Currently not working, needs gather and scatter to be implemented in the QReg class
    """
    def __init__(self, transformation_matrix, cbits, abit, qreg_size):
        self.c = cbits
        self.a = abit
        cmatrix_size = len(self.c)+1
        M = np.eye(cmatrix_size)
        M[-2:,-2:] = transformation_matrix
        self.inner_gate = VChain([MGate(M, cmatrix_size)] + [I]*(qreg_size - cmatrix_size))
    
    def apply(self, qreg):
        temp_qreg = qreg.copy()
        temp_qreg.gather(self.c + [self.a])
        self.inner_gate.apply(temp_qreg)
        temp_qreg.scatter()
        return temp_qreg

I = MGate(np.eye(2), 1)
H = MGate( 1 / np.sqrt(2) * np.array( [ [1,1] , [1,-1]] ), 1)

# not gates
not_array = np.array([
    [0,1],[1,0]
])
NOT = MGate(not_array, 1)

class CNOT(CGate):
    def __init__(self, cbit, abit, qreg_size):
        super(CNOT, self).__init__(not_array, cbit, abit, qreg_size)

class CCNOT(CCGate):
    def __init__(self,cbits, abit, qreg_size):
        super(CCNOT, self).__init__(not_array, cbits, abit, qreg_size)