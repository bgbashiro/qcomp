"""qregister module exposes QReg object and extra functions needed for register manupulation
"""

import numpy as np
from numpy.linalg import norm

def mk_reg(qbits):
    """Using kron product construct a quantum register from single bits

    Parameters
    ----------------------
    
    qbits : List of QBit objects

    Returns
    -------------------------
    qreg :: constructed QReg object
    """

    state = np.array([[1]])
    for q in qbits:
        state = np.kron(state, q.state)
    return QReg(len(qbits), state=state.ravel())


class QReg():
    """Object representing Quantum Register
    
    Instance variables
    ------------------------
    nbits: number of qbits in the register
    state: array of coefficients of each basis state

    Methods
    --------------------
    normalize

    """

    def __init__(self, nbits, state=None):
        """
        Parameters
        ------------------
        nbits: number of qbits in the register
        state: default None, if not given it will create register at state. It will be normalized automatically, so advance check is not necessary
        |0000..0> (ie state = [1,0,0....] in coeffs form). alternatively state can be passed
        """

        self.nbits = nbits
        if state is None:
            self.state = np.zeros(2**nbits, dtype=np.complex64)
            self.state[0] = 1.0
        else:
            self.state = np.array(state, dtype=np.complex64)
        self.normalize()

        # construct a base representation, useful for few things, perhaps
        self.bases = [bin(i_base)[2:] for i_base in range(2**self.nbits)]
        append_zeros = lambda bits : "0"*(self.nbits-len(bits)) + bits
        self.bases = [append_zeros(b) for b in self.bases]

    def normalize(self):
        """Set the magnitude of state vector to 1, keeps overall probability equal to 1
        Does not return anything
        """
        self.state = self.state / norm(self.state)

    def swap(self, q0, q1):
        """Swap the two bits of system (q0th and q1th). E.g. 
        010 -> 100 swapped 0th and 1st bit

        Swapping is done through bitwise xor with "state" that has 1s on desired bits. To swap 0th and 1st bit we xor with 110 (last one is 0 because we dont change that). E.g.
        110 xor 010 -> (1 xor 0)(1 xor 1)(0 xor 0) -> 100
        Problem arises when both swapped bits are same (00 or 11). Like,
        110 xor 001 -> (1 xor 0)(1 xor 0)(0 xor 0) -> 110 (we needed first two bits stay as 0)
        So extra check is done by making sure they are not same by using bitwise and with each individual bit which is the shifted to LSB and compared. If input is 110
        100 and 110 -> 100 - shift to LSB -> 001
        010 and 110 -> 010 -shift to LSB -> 001
        if results are same that means no swapping should be done
        """
        shift0 = self.nbits - q0 -1
        shift1 = self.nbits - q1 -1
        swapper_1 = 1 << shift0
        swapper_2 = 1 << shift1
        swapper = swapper_1 | swapper_2
        swappedstate = np.zeros_like(self.state)
        for i in range(len(self.state)):
            ex1 = (swapper_1 & i) >> shift0
            ex2 = (swapper_2 & i) >> shift1
            if ex1 == ex2:
                swappedstate[i] = self.state[i]
            else:
                opp_i = swapper ^ i
                swappedstate[opp_i] = self.state[i]
        self.state = swappedstate

    def gather(self, qbits):
        """Swap the desired bits to then LSB slots (end of register).
        E.g. we have stat |11010>
        if I call gathe([0,1]) that means I want first twi bits (11) in the end. So we should return |01011>.
        There is swap method which swaps two bits at a time. This one probably use it in some form of loop
        """
        raise NotImplementedError("This shit aint implemented!")
    
    def scatter(self, qbits):
        """Reverse what gather has done
        """
        raise NotImplementedError("Do not try it bro!")

    def copy(self):
        """Return copy of register
        """
        return QReg(self.nbits, self.state)

    def __len__(self):
        return self.nbits

    def __repr__(self):
        """print state of qreg as a|00..0> + b|00..1> + ...
        """
        representation = ""
        for coeff, base in zip(self.state, self.bases):
            representation += "{0:.2f} \t |{1}> \n".format(coeff, base)
        return representation

    def get_qbit(self, i_qbit):
        """Get state of ith qbit
        """
        state_1_coeffs = self.state[[b[i_qbit] == '1' for b in self.bases]]
        state_1_proba = (state_1_coeffs**2).sum()
        state_0_coeffs = self.state[[b[i_qbit] == '0' for b in self.bases]]
        state_0_proba = (state_0_coeffs**2).sum()
        return "{}|0> + {}|1>".format(state_0_proba, state_1_proba)