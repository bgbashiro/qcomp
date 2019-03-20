from .gate import H, I, construct_unitary_F
from .qregister import mk_reg
from .qbit import PLUS, MINUS
import numpy as np

class Grover:
    """Wrapper around gate and register modules that makes it easy to run grover algorithms. 
    """

    def __init__(self, nbits):
        """Create instance of Grover algorithm
        *nbits* : number of bits in input space
        """
        self.nbits = nbits
        self.HS = H**nbits * I
        self.U0 = construct_unitary_F(["0"*nbits])
        self.oracle = None

    def def_oracle(self, fdef):
        """Add oracle to the algorithm
        *fdef* : str, the input that will map to 1 
        """
        self.oracle = construct_unitary_F([fdef])
        self.circuit = self.oracle + self.HS + self.U0 + self.HS

    def run_iteration(self, n_iterations=None):
        """Run the algorithm
        *n_iterations* : number of iterations the circuit will be applied (defaule 2^sqrt(nbits))

        Returns:
        *qreg* : QReg, final form of register
        """
        if n_iterations is None:
            n_iterations = int(np.sqrt(2**self.nbits))
        register = mk_reg([PLUS]*self.nbits + [MINUS])
        for _ in range(n_iterations):
            register = self.circuit.apply(register)

        return register