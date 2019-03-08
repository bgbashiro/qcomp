from .gate import H, I, construct_unitary_F
from .qregister import mk_reg
from .qbit import ZERO, MINUS
import numpy as np

class Grover:

    def __init__(self, nbits):
        self.nbits = nbits
        self.HS = H**nbits * I
        self.U0 = construct_unitary_F(["0"*nbits])
        self.oracle = None

    def def_oracle(self, fdef):
        self.oracle = construct_unitary_F([fdef])
        self.circuit = self.oracle + self.HS + self.U0 + self.HS

    def run_iteration(self, n_iterations=None):
        if n_iterations is None:
            n_iterations = int(np.sqrt(2**self.nbits))
        register = self.HS.apply(mk_reg([ZERO]*self.nbits + [MINUS]))
        for _ in range(n_iterations):
            register = self.circuit.apply(register)

        return register