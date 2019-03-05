"""qbit module contains class QBit, and constants
ZERO, ONE which represent 1-qbit |0> and |1> states
"""

import numpy as np
from numpy.linalg import norm


class QBit:
    """QBit is represented as 2D complex vector
    """

    def __init__(self, state):
        """
        Instance variables
        -----------------
        state: numpy array (size 2)

        Methods
        -----------------
        normalize
        """
        self.state = np.array(state, dtype=np.complex64)
        self.normalize()

    def normalize(self):
        """Set the magnitude of state vector to 1, keeps overall probability equal to 1
        Does not return anything
        """
        self.state = self.state / norm(self.state)

ZERO = QBit([1,0])
ONE = QBit([0,1])
MINUS = QBit([1,-1])
PLUS = QBit([1,1])
