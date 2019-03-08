from qcomp.algorithms import Grover
import numpy as np

def main():
    g = Grover(6)
    fdef = "".join(['0' if np.random.random()>0.5 else '1' for _ in range(6)])
    g.def_oracle(fdef)
    print("We are looking for " + fdef)
    reg = g.run_iteration()
    final = ""
    for i in range(reg.nbits-1):
        _, p = reg.get_qbit(i)
        if p.state[0]>0.5:
            final += '0'
        else:
            final += '1'
    for i in range(reg.nbits):
        print(reg.get_qbit(i)[0])
        
    print("Predicted : "+final)

main()

        