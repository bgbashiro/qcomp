from qcomp.gate import CNOT, apply_gate_atbits
from qcomp.qregister import mk_reg
from qcomp.qbit import ZERO, ONE

q010101 = mk_reg([ZERO,ONE,ONE,ZERO,ZERO,ONE])
q0001 = apply_gate_atbits(CNOT, q010101, [3,5])
print(q0001)