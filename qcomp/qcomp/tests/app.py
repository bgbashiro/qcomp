from qcomp.gate import CNOT, apply_gate_atbits
from qcomp.qregister import mk_reg
from qcomp.qbit import ZERO, ONE

q010101 = mk_reg([ZERO,ZERO,ONE,ONE])
q0001 = apply_gate_atbits(CNOT, q010101, [2,3])
print(q0001)