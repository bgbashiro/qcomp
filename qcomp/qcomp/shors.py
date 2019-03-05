import numpy as np
import random
import math as M



class shors(object):
    def __init__(self,N):
        self.N = N
        self.Q = self.N**2
        self.numQBits = M.log(self.Q)/M.log(2)
        self.seqList = []


    def func(self):
        while True:
            self.a = random.randint(1,self.N)
            self.gcd = M.gcd(self.a,self.N)
            print(M.gcd(self.a,self.N))
            if M.gcd(self.a,self.N) == 1:
                self.QFT()
                if self.period % 2:
                    break
                else:
                    self.root1 = M.gcd(int(self.a**(self.period/2)+1),self.N)
                    self.root2 = M.gcd(int(self.a**(self.period/2)-1),self.N)
                    break
            else:
                self.root3 = self.N/self.gcd


    def generateSequence(self):
        for i in range(0,self.Q):
            self.seqList.append([i, self.a**i % self.N])


    def QRT(self):
        self.generateSequence()
        #Build qBIT register with self.numQbits
        #Apply Hadamard Gate to the entire register.
        #Apply QFT gate.
        reg = q.mk_register(self.numQBits,zero)
        reg.apply(Hadamard)
        reg.apply(QFT)

        #self.period = 3515

def main():
    A = shors(6)
    A.func()
    print(A.a)

    print(A.seqList)
main()
