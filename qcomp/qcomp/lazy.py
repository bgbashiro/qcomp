import numpy as np

class LazyMatrix:
    
    def __init__(self, sm, smdim, qbpos):
        self.sm = sm
        self.smdim = smdim
        self.qbpos = qbpos
        
    
    def apply(self, v):
        w = v.copy()
        
        for i in range(len(w)):
            r=self.gather(i)
            i0=i & ~self.scatter(r)
            for c in range(self.smdim):
                j = i0|self.scatter(c)
                w[i] += self.sm[r,c]*v[j]
        return w

        

    def gather(self,i):
        j = 0
        for k in range(len(self.qbpos)):
            j = j | ((i>>self.qbpos[k]) & 1) << k

        return j
    
    def scatter(self, j):
        i = 0
        for k in range(len(self.qbpos)):
            i = i | ((j>>k)&1)<<self.qbpos[k]

        return i


def main():

    lm = LazyMatrix(np.array([
        [1,0,0,0],[0,1,0,0],[0,0,0,1], [0,0,1,0]
    ]), 4, [0,2])

    res = lm.gather(0)
    print(res)

main()