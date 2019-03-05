import numpy as np

class matrix(object):
    def __init__(self,size):
        self.matrix = np.zeros(5)
   self.matrix[n] = [1,2]
ValueError: setting an array element with a sequence.
Daniels-MacBook-Pro:qcomp Dan$ vim matrix.py

  1 import numpy as np
  2 
  3 class matrix(object):
  4     def __init__(self,size):
  5         self.matrix = np.zeros(5)
  6         for n in range(0,len(self.matrix)):
  7             self.matrix[n] = [1,2]
  8 
  9 def main():
 10     n = matrix(3)
 11     print(n.matrix)
 12 main()
~                                                                                                                                                                                                                                                                                     
~                                  for n in range(0,len(self.matrix)):
            self.matrix[n] = [1,2]

def main():
    n = matrix(3)
    print(n.matrix)
main()
