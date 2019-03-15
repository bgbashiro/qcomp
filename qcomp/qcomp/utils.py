"""
This is module that contains mathematical manipulations needed for inner workings of Quantum Computer. 
* SparseMatrix for high performance linear algebra operations
* !Not written yet, unified matrix API, which can utilize np arrays or sparse matrices under the hood
"""

import numpy as np
import copy

def kron(b,a):
    c = np.empty([b.shape[0]*a.shape[0], b.shape[1]*a.shape[1]])
    for i in range(b.shape[0]):
        for j in range(b.shape[1]):
            c[i*a.shape[0]:(i+1)*a.shape[0],j*a.shape[1]:(j+1)*a.shape[1]] = b[i,j]*a
    return c


### Sparse Matrices

class element(object):
    def __init__(self,rowIndex,colIndex,value):
        self.rowIndex = rowIndex
        self.colIndex = colIndex
        self.value = value

class SparseMatrix(object):
    def __init__(self,row, col):
        self.rowSize = row
        self.colSize = col
        self.elements = []

    def convertToMatrix(self):
        normalMatrix = np.zeros([self.rowSize,self.colSize], dtype=np.complex)
        for n in self.elements:
            normalMatrix[n.rowIndex][n.colIndex] = n.value
        return normalMatrix

    def dot(self, vector):
        svector = SparseMatrix.convertToSparse(vector)
        sans = self.multiplySparseMatrix(svector)
        return sans.convertToMatrix()

    def multiplySparseMatrix(self,that):
        elements = []
        assert self.colSize == that.rowSize,"Matrix Sizes Are Incompatible!!!"
        for i in range(0,self.rowSize):
            for j in range(0,that.colSize):
                values = []
                for n in self.elements:
                    for m in that.elements:
                        if n.rowIndex == i and m.colIndex == j and n.colIndex == m.rowIndex:
                            values.append(n.value*m.value)
                            break
                            #if this happens we can skip to the next n in self.elements
                    value = sum(values)
                    elements.append(element(i,j,value))

        new = SparseMatrix(self.rowSize,that.colSize)
        new.elements = elements
        return new
    
    @staticmethod
    def convertToSparse(matrix):
        try:
            _ = matrix.shape[1]
        except IndexError:
            matrix = matrix.reshape([-1,1])
        new = SparseMatrix(matrix.shape[0],matrix.shape[1])        
        for i in range(0,matrix.shape[0]):
            for j in range(0,matrix.shape[1]):
                if matrix[i][j] != 0:
                    new.elements.append(element(i,j,matrix[i][j]))
        return new

    def scalarMultiply(self,scalar):
        for i in range(0,len(self.elements)):
            self.elements[i].value = scalar*self.elements[i].value
    @staticmethod
    def kronSparse(this,that):
        new1 = []
        for element in this.elements:
            new = copy.deepcopy(that)
            for i in range(0,len(new.elements)):
                new.elements[i].value = element.value*new.elements[i].value
            cornerRow = element.rowIndex*that.rowSize
            cornerCol = element.colIndex*that.colSize
            for i in range (0,len(new.elements)):
                new.elements[i].rowIndex += cornerRow
                new.elements[i].colIndex += cornerCol
            for element in new.elements:
                new1.append(element)
        kronProd = SparseMatrix(this.rowSize*that.rowSize, this.colSize*that.colSize)
        kronProd.elements = new1
        return kronProd
