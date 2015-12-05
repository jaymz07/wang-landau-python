###Module Containing Classes used for Wang-Landau Algorithm.
import random
import numpy as np
import bisect

###Ising model----------------------------------------------------
class IsingGrid:
    def __init__(self,dimension,N):
        self.dimension = dimension
        self.gridSize = N
        self.numSites = N**dimension
        self.MAX_E = dimension*N**dimension
        self.array = [0]*N**dimension
        self.fieldB = 0.0
        for i in xrange(0,N**dimension):
            self.array[i] = -1+2*int(random.random()*2)
    
    def getIndex(self,indexArray):
        index = 0;
        for i in range(0,self.dimension):
            index+=self.gridSize**(self.dimension - i - 1)*indexArray[self.dimension-i-1]
        return self.array[index]
    def flipIndex(self,indexArray):
        index = 0;
        for i in range(0,self.dimension):
            index+=self.gridSize**(self.dimension - i - 1)*indexArray[self.dimension-i-1]
        self.array[index]*=-1
    def neighborHalfEnergy(self,indexArray):
        sum=0.0
        val=self.getIndex(indexArray)
        dim = self.dimension
        N = self.gridSize
        for i in range(0,dim):
            ind = indexArray[0:i] + [ (indexArray[i] + 1) % N ] + indexArray[i+1:dim] #sum up each neigbor by adding one to each coordinate dimension
            sum+= -self.getIndex(ind)*val
        return sum
    def neighborEnergy(self,indexArray):
        sum=0.0
        val=self.getIndex(indexArray)
        dim = self.dimension
        N = self.gridSize
        for i in range(0,dim):
            ind1 = indexArray[0:i] + [ (indexArray[i] + 1) % N ] + indexArray[i+1:dim] #sum up each neigbor by adding one and subtracting one to each coordinate dimension
            ind2 = indexArray[0:i] + [ (indexArray[i] - 1) % N ] + indexArray[i+1:dim]
            sum+= -self.getIndex(ind1)*val
            sum+= -self.getIndex(ind2)*val
        return sum
    def calcEnergy(self,indexArrayPartial = []):
        if(len(indexArrayPartial)==self.dimension):
            return self.neighborHalfEnergy(indexArrayPartial) + self.getIndex(indexArrayPartial)*self.fieldB
        else:
            sum = 0.0
            for i in range(0,self.gridSize):
                sum += self.calcEnergy(indexArrayPartial+[i])
            return sum
    def deltaEnergy(self, indexArray):
        return -2*(self.neighborEnergy(indexArray) + self.getIndex(indexArray)*self.fieldB)
    def energyIndex(self,en):
        out = int((en+self.MAX_E)/4)
        ind = self.dimension
        sub = 0
        for i in range(0,self.dimension):
            dim = self.dimension - i
            if(out>=ind):
                sub+=dim-1
            ind+=dim-1
        for i in range(1,self.dimension):
            dim = i
            if(out>self.MAX_E/2-ind):
                sub+=i-1
            ind-=i
        return out - sub;

###Pots model-------------------------------
class PottsGrid:
    def __init__(self,dimension,N,numStates):
        self.numStates = numStates
        self.dimension = dimension
        self.gridSize = N
        self.numSites = N**dimension
        self.MAX_E = dimension*N**dimension
        self.array = [0]*N**dimension
        self.fieldB = 0.0
        for i in xrange(0,N**dimension):
            self.array[i] = int(random.random()*numStates)
    
    def getIndex(self,indexArray):
        index = 0;
        for i in range(0,self.dimension):
            index+=self.gridSize**(self.dimension - i - 1)*indexArray[self.dimension-i-1]
        return self.array[index]
    def setIndex(self,indexArray,val):
        index = 0;
        for i in range(0,self.dimension):
            index+=self.gridSize**(self.dimension - i - 1)*indexArray[self.dimension-i-1]
        self.array[index]=val
    def flipIndex(self,indexArray):
        index = 0;
        for i in range(0,self.dimension):
            index+=self.gridSize**(self.dimension - i - 1)*indexArray[self.dimension-i-1]
        val = self.array[index]
        self.array[index]=int(random.random()*self.numStates)
        return val
    def neighborHalfEnergy(self,indexArray):
        sum=0.0
        val=self.getIndex(indexArray)
        dim = self.dimension
        N = self.gridSize
        for i in range(0,dim):
            ind = indexArray[0:i] + [ (indexArray[i] + 1) % N ] + indexArray[i+1:dim] #sum up each neigbor by adding one to each coordinate dimension
            sum+= -np.cos(self.getIndex(ind)*np.pi*2/self.numStates-val*np.pi*2/self.numStates)
        return sum
    def neighborEnergy(self,indexArray):
        sum=0.0
        val=self.getIndex(indexArray)
        dim = self.dimension
        N = self.gridSize
        for i in range(0,dim):
            ind1 = indexArray[0:i] + [ (indexArray[i] + 1) % N ] + indexArray[i+1:dim] #sum up each neigbor by adding one and subtracting one to each coordinate dimension
            ind2 = indexArray[0:i] + [ (indexArray[i] - 1) % N ] + indexArray[i+1:dim]
            sum+= -np.cos(self.getIndex(ind1)*np.pi*2/self.numStates-val*np.pi*2/self.numStates)
            sum+= -np.cos(self.getIndex(ind2)*np.pi*2/self.numStates-val*np.pi*2/self.numStates)
        return sum
    def calcEnergy(self,indexArrayPartial = []):
        if(len(indexArrayPartial)==self.dimension):
            return self.neighborHalfEnergy(indexArrayPartial)
        else:
            sum = 0.0
            for i in range(0,self.gridSize):
                sum += self.calcEnergy(indexArrayPartial+[i])
            return sum
    def energyIndex(self,en):
        return (self.MAX_E + int(en))/8
        
###Ising model----------------------------------------------------
class SpinGlass:
    def __init__(self,dimension,N):
        self.dimension = dimension
        self.gridSize = N
        self.numSites = N**dimension
        self.MAX_E = dimension*N**dimension
        self.array = [0]*N**dimension
        self.j = [0]*N**dimension
        self.fieldB = 0.0
        for i in xrange(0,N**dimension):
            self.array[i] = -1+2*int(random.random()*2)
            self.j[i] = -1+2*int(random.random()*2)
    
    def getIndex(self,indexArray):
        index = 0;
        for i in range(0,self.dimension):
            index+=self.gridSize**(self.dimension - i - 1)*indexArray[self.dimension-i-1]
        return self.array[index]
    def getJ(self,indexArray):
        index = 0;
        for i in range(0,self.dimension):
            index+=self.gridSize**(self.dimension - i - 1)*indexArray[self.dimension-i-1]
        return self.j[index]
    def flipIndex(self,indexArray):
        index = 0;
        for i in range(0,self.dimension):
            index+=self.gridSize**(self.dimension - i - 1)*indexArray[self.dimension-i-1]
        self.array[index]*=-1
    def neighborHalfEnergy(self,indexArray):
        sum=0.0
        val=self.getIndex(indexArray)
        jVal=self.getJ(indexArray)
        dim = self.dimension
        N = self.gridSize
        for i in range(0,dim):
            ind = indexArray[0:i] + [ (indexArray[i] + 1) % N ] + indexArray[i+1:dim] #sum up each neigbor by adding one to each coordinate dimension
            sum+= -self.getIndex(ind)*val*self.getJ(ind)*jVal
        return sum
    def neighborEnergy(self,indexArray):
        sum=0.0
        val=self.getIndex(indexArray)
        jVal=self.getJ(indexArray)
        dim = self.dimension
        N = self.gridSize
        for i in range(0,dim):
            ind1 = indexArray[0:i] + [ (indexArray[i] + 1) % N ] + indexArray[i+1:dim] #sum up each neigbor by adding one and subtracting one to each coordinate dimension
            ind2 = indexArray[0:i] + [ (indexArray[i] - 1) % N ] + indexArray[i+1:dim]
            sum+= -self.getIndex(ind1)*val*jVal*self.getJ(ind1)
            sum+= -self.getIndex(ind2)*val*jVal*self.getJ(ind2)
        return sum
    def calcEnergy(self,indexArrayPartial = []):
        if(len(indexArrayPartial)==self.dimension):
            return self.neighborHalfEnergy(indexArrayPartial)
        else:
            sum = 0.0
            for i in range(0,self.gridSize):
                sum += self.calcEnergy(indexArrayPartial+[i])
            return sum
    def deltaEnergy(self, indexArray):
        return -2*(self.neighborEnergy(indexArray) + self.getIndex(indexArray)*self.fieldB)
    def energyIndex(self,en):
        out = int((en+self.MAX_E)/4)
        ind = self.dimension
        sub = 0
        for i in range(0,self.dimension):
            dim = self.dimension - i
            if(out>=ind):
                sub+=dim-1
            ind+=dim-1
        for i in range(1,self.dimension):
            dim = i
            if(out>self.MAX_E/2-ind):
                sub+=i-1
            ind-=i
        return out- sub;
        
###Dynamic Energy histogram

class EnergyHistogram:
    def __init__(self):
        self.N=0
        self.h = []
        self.logG = []
        self.energies = []
        self.energyTolerance = 0.0000001
    def addValue(self,energy, fFactor):
        if(self.N==0):
            self.N=1
            self.h.append(1)
            self.energies.append(energy)
            self.logG.append(fFactor)
            return
        index = bisect.bisect_left(self.energies, energy)
        if((index<self.N) & (index >= 0)):
            if((np.abs(self.energies[index] - energy) < self.energyTolerance)):
                self.logG[index]+=fFactor
                self.h[index]+=1
                return
            if(index > 0):
                if((np.abs(self.energies[index - 1] - energy) < self.energyTolerance)):
                    self.logG[index - 1]+=fFactor
                    self.h[index - 1]+=1
                    return
        self.logG = self.logG[0:index] + [fFactor] + self.logG[index:self.N]
        self.energies = self.energies[0:index] + [energy] + self.energies[index:self.N]
        self.h = self.h[0:index] + [1] + self.h[index:self.N]
        self.N+=1
    def gValue(self,energy):
        for i in range(0,self.N):
            if(np.abs(energy - self.energies[i]) < self.energyTolerance):
                return self.logG[i]
        return 0.0
    def resetHistogram(self):
        self.h = [0]*self.N