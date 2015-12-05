import numpy as np
import random, sys

from WLObjects import IsingGrid, PottsGrid   #diffent models which are Wang-Landaugarithmable
from Plot import plotArrayToFile
from fitCurve import fitParabola

dim = 2
periodicBoundaries=True
graphProgress = False
sizeN=10
fExp=0.1**7
minIterations = 30

fileOut = 'output'

for i in range(0,len(sys.argv)):
    if(sys.argv[i] == '-o'):
        fileOut=sys.argv[i+1]
    if(sys.argv[i] == '-d'):
        dim=int(sys.argv[i+1])
    if(sys.argv[i] == '-n'):
        sizeN=int(sys.argv[i+1])
    if(sys.argv[i] == '-a'):
        fExp=0.1**(float(sys.argv[i+1]))
    if(sys.argv[i] == '-g'):
        graphProgress = True
    if(sys.argv[i] == '--help'):
        print "\n-o\t[Filename]\tOutput file name\n-d\t[dimension]\tDimensions of problem.\n-n\t[Size]\t\tSize of grid.\n-a\t[accuracy]\tAccuracy control parameter. Iteration stop point order of magnitude\n-g\t\t\tShow graph being built\n"
        sys.exit()

if(graphProgress):
    import matplotlib; matplotlib.use('TKAgg'); import matplotlib.pyplot as plt
    plt.ion()

#Wang-Landau algorithm. Meat of the program...
def wangLandau(grid):
    logG = [0.0]*(grid.energyIndex(grid.MAX_E)+1)
    fFactor = np.e
    energy = grid.calcEnergy()
    while fFactor > np.exp(fExp):
        h = [0]*len(logG)
        count = 0
        while True:
            for iterations in range(0,grid.numSites):
                randIndex=[0.0]*grid.dimension
                for i in range(0,grid.dimension):
                    randIndex[i]=int(random.random()*grid.gridSize)
                tempVal = grid.getIndex(randIndex)
                enew = energy - grid.neighborEnergy(randIndex)*2
                if(np.abs(enew) > grid.MAX_E):
                    print 'Energy error\tE = ' + str(energy)
                    energy = grid.calcEnergy()
                    print ' Eactual = '+ str(energy)
                    enew = energy - grid.neighborEnergy(randIndex)*2
                    continue
                gNew, gOld = logG[grid.energyIndex(enew)], logG[grid.energyIndex(energy)]
                gRatioLog = gOld- gNew
                if( (gNew < gOld) | (np.log(random.random()) < (gRatioLog)) ):
                    logG[grid.energyIndex(enew)]+=np.log(fFactor)
                    h[grid.energyIndex(enew)]+=1
                    energy = enew
                    grid.flipIndex(randIndex)
                else:
                    logG[grid.energyIndex(energy)]+=np.log(fFactor)
                    h[grid.energyIndex(energy)]+=1
            if(graphProgress):
                plt.figure(1)
                plt.cla()
                plt.ylabel('Entropy (Log[state Density])')
                plt.plot(logG)
                plt.figure(2)
                plt.cla()
                plt.ylabel('Current Histogram')
                plt.plot(h)
                plt.draw()
            #if(isFlatQuad(h,.01**fExp) & (isFlatSlope(h,1.0**fExp)) & (count > 10) ):
            #if(isFlatStd(h,np.sqrt(np.log(fFactor))+fExp) & ((count > minIterations) | (0 not in h) ) ):
            if(isFlatAverage(h) & ((count > minIterations) | (0 not in h) ) ):
                minG=min(logG)
                for i in range(0,len(logG)):
                    logG[i]-=minG       #normalize logG vector to prevent overflows of exp()
                break
            count+=1
        fFactor=np.sqrt(fFactor)
        print fFactor
    out = [0.0]*len(logG)
    baseline = min(logG)
    for i in range(0,len(logG)):
        logG[i]-=baseline
        out[i]=np.exp(logG[i])
    #normConst = pow(float(dim**N)/sum(out),1.0/sumLogG)
    normConst = float(dim**grid.gridSize)/sum(out)
    normConstLog = np.log(normConst)
    #for i in range(0,len(logG)):
        #out[i]*=normConst
        #logG[i]+=normConstLog
       # pow(out[i],sumLogG)
    return [logG, out]

def isFlatAverage(histogram):
    avg = np.mean(histogram)
    for i in range(0,len(histogram)):
        if(histogram[i]<0.8*avg):
            return False
    return True
def isFlatStd(histogram, tolerance):
    avg = np.mean(histogram)
    std = np.std(histogram)
    return (std < (avg*tolerance))

def closeGraphs():
    plt.figure(1)
    plt.close()
    plt.figure(2)
    plt.close()

#main loop of program
grd = IsingGrid(dim,sizeN)

output = wangLandau(grd)
plotArrayToFile(output[0],fileOut)