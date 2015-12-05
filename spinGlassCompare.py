##standard libraries
import numpy as np
import random, sys

##User defined libraries
from WLObjects import IsingGrid, PottsGrid, EnergyHistogram, SpinGlass
from Plot import plotArrayToFile
import MetropolisModule

##default parameters
dim = 3
graphProgress = False
verbose = True
gridMode = 2
sizeN=4
fExp=0.1**7
minIterations = 1000
magField = 0.0
flatnessTestPeriod = 100
fileOut = 'SPINGLASS_D2N4'

##argument handling
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
    if(sys.argv[i] == '-m'):
        gridMode=int(sys.argv[i+1])
    if(sys.argv[i] == '-f'):
        magField=float(sys.argv[i+1])
    if(sys.argv[i] == '-v'):
        verbose = True
    if(sys.argv[i] == '--help'):
        print "\n-o\t[Filename]\tOutput file name\n-d\t[dimension]\tDimensions of problem.\n-n\t[Size]\t\tSize of grid.\n-a\t[accuracy]\tAccuracy control parameter. Iteration stop point order of magnitude\n-b\t\t\tUse periodic boundary conditions\n-g\t\t\tShow graph being built\n-m\t[Mode#]\tChoose Which model to compute density of states (0==Ising, 1==Potts)"
        sys.exit()

##if graphing option set, import appropriate libraries
if(graphProgress):
    import matplotlib; matplotlib.use('TKAgg'); import matplotlib.pyplot as plt
    plt.ion()

##Wang-Landau algorithm. Meat of the program...
def wangLandau(grid):
    hist = EnergyHistogram()
    fFactor = np.e
    energy = grid.calcEnergy()
    fSweeps=0
    while fFactor > np.exp(fExp):
        hist.resetHistogram()
        count = 0
        while True:
            for iterations in range(0,grid.numSites):
                randIndex=[0.0]*grid.dimension
                for i in range(0,grid.dimension):
                    randIndex[i]=int(random.random()*grid.gridSize)
                enew, tempVal = energy, grid.getIndex(randIndex)
                if((gridMode==0) | (gridMode==2)):
                    enew = energy + grid.deltaEnergy(randIndex)
                elif(gridMode==1):
                    delta = -grid.neighborEnergy(randIndex)
                    tempVal = grid.flipIndex(randIndex)
                    delta += grid.neighborEnergy(randIndex)
                    enew = energy + delta
                gNew, gOld = hist.gValue(enew), hist.gValue(energy)
                gRatioLog = gOld- gNew
                if( (gNew < gOld) | (np.log(random.random()) < (gRatioLog)) ):
                    hist.addValue(enew,fFactor)
                    if((gridMode==0) | (gridMode==2)):
                        grid.flipIndex(randIndex)
                    energy = enew
                else:
                    hist.addValue(energy,fFactor)
                    if(gridMode==1):
                        grid.setIndex(randIndex,tempVal)
            if(graphProgress):
                plt.figure(1)
                plt.cla()
                plt.ylabel('Entropy (Log[state Density])')
                plt.plot(hist.energies,hist.logG)
                plt.figure(2)
                plt.cla()
                plt.ylabel('Current Histogram')
                plt.plot(hist.h)
                plt.draw()
            count+=1
            #if(isFlatStd(h,np.sqrt(np.log(fFactor))+fExp) & ((count > minIterations) | (0 not in h) ) ):
            if((count%flatnessTestPeriod)==0):
                if(isFlatAverage(hist.h) & (isFlatStd(hist.h,np.sqrt(np.log(fFactor))+fExp)) & ((count > minIterations) | (0 not in hist.h) ) ):
                    minG=min(hist.logG)
    #               for i in range(0,hist.N):
    #                    logG[i]-=minG       #normalize logG vector to prevent overflows of exp()
                    break
                if(verbose):
                    print "Iterations:" + str(fSweeps) + "\tSweeps:" + str(count) + "\tNumber of energy states:" + str(hist.N)
        fFactor=np.sqrt(fFactor)
        fSweeps+=1
        print "F="+str(fFactor)+"\tSweeps:"+str(count)+"\tNumber of Energy States counted:"+str(hist.N)
    out = [0.0]*hist.N
    baseline = min(hist.logG)
    for i in range(0,hist.N):
        hist.logG[i]-=baseline
        out[i]=np.exp(hist.logG[i])
        hist.energies[i]/=grid.numSites
    #normConst = pow(float(dim**N)/sum(out),1.0/sumLogG)
    normConst = float(dim**grid.gridSize)/sum(out)
    normConstLog = np.log(normConst)
    #for i in range(0,len(logG)):
        #out[i]*=normConst
        #logG[i]+=normConstLog
       # pow(out[i],sumLogG)
    return [hist.logG, out, hist.energies]

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

##Main execution of program

#define model grid
grd = IsingGrid(dim,sizeN)
if(gridMode == 1):
    grd = PottsGrid(dim,sizeN,10)
if(gridMode==2):
    grd= SpinGlass(dim,sizeN)

#set magnetic field
grd.fieldB = magField

#run....
output = wangLandau(grd)

#write to file!
dataOut =[]
for i in range(0,len(output[0])):
    dataOut.append([output[2][i],output[0][i]])
plotArrayToFile(dataOut,fileOut)

thermalizationSweeps, measurementSweeps = 10000, 50000
plotRange = [0.0, 2.0, 30]

MetropolisModule.runAlgorithm(grd, plotRange, thermalizationSweeps, measurementSweeps, 5, fileOut+"_Metropolis")
