import numpy as np
import random

from WLObjects import IsingGrid, PottsGrid, SpinGlass   #diffent models which are Wang-Landaugarithmable
from Plot import plotToFile, averageFunc

fileOut = 'outputMetropolis'
thermalizationSweeps, measurementSweeps = 3000, 3000
numAverages = 3
plotRange = [2.0,2.3,50]
grid = None
    
def metropolis(grid,numSteps,temperature):
    N = grid.gridSize
    dim=grid.dimension
    energy = grid.calcEnergy()
    magnetization = sum(grid.array)
    
    outE=[energy]*(numSteps+1)
    outESqr=[energy**2]*(numSteps+1)
    outMSqr=[magnetization**2]*(numSteps + 1)
    outMQuad=[magnetization**4]*(numSteps + 1)
    
    def transitionProb(time,temp,deltaH,spinI):
        gamma = 1.0#/(time+1)
        if(np.sign(spinI) == np.sign(deltaH)):
            gamma*=np.exp( -deltaH/temp )
        return gamma
    
    for iteration in range(0,numSteps):
        for i in range(0,N**dim):
            randPos=[0.0]*dim
            for j in range(0,dim):
                randPos[j]=int(random.random()*N)
            delta = -2*grid.neighborEnergy(randPos)
            sI=grid.getIndex(randPos)
            pTransition = transitionProb(iteration,temperature,delta,sI)
            if( pTransition > random.random()):   #accept change
                grid.flipIndex(randPos)
                energy += delta
                magnetization -= sI*2
                if(i%(N-1)==0):
                    outE[iteration+1] = energy
                    outESqr[iteration+1] = energy**2
                    outMSqr[iteration+1] = magnetization**2
                    outMQuad[iteration+1] = magnetization**4
                  #  print str(magnetization) + "\tM=" + str(sum(grid.array))
    return [outE, outESqr, outMSqr, outMQuad, grid]

###Define functions to pass to averaging and plotting modules
def mainLoop(temperature):
    #iterate pre-defined number of times
    global grid
    grid  = (metropolis(grid,thermalizationSweeps,temperature))[4]
    #then measure over a pre-defined number of sweeps
    output = metropolis(grid,measurementSweeps,temperature)
    numParticles=grid.numSites
    return [float(sum(output[0]))/len(output[0])/numParticles, float(sum(output[1]))/len(output[1])/numParticles**2, float(sum(output[2]))/len(output[2])/numParticles**2, float(sum(output[3]))/len(output[3])/numParticles**4]


def plotFunc(temperature):
    return averageFunc(mainLoop,numAverages,temperature) #average 3 times....

###execute everything and write file....
def runAlgorithm(model, temperatureRange,tSweeps,mSweeps,averages,outputFile):
    global grid, fileOut, thermalizationSweeps, measurementSweeps, numAverages, plotRange
    grid = model
    fileOut = outputFile
    thermalizationSweeps, measurementSweeps = tSweeps, mSweeps
    numAverages = averages
    plotRange = temperatureRange
    
    header = '#Thermalization Sweeps:' + str(thermalizationSweeps) + '\tMeasurement Sweeps:'+ str(measurementSweeps) +'\n#' + str(model.dimension) + 'dimensional grid with N=' + str(model.gridSize) + '\n#Pairs of y values (yval, yerr)\n#Temp\tEnergy\tE^2\tMagnetization^2\tMagnetization^4\n'
    
    plotToFile(plotFunc,plotRange,fileOut,header,False,True)


