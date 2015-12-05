from WLObjects import IsingGrid, PottsGrid, SpinGlass   #diffent models which are Wang-Landaugarithmable
import MetropolisModule
import sys

fileOut = 'outputMetropolis'
dim, gridSize = 2, 10
thermalizationSweeps, measurementSweeps = 3000, 3000
numAverages = 3
plotRange = [2.0,2.3,50]
gridMode = 0
grid = None

for i in range(0,len(sys.argv)):
    if(sys.argv[i] == '-o'):
        fileOut=sys.argv[i+1]
    if(sys.argv[i] == '-d'):
        dim=int(sys.argv[i+1])
    if(sys.argv[i] == '-n'):
        gridSize=int(sys.argv[i+1])
    if(sys.argv[i] == '-t'):
        thermalizationSweeps=int(sys.argv[i+1])
    if(sys.argv[i] == '-r'):
        plotRange[0]=float(sys.argv[i+1])
        plotRange[1]=float(sys.argv[i+2])
        plotRange[2]=int(sys.argv[i+3])
    if(sys.argv[i] == '-m'):
        measurementSweeps=int(sys.argv[i+1])
    if(sys.argv[i] == '-a'):
        numAverages=int(sys.argv[i+1])
    if(sys.argv[i] == '--mode'):
        gridMode=int(sys.argv[i+1])
    if(sys.argv[i] == '--help'):
        print "\n-o\t[Filename]\tOutput file name\n-d\t[dimension]\tDimensions of problem.\n-n\t[Size]\t\tSize of grid.\nIteration stop point order of magnitude\n-m\t[sweeps]\tMeasurement Sweeps\n-t\t[sweeps]\tThermalization Sweeps"
        sys.exit()
    
if(gridMode == 0):
    grid = IsingGrid(dim,gridSize)
if(gridMode == 1):
    grid = PottsGrid(dim,gridSize,10)
if(gridMode == 2):
    grid = SpinGlass(dim,gridSize)
else:
    print "Unrecognized grid mode"
    sys.exit()
    
print '-------------------------------------------------\Running Metropolis Algoritn\n' + 'Thermalization Sweeps:' + str(thermalizationSweeps) + '\tMeasurement Sweeps:'+ str(measurementSweeps)
MetropolisModule.runAlgorithm(grid, plotRange, thermalizationSweeps, measurementSweeps, numAverages, fileOut)