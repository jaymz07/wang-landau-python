#module to automate averaging and saving data to file in format compatible with gnuplot
import numpy
def plotToFile(plotFunction,plotRange,fileOut,header='',log10plot=False,verbose=False):
    from subprocess import call
    file = open(fileOut, 'w')
    
    file.write("#" + header)
    plotVals = numpy.linspace( plotRange[0],plotRange[1], plotRange[2] )
    for k in range(0,len(plotVals)):
        x = plotVals[k]
        if(log10plot):
            x=10**plotVals[k]
        out=plotFunction(x)
        file.write( str(x) )
        if( isinstance(out,list) ):
            for i in range(0,len(out)):
                file.write('\t' + str(out[i]) )
        else:
            file.write('\t' + str(out) )
        file.write('\n')
        if(verbose):
            print str(float(k)*100.0/len(plotVals)) + " % done"
    file.close()

#used for plotting a matrix of data points (to a file for gnuplot)
def plotArrayToFile(plotData,fileOut,header='',writeIndex=True):
    from subprocess import call
    file = open(fileOut, 'w')
    
    file.write("#" + header + '\n')
    for k in range(0,len(plotData)):
        if(writeIndex):
            file.write(str(k) + '\t')
        if( isinstance(plotData[k],list)):
            for i in range(0,len(plotData[k])):
                if(i!=0):
                    file.write('\t')
                file.write(str(plotData[k][i]) + '\t')
        else:
            file.write(str(plotData[k]))
        file.write('\n')
    file.close()
    
#averages multiple runs of a function with the same arguments (good for monte carlo functions). Also produces standard deviations sqrt(<X^2> - <x>^2)
###Now compatible with functions that give arrays as output. Averages each element & produces deviations
def averageFunc(func,numAverages,funcArg):
    f=func(funcArg)
    if(not isinstance(f,list) ):
        sum, sumSqr = 0.0, 0.0
        for i in range(0,numAverages):
            sum += f
            sumSqr += f**2
            f=func(funcArg)
        return [sum/numAverages,numpy.sqrt(sumSqr/numAverages-(sum/numAverages)**2)]
    else:
        dim=len(f)
        sum, sumSqr = [0.0]*dim, [0.0]*dim
        for i in range(0,numAverages):
            for j in range(0,dim):
                sum[j] += f[j]
                sumSqr[j] += f[j]**2
            f=func(funcArg)
        out = []
        for j in range(0,dim):
            out.append(sum[j]/numAverages)
            out.append(numpy.sqrt(sumSqr[j]/numAverages - (sum[j]/numAverages)**2))
        return out