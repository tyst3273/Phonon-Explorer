
from TextFile import *
from Data import *
from FittingFunction import *
from FittingData import *
from Display import *
import re
import os
import math
import time
from numpy import *
from plotDataWithFit import *
from Background import *
import matplotlib.pyplot as plt
from RSE_Constants import *


def Adjust(H,K,L,Intercept,slope, maxY, inFolder,T=0.1):
    outFolder=params.folderForBkgSubtractedFiles
    filename=str(RSE_Constants.FILENAME_FORMAT % (H,K,L))
    data=Dataset(inFolder,[filename])
    data.DivideByBoseFactorNorm(T)
    data.SubtractLine(Intercept,slope)
    dataFile=DataTextFile(params.folderForBkgSubtractedFiles,filename)
    dataFile.Write(data.Energy,data.Intensity,data.Error)
    plt.errorbar(data.Energy,data.Intensity,data.Error,fmt='o')
#    plt.plot(data.Energy, data.Intensity, '*', label='data')
    plt.xlabel(RSE_Constants.X_LABEL)
    plt.ylabel(RSE_Constants.Y_LABEL)
    plt.title(filename)
                #plt.axis([0, 40, 0, 0.0005])
    plt.ylim((-0.1*maxY,maxY)) 
    plt.grid() 
    plt.savefig(outFolder+filename+'.pdf')
                #plt.show()
    plt.close()
    return

#factor=input('Enter factor (1 or -1):')

# the filename is now gotten INSIDE parameters; can pass thru cmd line args as '-i input_file_path'
params = Parameters()

inFolder=params.folderForBkgSubtractedFiles
Qs=np.genfromtxt(inFolder+"BackgroundAdjustment.txt")
print(inFolder+"BackgroundAdjustment.txt")
QHlist=Qs[:,][:,0]
QKlist=Qs[:,][:,1]
QLlist=Qs[:,][:,2]
InterceptList=Qs[:,][:,3]
SlopeList=Qs[:,][:,4]
try:
    T=float(sys.argv[1])
except:
    T=1
    
for i in range(0,len(QHlist)):
    Adjust(QHlist[i],QKlist[i],QLlist[i],float(InterceptList[i]),float(SlopeList[i]),params.maxY,inFolder,T)

folder=params.folderForBkgSubtractedFiles
#PlotDataWithFitting=PlotDataWithFitParamCustomFolder(params,folder,folder,folder)
Disp=Display()
Disp.MakePlotSummary(folder,params.ProcessedDataName)


