#________________________________________________________________________%
#              University of Colorado Boulder                            %
#                                                                        %
#              Dmitry Reznik                             %
#                                                                        %
#    Work supported by the DOE, Office of Basic Energy Sciences,         %
#    Office of Science, under Contract No. DE-SC0006939                  %                  
#________________________________________________________________________%


import os
import numpy as np
from collections import namedtuple
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from TextFile import *
from Display import *
from plotDataWithFit import *
from RSE_Constants import *

# define the quadratic function to fit


def delete_points_at_indices(lst, indices):
        indices = set(indices)
        return [x for i, x in enumerate(lst) if i not in indices]


# Define the parabolic function
def parabola(x, a, b, c):
#    print(RSE_Constants.Temperature)
#    return (1+1/(np.exp(x/(0.08617*RSE_Constants.Temperature))-1))*((a*x)**2.5 + b * x + c)
    return ((a*x)**2.5 + b * x + c)

def trim(background,x,y,err):
    IndecesOfPointsToBeDeleted = []
    i=0
    while y[i] - background[i]>0:
        IndecesOfPointsToBeDeleted.append(i)
        i=i+1

    i=len(y)-1
    while y[i] - background[i]>0:
        IndecesOfPointsToBeDeleted.append(i)
        i=i-1

    xnew=delete_points_at_indices(x,IndecesOfPointsToBeDeleted)
    ynew=delete_points_at_indices(y,IndecesOfPointsToBeDeleted)
    errnew=delete_points_at_indices(err,IndecesOfPointsToBeDeleted)
    return np.array(xnew),np.array(ynew),np.array(errnew)

def DeletePointsAboveLine(a,b,c,x,y,err):
    curve = parabola(x, a, b, c)
    IndecesOfPointsToBeDeleted = []
    for i in range(0,len(x)):
        if y[i] - 1.5*err[i]-curve[i]>0:
            IndecesOfPointsToBeDeleted.append(i)
 #   print(IndecesOfPointsToBeDeleted)
    xnew=delete_points_at_indices(x,IndecesOfPointsToBeDeleted)
    ynew=delete_points_at_indices(y,IndecesOfPointsToBeDeleted)
    errnew=delete_points_at_indices(err,IndecesOfPointsToBeDeleted)
    return np.array(xnew),np.array(ynew),np.array(errnew)

# fit the parabola to the data using curve_fit from SciPy
def generateFit(x,y,errors):
    
#    popt, pcov = curve_fit(weighted_parabola, x, y, sigma=error)

    popt, pcov = curve_fit(parabola, x, y, sigma=errors)
    a, b, c = popt
    a_err, b_err, c_err = np.sqrt(np.diag(pcov))

    return a,b,c

def SutractPolyBackground1File(energy,intensity,error,Trim):
    NPoints=len(energy)+1
    
    energyTemp=energy
    intensityTemp=intensity
    errorTemp=error

    ii = 0
    while len(energy)-NPoints<0:
        if ii>10:
            break
        NPoints=len(energyTemp)
        a,b,c=generateFit(energyTemp,intensityTemp,errorTemp)
        energyTemp,intensityTemp,errorTemp=DeletePointsAboveLine(a,b,c,energyTemp,intensityTemp,errorTemp)
        ii += 1
    fit=parabola(energy, a, b, c)
    if Trim==1:
        energy,intensity,error=trim(fit,energy,intensity,error)
        fitNew=parabola(energy, a, b, c)
    else:
        fitNew=fit
    DataWithoutBackground=intensity-fitNew

    return energy,DataWithoutBackground,error,fitNew

def SubtractPolyBackground():

    params = Parameters()
    params.ReadPolyBackgroundParams()
#    RSE_Constants.Temperature=params.Temperature

    folderForBkgSubtractedFiles=[subdir for subdir in os.listdir(params.path_data) if subdir.startswith(RSE_Constants.BACKGROUND_SUBTRACTED_FOLDER)]

    folder=params.path_data

    files=[file for file in os.listdir(folder) if file.startswith(RSE_Constants.STARTS_WITH) and not file.endswith(RSE_Constants.ENDS_WITH)]

    print(RSE_Constants.BACKGROUND_SUBTRACTED_FOLDER, os.path.isdir(params.folderForBkgSubtractedFiles))
    if not os.path.isdir(params.folderForBkgSubtractedFiles):
        os.makedirs(params.folderForBkgSubtractedFiles)


    for i in range (0,len(files)):

        try:
            font = {'family': 'serif',
            'color':  'darkred',
            'weight': 'normal',
            'size': 16,
            }

        #print(filenames1[i],filenames2[j])
            n=params.PolyFitNumberIgnoredPointsAtEnd #allows deleting n points at the end of the file
            filename=folder+files[i]
            AllData=np.genfromtxt(filename)
            AllData1=np.array(AllData[:-n])
    #       print(AllData)
            intensity=AllData1[:,1]
            energy=AllData1[:,0]
            error=AllData1[:,2]
            energy,DataWithoutBackground,error,background = SutractPolyBackground1File(energy,intensity,error,params.Trim)
            print(files[i])
#    print(params.folderForBkgSubtractedFiles)
#    if len(AllData1[:,0])-len(energy)>0.6*len(AllData1[:,0]):
            TxtFile=open(params.folderForBkgSubtractedFiles+files[i],'w+')
            for j in range (0,len(energy)):
                TxtFile.write(str(energy[j])+'  '+str(DataWithoutBackground[j])+'  '+str(error[j])+'\n')
            TxtFile.close()
    
            BackgroundCurve=[]
#        BackgroundCurve.append(AllData1[:,0])
            BackgroundCurve.append(energy)
            BackgroundCurve.append(background)
            plt=Plot(params)
            plt.plotSingle(params.path_data,params.path_data,[files[i]],BackgroundCurve)
            plt.plotSingle(params.folderForBkgSubtractedFiles,params.folderForBkgSubtractedFiles,[files[i]])
        except:
            a=1
    Disp=Display()
    Disp.MakePlotSummary(params.folderForBkgSubtractedFiles,params.ProcessedDataName)
    Disp.MakePlotSummary(params.path_data,params.ProcessedDataName)

#SubtractPolyBackground()
#print("here")
