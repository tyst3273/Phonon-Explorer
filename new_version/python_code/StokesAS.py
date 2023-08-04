#________________________________________________________________________%
#              University of Colorado Boulder                            %
#                                                                        %
#              Dmitry Reznik, Irada Ahmadova                             %       
#                                                                        %
#    Work supported by the DOE, Office of Basic Energy Sciences,         %
#    Office of Science, under Contract No. DE-SC0006939                  %                  
#________________________________________________________________________%


import os
import numpy as np
from Display import *
from collections import namedtuple
import matplotlib.pyplot as plt
from SubtractPolyBackgr import *

#ProjectDir='StretchingTests100'
#ProjectDir='0.5,0,0/'

location_ForPlots='/home/ty/Desktop/phonon_explorer/test/test/stokes/'
#location_ForPlots='/users/dmitryreznik/To be backed up/working files/LSr0x2MnO3JPark/StokesAS/'+ProjectDir
if not os.path.isdir(location_ForPlots):
    os.makedirs(location_ForPlots)

#datasetdir1='E:/Data_LSNO_120meV_240K/'+ProjectDir + '/good_slices/'
#datasetdir2='E:/Data_LSNO_120meV_450K/'+ProjectDir + '/good_slices/'
#datasetdir1='/users/dmitryreznik/To be backed up/working files/LSr0x2MnO3JPark/nxspe054Ei_010K_Ebin0p5ZB/'+ProjectDir + 'good_slices/'
datasetdir1='/home/ty/Desktop/phonon_explorer/test/test/good_slices/'

filenames1=[file for file in os.listdir(datasetdir1) if file.startswith("H")and not file.endswith("pdf")]

#print (filenames1)
#print (filenames2)
T1=300

M1=1

font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }
for i in range(0,len(filenames1)):
    
    print(filenames1[i])
    elasticLineWidth=4 #FWHM in meV
    offset=3     #Make sure that the offset is such that it would ensure that both stokes and antistokes sides have at least3 points
    AllData=np.genfromtxt(datasetdir1+filenames1[i])
    intensity=AllData[:,][:,1]
    energy=AllData[:,][:,0]
    error=AllData[:,][:,2]
    
    
    indices0 = [index for index, value in enumerate(energy) if -elasticLineWidth <= value <= elasticLineWidth]
    
    if len(indices0)>0:
#Make sure that the offset is such that it would ensure that both stokes and antistokes sides have at least 3 points
        if min(energy)<-elasticLineWidth-offset and max(energy)>elasticLineWidth+offset:
            print(min(energy),max(energy))
            energy1=energy[max(indices0):]
            energy2=-energy[:min(indices0)]
            intensity1=intensity[max(indices0):]
            intensity2=intensity[:min(indices0)]
            error1=error[max(indices0):]
            error2=error[:min(indices0)]
            energyS,DataWithoutBackgroundS,error1S,backgroundS = SutractPolyBackground1File(energy1,intensity1,error1,0)
            energyA,DataWithoutBackgroundA,error1A,backgroundA = SutractPolyBackground1File(energy2,intensity2,error2,0)
        
            intensityS=np.zeros(len(energyS))
            errorS=np.zeros(len(energyS))
            intensityA=np.zeros(len(energyA))
            errorA=np.zeros(len(energyA))

        
            for ii in range(0,len(energyS)):
                intensityS[ii]=M1*DataWithoutBackgroundS[ii]/(1+1/(np.exp(energyS[ii]/(0.08617*T1))-1))
                errorS[ii]=M1*error1S[ii]/(1+1/(np.exp(energyS[ii]/(0.08617*T1))-1))
    
            for ii in range(0,len(energyA)):
#                intensityA[ii]=M1/(1+1/(np.exp(energyA[ii]/(0.08617*T1))-1))
                intensityA[ii]=M1*DataWithoutBackgroundA[ii]/(1/(np.exp(energyA[ii]/(0.08617*T1))-1))
                errorA[ii]=M1/error1A[ii]/(1+1/(np.exp(energyA[ii]/(0.08617*T1))-1))


#            plt.errorbar(energyS,DataWithoutBackgroundS,error1,fmt='o')
            plt.errorbar(energyS,intensityS,errorS,fmt='o')
            plt.errorbar(energyA,intensityA,error1A,fmt='*')
            plt.xlabel('Energy')
            plt.ylabel('Intensity')
            plt.grid()
            plt.title(filenames1[i])
            plt.ylim(-0.5,1)
            #plt.text(20, 0.001,filenames1[i],fontdict=font)

            plt.savefig(location_ForPlots+'/'+filenames1[i]+'a4'+'portrait'+'.pdf')
            print(location_ForPlots+'/'+filenames1[i]+'a4'+'portrait'+'.pdf')
#            plt.show()
            plt.close()
Disp=Display()
Disp.MakePlotSummary(location_ForPlots,'')

