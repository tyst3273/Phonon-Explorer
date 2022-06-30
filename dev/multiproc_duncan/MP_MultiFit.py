#________________________________________________________________________%
#              University of Colorado Boulder                            %
#                                                                        %
#              Dmitry Reznik, Irada Ahmadova, Aaron Sokolik                             %       
#                                                                        %
#    Work supported by the DOE, Office of Basic Energy Sciences,         %
#    Office of Science, under Contract No. DE-SC0006939                  %                  
#________________________________________________________________________%

from TextFile import *
from Data import *
from FittingFunction import *
from FittingData import *
import re
import os
import math
from numpy import *
from FitParameters import *
from plotDataWithFit import *
import datetime
from RSE_Constants import *
from Display import *
import time
from multiprocessing import Process

starttime = time.time()

#my_file = open(input_path)
#string_list = my_file.readlines() #This reads the input parameters file and turns in into a readable list of strings
#my_file.close()

def MultiFit(direction,q):
#    string_list[25]= f"projectRootDir=C:/Users/Duncan/Documents/TiO2/327K/{direction}/{q}/\n"
#    string_list[188]= f"fileWithGuesses=C:/Users/Duncan/Documents/TiO2/327K/{direction}/{q}/Guesses.txt\n"
#    my_file = open(input_path, "w")
#    new_file_contents = "".join(string_list)
#    my_file.write(new_file_contents)
#    my_file.close()
    print("Parameters")
    root_run = 'C:/Users/Duncan/Documents/TiO2/327K'
    params=Parameters(f'{root_run}/{direction}/{q}/',RSE_Constants.INPUTS_FILENAME_MAIN)
    params.ReadMultizoneFitParams()
    params.locationForOutputParam=params.folderForBkgSubtractedFiles
    for i in range (0,len(params.positionGuessesList[:])):
        print("DataSmall_q")
        if len(params.reducedQlist[i:i+1])==0:
            data=DataSmall_q(params,params.folderForBkgSubtractedFiles)
#            params.qh=1000  #1000 means that reduced q is not specified
#            params.qk=1000
#            params.ql=1000
            paramFileName='_'+RSE_Constants.FITTING_PARAM_FILE
    
        else:
            data=DataSmall_q(params,params.folderForBkgSubtractedFiles,params.reducedQlist[i:i+1][0])
            params.qh=params.reducedQlist[i:i+1][0][0]
            params.qk=params.reducedQlist[i:i+1][0][1]
            params.ql=params.reducedQlist[i:i+1][0][2]
            paramFileName='_'+RSE_Constants.FITTING_PARAM_FILE+'_'+str(params.qh)+'_'+str(params.qk)+'_'+str(params.ql)+'.txt'
    
        data.Read()
        print("InitialGuesses")
        params.positionGuesses=params.positionGuessesList[i]
        params.NumberofPeaks=len(params.positionGuesses)
        InitialGuess=InitialGuesses(params,data)
    
        #print("Fitting")
        Fitting=FittingData(params,InitialGuess,data)
    
        popt, pcov=Fitting.doFitting()
    
        WriteToText=FitParameters(popt,data.filenames)
        WriteToText.writeToFile(params,params.locationForOutputParam+paramFileName+'.txt')
        
        WriteToText=FitParameters(numpy.sqrt(numpy.diag(pcov)),data.filenames)
        WriteToText.writeToFile(params,params.locationForOutputParam+'err'+paramFileName+'.txt')
    
#        WriteToText=WriteFinalParamToFile(params,data)
    
 #       WriteToText.writeToFile(params.locationForOutputParam,paramFileName,popt)
 #       WriteToText.writeToFile(params.locationForOutputParam,"err"+paramFileName,numpy.sqrt(numpy.diag(pcov)))
        folder=params.locationForOutputParam
        PlotDataWithFitting=PlotDataWithFitParamCustomFolder(params,folder,folder,folder)
        Disp=Display()
        #from FitParameters import *
        Disp.MakePlotSummary(params.folderForBkgSubtractedFiles,params.ProcessedDataName)

if __name__ == '__main__':
    direction = ["GM", "GZ", "GX", "MA", "XR", "ZR"]
    #direction1 = ["GM", "GZ", "GX"]
    #direction2 = ["MA", "XR", "ZR"]
    jobs = []
    
    for i in direction:
        for q in (0,0.1,0.2,0.3,0.4,0.5): #Iterates over all reduced q vectors in each direction
        #q=0
            print(f"Fitting {i} {q}")
            fitjob = Process(target=MultiFit, args=(i,q)) #Runs cuts() function using i and q as input arguments
            jobs.append(fitjob)
            #time.sleep(1)
            fitjob.start() #Initializes each job
            print(f"Fit of {i} {q} completed")
            #time.sleep(1) #I included these so it won't skip some iterables
            
    #Makes each job run    
    for fitjob in jobs:
        print("fitjobs")
        fitjob.join()
        time.sleep(1)
                
print("\nFinished")
print(f'Time taken = {time.time() - starttime} seconds')

