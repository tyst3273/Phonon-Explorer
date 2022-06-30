#________________________________________________________________________%
#              University of Colorado Boulder                            %
#                                                                        %
#              Dmitry Reznik, Irada Ahmadova, Aaron Sokolik              %       
#                                                                        %
#    Work supported by the DOE, Office of Basic Energy Sciences,         %
#    Office of Science, under Contract No. DE-SC0006939                  %                  
#________________________________________________________________________%


#This code makes cuts in an sqw file at specified reduced q vectors along specific reciprocal lattic vectors.

from TextFile import *  #Provides definitions for input parameters
from Data import *
#import matlab.engine
from FittingFunction import *
from FittingData import *
from RSE_Constants import * #Tell where input parameters are located, among other things
from plotDataWithFit import *   #This tells the program how to plot everything
from Display import *   #Creates merged PDF at the end of the program
import re
import os
import math
import time
from numpy import *
from FitParameters import *
from Background import *
import matplotlib.pyplot as plt

import gc   #garbage collection
#import multiprocessing as mp
from shutil import copy
from multiprocessing import Process, Pool, set_start_method, Manager   #Currently I only use Process and not the others

starttime = time.time() #Just to keep track of how long the program takes

root_path = 'C:/Users/Duncan/Documents'
dir_path = f'{root_path}/TiO2' #This is where the program and output folders are
input_path = f'{root_path}/Phonon-Explorer-master/Input_Files/InputParameters.txt'
run = '327K'

my_file = open(input_path)
string_list = my_file.readlines() #This reads the input parameters file and turns in into a readable list of strings
my_file.close()

gc.enable() #Turns on garbage collection
Disp=Display()
def CopyInput(directions):
    for direction in directions: 
        for q in (0,0.1,0.2,0.3,0.4,0.5):
            print(f'Direction: {direction}')
            output = f'{run}/{direction}/{q}/'
            src_1 = f'{dir_path}/Input/guru99.py'
            src_2 = input_path
            dest = f'{dir_path}/{output}'
            string_list[25]= f'projectRootDir={dir_path}/{output}\n' #Changes output directory
            string_list[43]= f'InputFilesDir={dir_path}/{output}\n'
            if direction == "GM":
                string_list[92]= f"qh={q}\n"
                string_list[93]= f"qk={q}\n"
                string_list[94]= "ql=0\n"
            if direction == "GZ":
                string_list[92]= "qh=0\n"
                string_list[93]= "qk=0\n"
                string_list[94]= f"ql={q}\n"
            if direction == "GX":
                string_list[92]= f"qh={q}\n"
                string_list[93]= "qk=0\n"
                string_list[94]= "ql=0\n"
            if direction == "MA":
                string_list[92]= "qh=0.5\n"
                string_list[93]= "qk=0.5\n"
                string_list[94]= f"ql={q}\n"
            if direction == "XR":
                string_list[92]= "qh=0.5\n"
                string_list[93]= "qk=0\n"
                string_list[94]= f"ql={q}\n"
            if direction == "ZR":
                string_list[92]= f"qh={q}\n"
                string_list[93]= "qk=0\n"
                string_list[94]= "ql=0.5\n"
            if direction == "GA":
                string_list[92]= f"qh={q}\n"
                string_list[93]= f"qk={q}\n"
                string_list[94]= f"ql={q}\n"
            my_file = open(input_path, "w")
            new_file_contents = "".join(string_list)
            my_file.write(new_file_contents)
            my_file.close()
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            copy(src_1,dest) #Copies guru99.py and InputParameters.txt to each q subfolder. Prevents errors from multiple processes reading the same file at the same time.
            copy(src_2,dest)

def cuts(direction,q):
    params=Parameters(f'{dir_path}/{run}/{direction}/{q}/', RSE_Constants.INPUTS_FILENAME)
    #time.sleep(1)
    print(f'{direction} {q}')
    #Main function that makes cuts
    if params.QMode==0:
        dd=1
        print(f"Current save path: {params.path_data}")
        testData=DataSmall_q(params, params.path_data)
    if params.QMode==1:
        testData=CollectionOfQs(params)
    testData.Generate()
    if params.BkgMode==0:
        plot=Plot(params)
        plot.Plot(direction,q)
        Disp=Display()
        Disp.MakePlotSummary(params.path_data,params.ProcessedDataName) 
    print("Plots Completed")
    #print("line 21, explore") #I don't know what this or the later ones mean.
    if params.BkgMode==1:
        randomFiles=DataBackgroundQs(params)
        print("line 23, explore")
        randomFiles.GenerateAllFiles()
        print("line 25, explore")   

def subBack(direction,q):
    params=Parameters(f'{dir_path}/{run}/{direction}/{q}/', RSE_Constants.INPUTS_FILENAME)
    params.ReadBackgroundParams()
    folderForBkgSubtractedFiles=[subdir for subdir in os.listdir(params.path_data) if subdir.startswith(RSE_Constants.BACKGROUND_SUBTRACTED_FOLDER)]
    for i in range (0,len(folderForBkgSubtractedFiles)):
        folder=params.path_data+folderForBkgSubtractedFiles[i]+'/'
        print(folder)
        files=[file for file in os.listdir(folder) if file.startswith(RSE_Constants.STARTS_WITH) and not file.endswith(RSE_Constants.ENDS_WITH)] 
        rawFileName=folder[folder.index('_' + RSE_Constants.STARTS_WITH)+1:-1]
    #    print(rawFileName)
        rawData=Dataset(params.path_data,[rawFileName])
    #    print(rawData.Intensity)
        for j in range (0,len(files)):
            try:
    #        if 1==1: 
                Fitting=FittingBackgroundData(params,folder,files[j])
                popt, pcov=Fitting.doFitting()
                data=Dataset(folder,[files[j]])
    #            WriteToText=WriteFinalParamToFile(params,data)
    #            WriteToText.writeToFile(folder,files[j],popt)
                WriteToText=FitParameters(popt,data.filenames)
                WriteToText.writeToFile(params,folder+'_'+files[j]+'.txt')
            except Exception as e:
                print("fit failed:"+folder+" "+files[j])
                print(e)
    #            time.sleep(5)
        PlotDataWithFitting=PlotDataWithFitParamCustomFolder(params,folder,folder,folder)
        Backgr=Background(params,folder,rawData)
        Backgr.DisplayAllFiles(folder,folder,rawFileName)
        Backgr.DisplayOrigFile(params.path_data,params.location_ForPlots,rawFileName)
        Disp.MakePlotSummary(folder,RSE_Constants.BACKGR_PREFIX+rawFileName) #Make single PDF of all plots in the background folder
    Disp.MakePlotSummary(params.path_data,params.ProcessedDataName)
    #Backgr.Adjust(rawFileName,params.folderForBkgSubtractedFiles)
    Disp.MakePlotSummary(params.folderForBkgSubtractedFiles,params.ProcessedDataName)

#These are my attempts at using the Pool class. For some reason it doesn't work. If it works at all, it only does one process at a time. 
#It's in shambles right now.

#def main():
#    man = Manager()
#    direction = ["GM", "GX", "ZR"]
#    q_dir = man.list([])
#    for d in direction:
#        for q in (0,0.1,0.2,0.3,0.4,0.5):
#            r = (d,q)
#            q_dir.append(tuple(r))
#    print(q_dir)
#    with Pool(12) as pool:
#        pool.starmap(cuts, q_dir)
#        
#        #time.sleep(1)
#        #params=Parameters(RSE_Constants.INPUTS_PATH, RSE_Constants.INPUTS_FILENAME)
#    print("\nFinished")
#    print('Time taken = {} seconds'.format(time.time() - starttime))
    
#if __name__ == '__main__':
#    set_start_method("spawn")
#    main() #Either uncomment this line and the above function or just the code below.
    #direction = ["GM", "GZ", "GX", "MA", "XR", "ZR"]
    #direction = ["GM", "ZR"]
    #p = Pool(processes=4) #Here, you can explicitly state how many cores/logical processors to use. The Process class does not seem to have this ability.
    #for i in direction:
    #    for q in (0,0.1,0.2,0.3,0.4,0.5):    
    #        p.map(cuts, i,q) #I use starmap above since map cannot iterate over more than 1 iterable.
    #print("\nFinished")
    #print('Time taken = {} seconds'.format(time.time() - starttime))

#Main multiprocessing code. Currently uses the Process class instead of the Pool class.
#I attemped to add in a Queue, but no luck yet.
#Works fine with 3-4 directions. If I change to all directions, it is too much for my PC to handle.
if __name__ == '__main__':
    #set_start_method("spawn")
    AllDirections = ["GM", "GZ", "GX", "MA", "XR", "ZR", "GA"]
    direction1 = ["GM", "GZ", "GX"]
    direction2 = ["MA", "XR", "ZR"]
    direction3 = ["GA"]
    direction_sets = [direction1,direction2,direction3]
    jobs = []
    CopyInput(AllDirections)
    
    for subset in direction_sets:
        for i in subset:
            for q in (0,0.1,0.2,0.3,0.4,0.5): #Iterates over all reduced q vectors in each direction
            #q=0
                print("Make_cut")
                cutjob = Process(target=cuts, args=(i,q)) #Runs cuts() function using i and q as input arguments
                jobs.append(cutjob)
                #time.sleep(1)
                cutjob.start() #Initializes each job
                #time.sleep(1) #I inclued these so it won't skip some iterables
                
        for cutjob in jobs: 
            print("cutjobs")
            cutjob.join()   #Prevents later jobs from starting before these are finished
            time.sleep(1)
        
        for i in subset:
            for q in (0,0.1,0.2,0.3,0.4,0.5): #Iterates over all reduced q vectors in each direction
                print("Make_sub")
                subjob = Process(target=subBack, args=(i,q)) #Runs cuts() function using i and q as input arguments
                jobs.append(subjob)
                #time.sleep(1)
                subjob.start() 
                time.sleep(1) 
        for subjob in jobs:
            print("subjobs")
            subjob.join()
            time.sleep(1)
    
#    for i in direction2:
#        for q in (0,0.1,0.2,0.3,0.4,0.5): #Iterates over all reduced q vectors in each direction
#        #q=0
#            print("Make_cut")
#            cutjob = Process(target=cuts, args=(i,q)) #Runs cuts() function using i and q as input arguments
#            jobs.append(cutjob)
#            #time.sleep(1)
#            cutjob.start() 
#            time.sleep(1) 
#            
#    #Makes each job run    
#    for cutjob in jobs:
#        print("cutjobs")
#        cutjob.join()
#        time.sleep(1)
#    
#    for i in direction2:
#        for q in (0,0.1,0.2,0.3,0.4,0.5): #Iterates over all reduced q vectors in each direction
#        #q=0
#            print("Make_sub")
#            subjob = Process(target=subBack, args=(i,q)) #Runs cuts() function using i and q as input arguments
#            jobs.append(subjob)
#            #time.sleep(1)
#            subjob.start() #Initializes each job
#            time.sleep(1) #I inclued these so it won't skip some iterables
#    #Makes each job run    
#    for subjob in jobs:
#        print("subjobs")
#        subjob.join()
#        time.sleep(1)    
        
    print("\nFinished")
    print(f'Time taken = {time.time() - starttime} seconds')
