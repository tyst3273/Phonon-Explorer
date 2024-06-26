#________________________________________________________________________%
#              University of Colorado Boulder                            %
#                                                                        %
#              Dmitry Reznik, Irada Ahmadova, Aaron Sokolik                             %       
#                                                                        %
#    Work supported by the DOE, Office of Basic Energy Sciences,         %
#    Office of Science, under Contract No. DE-SC0006939                  %                  
#________________________________________________________________________%



import numpy as np
import math
import os
from RSE_Constants import RSE_Constants
from Utils import arg_parser


# --------------------------------------------------------------------------------------------------

class TextFile():

    def __init__(self,FolderName,FileName):
        self.filename = FileName
        self.foldername = FolderName

    # ----------------------------------------------------------------------------------------------
            
    def getFileName(self):
        return self.filename

    # ----------------------------------------------------------------------------------------------

    def getFolderName(self):
        return self.foldername

    # ----------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------

class Parameters(TextFile):

    def __init__(self):
        """
        a new interface to parameters that optionally gets the filename from cmd line args. if none
        are give, gets them from 'RSE_constants' ...
        """

        _args = arg_parser()
        _args.get_input_file()
        if _args.input_file is None:
            self._init_parameters(RSE_Constants.INPUTS_PATH,RSE_Constants.INPUTS_FILENAME)
        else:
            self._init_parameters(_args.input_path,_args.input_file)

    # ----------------------------------------------------------------------------------------------

    def _init_parameters(self,FolderName,FileName):

        TextFile.__init__(self,FolderName,FileName)

        with open(os.path.join(self.foldername,self.filename)) as f:
             parameters = f.read().splitlines()
        f.close()
        
        self.keyword=''
        #self.num_processes = self.evalIntWarning(self.ParseByKeyword('num_processes',parameters),1)
        #self.horace_threads = self.evalIntWarning(self.ParseByKeyword('horace_threads',parameters),4) # 4 is probably safe 
        self.BkgMode=self.evalIntWarning(self.ParseByKeyword("BkgMode",parameters),0)
        self.QMode=self.evalIntWarning(self.ParseByKeyword("QMode",parameters),0)
        self.sqw_path=self.evalError(self.ParseByKeyword("sqw_path",parameters))
        self.dataFileType=self.sqw_path[self.sqw_path.rfind('.')+1:].lower()
        self.projectRootDir=self.evalError(self.ParseByKeyword("projectRootDir",parameters))

        if self.dataFileType=="sqw":
            self.rawDataClassFile="SQWAccess"
            from sys import platform
            if platform == "linux" or platform == "linux2":
                print(platform)
                self.HoracePath=self.evalError(self.ParseByKeyword("HoracePath",parameters))
        elif self.dataFileType=="nxs":
            self.rawDataClassFile="NXSAccess"
        elif self.dataFileType == 'hdf5':
            self.rawDataClassFile='hdf5Access'
        else:
            msg = "Raw data file extension is not valid. Must be either .sqw, .nxs, or .hdf5"
            print(msg)
            raise Exception(msg)
            return

        self.ProcessedDataName=self.evalError(self.ParseByKeyword("ProcessedDataName",parameters))
        self.path_data=self.evalError(self.ParseByKeyword("projectRootDir",parameters))+self.ProcessedDataName+'/good_slices/'
        print(os.getcwd())
        print(self.path_data)
        if not os.path.isdir(self.path_data):
            os.makedirs(self.path_data)

        self.Projection_u=np.asarray(self.evalError(self.ParseByKeyword("Projection_u",parameters)).split(",")).astype(float)
        self.Projection_v=np.asarray(self.evalError(self.ParseByKeyword("Projection_v",parameters)).split(",")).astype(float)
        self.path_InputFiles=self.evalError(self.ParseByKeyword("InputFilesDir",parameters))
        self.ErrorToIntensityMaxRatio=self.evalRealWarning(self.ParseByKeyword("ErrorToIntensityMaxRatio",parameters),1)
        if self.QMode==1:    
            self.textfile_for_selectedQs=self.evalError(self.ParseByKeyword("textfile_for_selectedQs",parameters))
        if self.QMode==0:
            self.qh=self.evalRealWarning(self.ParseByKeyword("qh",parameters),0)
            self.qk=self.evalRealWarning(self.ParseByKeyword("qk",parameters),0)
            self.ql=self.evalRealWarning(self.ParseByKeyword("ql",parameters),0)
            self.h_start=self.evalIntWarning(self.ParseByKeyword("h_start",parameters),0)
            self.h_end=self.evalIntWarning(self.ParseByKeyword("h_end",parameters),0)
            self.k_start=self.evalIntWarning(self.ParseByKeyword("k_start",parameters),0)
            self.k_end=self.evalIntWarning(self.ParseByKeyword("k_end",parameters),0)
            self.l_start=self.evalIntWarning(self.ParseByKeyword("l_start",parameters),0)
            self.l_end=self.evalIntWarning(self.ParseByKeyword("l_end",parameters),0)
        self.e_start=eval(self.evalError(self.ParseByKeyword("e_start",parameters)))
        self.e_end=eval(self.evalError(self.ParseByKeyword("e_end",parameters)))
        self.e_step=eval(self.evalError(self.ParseByKeyword("e_step",parameters)))
        self.Deltah=eval(self.evalError(self.ParseByKeyword("Deltah",parameters)))
        self.Deltak=eval(self.evalError(self.ParseByKeyword("Deltak",parameters)))
        self.Deltal=eval(self.evalError(self.ParseByKeyword("Deltal",parameters)))

        self.Offset_H=self.evalReaNoWarning(self.ParseByKeyword("Offset_H",parameters),0)
        self.Offset_K=self.evalReaNoWarning(self.ParseByKeyword("Offset_K",parameters),0)
        self.Offset_L=self.evalReaNoWarning(self.ParseByKeyword("Offset_L",parameters),0)
        self.Offset_E=self.evalReaNoWarning(self.ParseByKeyword("Offset_E",parameters),0)

        if not self.Offset_H==0 or not self.Offset_K==0 or not self.Offset_L==0 or not self.Offset_E==0:
            print ("WARNING: Nonzero offsets for H,K,L,or E. Are you sure?????????")
            print ("Offset_H:  ",self.Offset_H)
            print ("Offset_K:  ",self.Offset_K)
            print ("Offset_L:  ",self.Offset_L)
            print ("Offset_E:  ",self.Offset_E)
            
        self.MinPointsInDataFile=self.evalIntWarning(self.ParseByKeyword("MinPointsInDataFile",parameters),10)
        self.location_ForPlots=self.path_data
        self.maxY=eval(self.evalError(self.ParseByKeyword("maxY",parameters)))
        self.dataFileNameStart=self.evalWarning(self.ParseByKeyword("dataFileNameStart",parameters),"H")
        self.folderForBkgSubtractedFiles=self.projectRootDir+self.ProcessedDataName+'/subtr_background/'
        self.SmallqAlgorithm=self.evalError(self.ParseByKeyword("SmallqAlgorithm",parameters))
        
    # ----------------------------------------------------------------------------------------------

    def ReadBackgroundParams(self):
        with open(os.path.join(self.foldername,self.filename)) as f:
             parameters = f.read().splitlines()

        self.a=eval(self.evalError(self.ParseByKeyword("a",parameters)))
        self.b=eval(self.evalError(self.ParseByKeyword("b",parameters)))
        self.c=eval(self.evalError(self.ParseByKeyword("c",parameters)))
        self.phiRange=0#eval(self.ParseByKeyword("phiRange",parameters)) Not used now
        self.thetaRange=0#eval(self.ParseByKeyword("thetaRange",parameters)) Not used now
        self.NumberOfTries=self.evalIntWarning(self.ParseByKeyword("NumberOfTries",parameters),10)
        self.maxFiles=self.evalIntWarning(self.ParseByKeyword("maxFiles",parameters),10)
#        self.Resolution=eval(self.evalError(self.ParseByKeyword("Resolution",parameters)))
        self.MinPointsInDataBackgroundFile=self.evalIntWarning(self.ParseByKeyword("MinPointsInDataBackgroundFile",parameters),10)
        self.Resolution=eval(self.evalError(self.ParseByKeyword("Resolution",parameters)))
        self.BackgroundAlgorithm=self.evalWarning(self.ParseByKeyword("BackgroundAlgorithm",parameters),"Standard")
#        if self.BackgroundAlgorithm!="Standard":
#            self.maxFiles=2
        self.NumberofPeaks=0
        self.ReadSharedParams(parameters)

    # ----------------------------------------------------------------------------------------------
    def ReadPolyBackgroundParams(self):
        with open(os.path.join(self.foldername,self.filename)) as f:
             parameters = f.read().splitlines()

        self.PolyFitNumberIgnoredPointsAtEnd=self.evalIntWarning(self.ParseByKeyword("PolyFitNumberIgnoredPointsAtEnd",parameters),0)
        
        if self.PolyFitNumberIgnoredPointsAtEnd==1:
             print("Last "+str(self.PolyFitNumberIgnoredPointsAtEnd)+" points will be ignored in background determination to remove effect of detector edge")

             

        # ----------------------------------------------------------------------------------------------
        
        self.Trim=self.evalIntWarning(self.ParseByKeyword("Trim",parameters),0)
        
        if self.Trim==1:
             print("Last "+str(self.PolyFitNumberIgnoredPointsAtEnd)+" points will be ignored in background determination to remove effect of detector edge")
        else:
             self.Trim=0
        
        self.Temperature=self.evalReaNoWarning(self.ParseByKeyword("Temperature",parameters),3.0)
        
    def ReadMultizoneFitParams(self):
        with open(os.path.join(self.foldername,self.filename)) as f:
             parameters = f.read().splitlines()
        f.close()
        self.reducedQlist=[]
        self.positionGuessesList=[]
        self.WidthLowerBound=eval(self.ParseByKeyword("WidthLowerBound",parameters))
        self.fileWithGuesses=self.evalError(self.ParseByKeyword("fileWithGuesses",parameters))
        self.InitWidthsFinal=self.WidthLowerBound
#        try:
#            f=open(self.path_InputFiles+self.fileWithGuesses)
#            f.close()
#        except:
#            print ("WARNING: Position Guesses not specified")
        
        try:
#        if 1==1:
            with open(self.path_InputFiles+self.fileWithGuesses) as f:
                self.positionGuessesList.append([float(x) for x in next(f).split()])
                while True:
                    self.SmallqAlgorithm=self.evalError(self.ParseByKeyword("SmallqAlgorithm",parameters))
                    self.reducedQlist.append([float(x) for x in next(f).split()])
                    self.positionGuessesList.append([float(x) for x in next(f).split()])
        except:
            f.close()
            d=1   #dummy statement
#
        print(self.positionGuessesList)
        self.ReadSharedParams(parameters)

    # ----------------------------------------------------------------------------------------------

    def ReadSharedParams(self,parameters): #multizone and background
        self.locationForOutputParam=self.path_data
        self.InitAmplitudes=eval(self.ParseByKeyword("InitialAmplitude",parameters))
        self.NumberofFitIter=self.evalIntWarning(self.ParseByKeyword("NumberofIter",parameters),3)
    
    # ----------------------------------------------------------------------------------------------

    def getIndex(self,keyword,parameters):
        for i in range(0,len(parameters)):
            try:
                if keyword==parameters[i][:parameters[i].index('=')]:
                    return i
            except:
                d=1
        return -1

    # ----------------------------------------------------------------------------------------------

    def ParseByKeyword(self,keyword,parameters):
        Index=self.getIndex(keyword,parameters)
        self.keyword=keyword
        if Index==-1:
            return "None"
        return self.Parse(parameters[Index])
    def Parse(self,line):
        return line[line.index('=')+1:]
    def evalRealWarning(self,string,default):
        if string=="None":
            print("WARNING: "+self.keyword+" not specified")
            return default
        else:
            return eval(string)

    # ----------------------------------------------------------------------------------------------

    def evalReaNoWarning(self,string,default):
        if string=="None":
            return default
        else:
            return eval(string)
 
    # ----------------------------------------------------------------------------------------------

    def evalIntWarning(self,string,default):
        if string=="None":
            print("WARNING: "+self.keyword+" not specified")
            return default
        else:
            return int(string)

    # ----------------------------------------------------------------------------------------------

    def evalWarning(self,string,default):
        if string=="None":
            print("WARNING: "+self.keyword+" not specified")
            return default
        else:
            return string

    # ----------------------------------------------------------------------------------------------

    def evalError(self,string):
        if string=="None":
            strErr="Parameter ERROR: "+ self.keyword+" not specified"
            print(strErr)
            raise Exception(strErr)
        else:
            return string
    
# ----------------------------------------------------------------------------------------------

class DataTextFile(TextFile):
    
    # ----------------------------------------------------------------------------------------------

    def __init__(self,FolderName,FileName):
        TextFile.__init__(self,FolderName,FileName)
    
    # ----------------------------------------------------------------------------------------------

    def Read(self):
#        print (self.foldername+self.filename)
        data=np.genfromtxt(self.foldername+self.filename)
        return data
    
    # ----------------------------------------------------------------------------------------------

    def Write(self, Energy, Intensity, Error):

#        print os.path.isdir("E:/")
        if not os.path.isdir(self.foldername):
            os.makedirs(self.foldername)

        TxtFile=open(self.foldername+self.filename,'w+')    
        for i in range (0,len(Energy)):
            TxtFile.write(str(Energy[i])+'  '+str(Intensity[i])+'  '+str(Error[i])+'\n')
        TxtFile.close()

    # ----------------------------------------------------------------------------------------------



# --------------------------------------------------------------------------------------------------

class FitTextFile(TextFile):

    def __init__(self,FolderName,FileName):
        TextFile.__init__(self,FolderName,FileName)

    # ----------------------------------------------------------------------------------------------
        
    def Read(self):
#        print ('HERE  '+self.foldername+self.filename)
        data=np.genfromtxt(self.foldername+self.filename)
        return data

    # ----------------------------------------------------------------------------------------------
    
    def Write(self, FitResultsArray):

#        print os.path.isdir("E:/")
        if not os.path.isdir(self.foldername):
            os.makedirs(self.foldername)
        np.savetxt(self.foldername+self.filename,FitResultsArray)
        '''TxtFile=open(self.foldername+self.filename,'w+')
        writeString=""
        for i in range (0,len(FitResultsArray[0])):
            for j in range (0,len(FitResultsArray)):
                writeString=writeString+str(FitResultsArray[j][i]        
            TxtFile.write(writeString+'\n')
        TxtFile.close()'''

    # ----------------------------------------------------------------------------------------------

