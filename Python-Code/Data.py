#________________________________________________________________________%
#              University of Colorado Boulder                            %
#                                                                        %
#              Dmitry Reznik, Irada Ahmadova, Aaron Sokolik                             %       
#                                                                        %
#    Work supported by the DOE, Office of Basic Energy Sciences,         %
#    Office of Science, under Contract No. DE-SC0006939                  %                  
#________________________________________________________________________%


# --------------------------------------------------------------------------------------------------

from TextFile import *
import numpy
import itertools
from numpy import *
import math
import os
import random
import io
from RSE_Constants import *
import time
import timeit
import sys
from decimal import *
getcontext().prec=6
#from nexusformat.nexus import *


class Dataset: 
    
    """
    Dataset can be either a single cut at one Q or several cuts put together for the purposes of multizone fitting,
    In the former case it is initialized with an array of filenames that has one filename, in the latter case with 
    an array containing many filenames. 
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self,folder,filenames):

#        print(filenames)

        self.NumberofDatasets = len(filenames)
        self.Energy=[]
        self.Intensity=[]
        self.Error=[]
        self.Nint=[]
        self.Neng=[]
        self.filenames=filenames
        self.dataDirectory=folder
        self.appendDataset()
        self.tempIndex=1
        #self.eng = matlab.engine.start_matlab()

    # ----------------------------------------------------------------------------------------------

    def initialize(self):
        
        sys.path.insert(0,"Tools to access raw data")
        R=__import__(self.params.rawDataClassFile)
        RSE_Constants.rawData=R.RawData(self.params)
        RSE_Constants.FLAG=1

    # ----------------------------------------------------------------------------------------------

    def appendDataset(self):

        for m in range (0,self.NumberofDatasets):

#            print(self.dataDirectory,self.filenames[m])
            fileData=DataTextFile(self.dataDirectory,self.filenames[m])
            AllData=fileData.Read()
            energy=AllData[:,][:,0]+100*m
            intensity=AllData[:,][:,1] 
            error=AllData[:,][:,2]
            self.Energy=numpy.append(self.Energy,energy) 
            self.Intensity=numpy.append(self.Intensity,intensity)
            self.Error=numpy.append(self.Error,error)
    
    # ----------------------------------------------------------------------------------------------

    def clean(self):

        """
        removes data points where Intensity is zero or where error/intensity>self.params.ErrorToIntensityMaxRatio
        """

        index=[]
        for i in range (0,len(self.Energy)):
             if numpy.isnan(self.Intensity[i]):
                 self.Intensity[i]=0
        for i in range (0,len(self.Energy)):
            if self.Intensity[i]>0:
                if self.Error[i]/self.Intensity[i]>=self.params.ErrorToIntensityMaxRatio:
                    d=0
#                    self.Intensity[i]=0
        index.extend(np.nonzero(self.Intensity))
        Nint=zeros(len(index[0]))
        Neng=zeros(len(index[0]))
        Nerr=zeros(len(index[0]))
        for i in range (0,len(index[0])):
            Nint[i]=self.Intensity[index[0][i]]
            Neng[i]=self.Energy[index[0][i]]
            Nerr[i]=self.Error[index[0][i]]
        self.Intensity=Nint
        self.Energy=Neng
        self.Error=Nerr
        
    # ----------------------------------------------------------------------------------------------

    def smooth(self):

        firstEnergyIndex=0
        lastEnergyIndex=len(self.Energy)-1
#        print(params.Resolution, params.e_step,0.8*params.Resolution/params.e_step)
#        print(data.Energy[firstEnergyIndex:lastEnergyIndex:int(0.8*params.Resolution/params.e_step)])
        positions=self.Energy[firstEnergyIndex:lastEnergyIndex:int(params.Resolution/params.e_step)]
#        print(positions)
#        print([data.Energy[lastEnergyIndex]])
        if positions[len(positions)-1]<self.Energy[lastEnergyIndex]:
            positions=numpy.append(positions,[self.Energy[lastEnergyIndex]])
#        print(positions)
        params.positionGuesses=positions
        params.NumberofPeaks=len(positions)        
        params.InitWidthsFinal=params.MinPeakWidthForSmoothing
        params.WidthLowerBound=params.MinPeakWidthForSmoothing
        InitGuess=InitialGuesses(params,self)
        InitGuess.NumberofPeaks=len(positions) #overwrite
        fData=FittingData(params,InitGuess,self)
        self.fitParamsFofSmoothing=fData.doFitting()
    
    # ----------------------------------------------------------------------------------------------
        
    def removeNAN(self):
        Neng=np.zeros(1)
        Nint=np.zeros(1)
        Nerr=np.zeros(1)
        for i in range (0,len(self.Energy)):
             if not numpy.isnan(self.Intensity[i]):
                 Neng=numpy.append(Neng,self.Energy[i])
                 Nint=numpy.append(Nint,self.Intensity[i])
                 Nerr=numpy.append(Nerr,self.Error[i])
        self.Intensity=Nint[1:len(Nint)]
        self.Energy=Neng[1:len(Nint)]
        self.Error=Nerr[1:len(Nint)]

    # ----------------------------------------------------------------------------------------------

    def makeRawSlice(self, bin_h, bin_k, bin_l, bin_e,folder,file,minPoints,mult=1):

        bin_h[0]=bin_h[0]+self.params.Offset_H
        bin_h[1]=bin_h[1]+self.params.Offset_H
        bin_k[0]=bin_k[0]+self.params.Offset_K
        bin_k[1]=bin_k[1]+self.params.Offset_K
        bin_l[0]=bin_l[0]+self.params.Offset_L
        bin_l[1]=bin_l[1]+self.params.Offset_L
        bin_e[0]=bin_e[0]+self.params.Offset_E
        bin_e[2]=bin_e[2]+self.params.Offset_E

        RSE_Constants.rawData.GetSlice(bin_h, bin_k, bin_l, bin_e, self.params.Projection_u, self.params.Projection_v)

        bin_h[0]=bin_h[0]-self.params.Offset_H
        bin_h[1]=bin_h[1]-self.params.Offset_H
        bin_k[0]=bin_k[0]-self.params.Offset_K
        bin_k[1]=bin_k[1]-self.params.Offset_K
        bin_l[0]=bin_l[0]-self.params.Offset_L
        bin_l[1]=bin_l[1]-self.params.Offset_L
        bin_e[0]=bin_e[0]-self.params.Offset_E
        bin_e[2]=bin_e[2]-self.params.Offset_E

        self.Energy=RSE_Constants.rawData.Energy
        self.Intensity=RSE_Constants.rawData.Intensity
        self.Error=RSE_Constants.rawData.Error
        
        self.clean()
        FileIsGood=self.dataIsValid(minPoints)
        if FileIsGood:
            fileForSlice=DataTextFile(folder,file)
#            print ('mult  '+ str(mult))
            fileForSlice.Write(self.Energy,mult*self.Intensity,mult*self.Error)
            return 0
        return 1

    # ----------------------------------------------------------------------------------------------

    def dim2array(self,d):
        
        """
        Create a numpy array containing bin centers along the dimension d
        input: d - IMDDimension
        return: numpy array, from min+st/2 to max-st/2 with step st  
        """
        dmin=d.getMinimum()
        dmax=d.getMaximum()
        dstep=d.getX(1)-d.getX(0)
        return np.arange(dmin+dstep/2,dmax,dstep)

    # ----------------------------------------------------------------------------------------------

    def dataIsValid(self,minPointsInDataFile):

        nn=0
        FileIsGood=False
        for i in range (0,len(self.Energy)):
            if self.Intensity[i]>0:
                if self.Error[i]/self.Intensity[i]<self.params.ErrorToIntensityMaxRatio:
                   nn=nn+1
        if nn>minPointsInDataFile:
            FileIsGood=True

        return FileIsGood  

    # ----------------------------------------------------------------------------------------------

    def ExtractQfromFileName(self,FileName):

        h=FileName[1:FileName.find("K")]
        k=FileName[FileName.find("K")+1:FileName.find("L")]
        l=FileName[FileName.find("L")+1:21] 

        H=float(h);
        K=float(k);
        L=float(l);

        Q=[H,K,L]
        return Q

    def SubtractConstant(self,Const):
        for i in range(len(self.Intensity)):
            self.Intensity[i]=self.Intensity[i]-Const

    def SubtractLine(self,Intercept,Slope):
        for i in range(len(self.Intensity)):
            self.Intensity[i]=self.Intensity[i]-Intercept-Slope*self.Energy[i]

    def DivideByBoseFactorNorm(self,T,M=1.0):  #M normalization
        for i in range(0,len(self.Energy)):
           self.Intensity[i]=M*self.Intensity[i]/(1+1/(np.exp(self.Energy[i]/(0.08617*T))-1))
           self.Error[i]=M*self.Error[i]/(1+1/(np.exp(self.Energy[i]/(0.08617*T))-1))

    # ----------------------------------------------------------------------------------------------

    def Generate(self):

        """
        this is a wrapper to _generate.  it comes up with an array of Qpts to get data for, 
        then goes and gets the data
        """

        # get the Q-points
        self.get_Q_points_arr()
        self._generate()

    # ----------------------------------------------------------------------------------------------

    def _generate(self):

        """
        modified by T.S. Apr 2022

        generate the SQW data for the Qpts stored in self.Qpts 
        """

        bin_e=[self.params.e_start,self.params.e_step,self.params.e_end]        

        Qs = self.Qpts

        # get metadata for Q-pts
        try:
            QHlist=Qs[:,][:,0]
            QKlist=Qs[:,][:,1]
            QLlist=Qs[:,][:,2]
        except:
            QHlist=[Qs[0]]
            QKlist=[Qs[1]]
            QLlist=[Qs[2]]

        # matlab engine has to be restarted on every process! cant broadcast it in memory for some reason.
        self.fileHandle = self.initialize()

        # now go and cut data from files
        for i in range (0,len(QHlist)):

            fileName=str(RSE_Constants.FILENAME_FORMAT % (QHlist[i],QKlist[i],QLlist[i]))
            bin_h=[QHlist[i]-self.params.Deltah, QHlist[i]+self.params.Deltah]
            bin_k=[QKlist[i]-self.params.Deltak, QKlist[i]+self.params.Deltak]
            bin_l=[QLlist[i]-self.params.Deltal, QLlist[i]+self.params.Deltal]

            try:
                print(fileName)
                self.makeRawSlice(bin_h, bin_k, bin_l, bin_e, self.dataDirectory, \
                                    fileName, self.params.MinPointsInDataFile)

            except Exception as e:
                print(e)
                print("no slice CollectionOfQs")

    # ----------------------------------------------------------------------------------------------


# =======================================================================================================================
# -----------------------------------------------------------------------------------------------------------------------
# =======================================================================================================================


class DataSmall_q(Dataset):

    # ----------------------------------------------------------------------------------------------

    def __init__(self,params,dataDirectory,q=[1000,1000,1000]):

        self.params=params
        self.dataDirectory=dataDirectory
        self.filenames=[]
        self.q=q

    # ----------------------------------------------------------------------------------------------

    def Read(self):

        try:   #if folder with data slices is not there, just skip this

            self.filenames = self.Filterlist([file for file in os.listdir(self.dataDirectory) \
                    if file.startswith("H") and not file.endswith(".pdf")])
            
#            print(self.filenames)
            self.NumberofDatasets = len(self.filenames)
            self.dataset=Dataset(self.dataDirectory,self.filenames)
            self.Energy=self.dataset.Energy
            self.Intensity=self.dataset.Intensity
            self.Error=self.dataset.Error
#            print(self.NumberofDatasets)
        except:
            print("data 116, except");
#            return
        return
    
    # ----------------------------------------------------------------------------------------------

    def Filterlist(self,filenames):
        filteredFilenames=[]
        if self.q==[1000,1000,1000]:
            return filenames
        else:
           sys.path.insert(0,"reduced q Algorithms")
           Sq=__import__(self.params.SmallqAlgorithm)
           for i in range(0,len(filenames)):
               filename=filenames[i]
               Q=self.ExtractQfromFileName(filename)
               qq=Sq.Smallq(Q)
               q=qq.q

#               q=self.ConvertToSmallQ(Q)
               if (abs(q[0]-self.q[0])<0.0001 and abs(q[1]-self.q[1])<0.0001 and abs(q[2]-self.q[2])<0.0001):
                   filteredFilenames.append(filename)
                   
        return filteredFilenames               
    
    # ----------------------------------------------------------------------------------------------

    def get_Q_points_arr(self):

        """
        added by T.S. Apr 2022

        need to know what Qpoints to do in advance of looping over them. this is sensible in general but is
        required to do them in parallel

        in the case of 'DataSmall_q', i cut the qpt generation part of 'Generate' out and put it here. now
        we can use the same Generate method as 'CollectionOfQs'
        """

        _Qpts = []

        qh=self.params.qh
        qk=self.params.qk
        ql=self.params.ql
        sys.path.insert(0,"reduced q Algorithms")
        Sq=__import__(self.params.SmallqAlgorithm)
        qq=Sq.Smallq([qh,qk,ql])

        signVarL=[]
        signVarL=self.signVar(ql)
        signVarH=self.signVar(qh)
        signVarK=self.signVar(qk)
        for jjj in range(0,len(qq.qlist)):
            qh=qq.qlist[jjj][0]
            qk=qq.qlist[jjj][1]
            ql=qq.qlist[jjj][2]

            for l in range (self.params.l_start,self.params.l_end+1):
                for ii in range (0,len(signVarL)):
                    L=l+signVarL[ii]*ql

                    for h in range (self.params.h_start,self.params.h_end+1):
                        for jj in range (0,len(signVarH)):
                           H=h+signVarH[jj]*qh

                           for k in range (self.params.k_start,self.params.k_end+1):
                               for kk in range (0,len(signVarK)):
                                   K=k+signVarK[kk]*qk

                                   _Qpts.append([H,K,L])

        self.Qpts = np.array(_Qpts)
        self.num_Qpts = self.Qpts.shape[0]
        print(f'\n there are {self.num_Qpts} Qpts\n')

    # ----------------------------------------------------------------------------------------------

    def signVar(self,q):
        if q==0 or q==0.5:
            return [1]
        else:
            return [-1,1]





# =======================================================================================================================
# -----------------------------------------------------------------------------------------------------------------------
# =======================================================================================================================




class CollectionOfQs(Dataset):

    # ----------------------------------------------------------------------------------------------

    def __init__(self,params):

        self.params=params
        self.dataDirectory=self.params.path_data

    # ----------------------------------------------------------------------------------------------

    def get_Q_points_arr(self):

        """
        added by T.S. Apr 2022

        need to know what Qpoints to do in advance of looping over them. this is sensible in general but is
        required to do them in parallel

        in this case, just read the specified file
        """

        self.Qpts=np.genfromtxt(self.params.path_InputFiles+self.params.textfile_for_selectedQs)
        if self.Qpts.shape[1] != 3:
            exit('\n ERROR!\n  Qpts in file should have shape [N_Q]x[3]\n')
        self.num_Qpts = self.Qpts.shape[0]
        print(f'\n there are {self.num_Qpts} Qpts\n')

    # ----------------------------------------------------------------------------------------------






# =======================================================================================================================
# -----------------------------------------------------------------------------------------------------------------------
# =======================================================================================================================


class DataBackgroundQs(Dataset):

    # ----------------------------------------------------------------------------------------------

    def __init__(self,params):

        self.params=params
        self.params.ReadBackgroundParams() # gets input opts for background cuts
        self.dataDirectory=self.params.path_data

    # ---------------------------------------------------------------------------------------------

    def GenerateAllFiles(self):
       
        # loop over the extant Q-files
        dataFileNames=[file for file in os.listdir(self.dataDirectory) if file.startswith("H") \
                                                    and not file.endswith("pdf")]

        for ii in range(0,len(dataFileNames)):

            raw_data_file = dataFileNames[ii]
            raw_data=Dataset(self.dataDirectory,[raw_data_file]) # get the data in the file

            Q=self.ExtractQfromFileName(raw_data_file) # get the Q-pt from the file name
            
            # need these later
            self.raw_data_file = raw_data_file 
            self.raw_data = raw_data 
            self.Q =  Q

            # get the energy binning from the Q-pt data
            self.energy_start=raw_data.Energy[0]
            deltaEnergy=self.params.e_step
            lastEnergyIndex=len(raw_data.Energy)-1
            self.energy_end=raw_data.Energy[lastEnergyIndex]
            Offset=floor(self.params.Resolution/deltaEnergy)+1
            self.bin_e=[self.energy_start-Offset*deltaEnergy, \
                                    deltaEnergy,self.energy_end+Offset*deltaEnergy]

            # get the background data for this Q
            Q=self.ExtractQfromFileName(raw_data_file) # get the Q-pt from the file name
            self.GenerateFolder(Q) # now get the BG cuts

    # ----------------------------------------------------------------------------------------------

    def GenerateFolder(self,Q):

        self.H=Q[0]
        self.K=Q[1]
        self.L=Q[2]

        # get the path to the folder
        folder=self.dataDirectory+RSE_Constants.GEN_FOLDER + \
                                str(RSE_Constants.DIR_FORMAT%(self.H,self.K,self.L))

        # we want a ref for this later
        self.folder = folder
        
        # enforce that a BG folder exists
        if not os.path.isdir(folder):
            os.makedirs(folder)

        self.bin_h=[self.H-self.params.Deltah, self.H+self.params.Deltah]
        self.bin_k=[self.K-self.params.Deltak, self.K+self.params.Deltak]
        self.bin_l=[self.L-self.params.Deltal, self.L+self.params.Deltal]

        # ------------- added by Tyler Sterling May 2025 -----------------------
        # gets the data for the non-BG Q-pt, find energy window along the way
        #self.get_adaptive_energy_window() 
        # ----------------------------------------------------------------------

        # ------------- this stuff added by Tyler Sterling Apr 2022 ---------------------------
        self.prep_BG_cuts() # this gets the set of Qpts to do
        self.generate_BG_cuts() # this gets the data from the sqw/nxs file
        # ------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------------------

    def get_adaptive_energy_window(self):
        
        """
        added by Tyler Sterling May 2022

        this uses an 'adaptive' method to cut the data on a larger energy window. the window in bin_e
        is padded by an "offset" that is added to remove end effects from the Gaussian convolution

        the window is then checked by "adaptively" cutting the raw data, checking that its good, and if 
        not, removing some of the padding. 

        I think this is unnecessary ... in my BG algorithm, the data are subtracted only for the window
        that the raw data is cut on. they don't even have to have the same energy step since interpolation
        is done. if there are NANS at the ends, they are removed from the backgroud subtraction by setting
        thier value to something very large, thus removing those curves from the point-by-point minumum 
        method

        """

        # need to start a matlab engine 
        self.initialize()

        # file name for non-BG Qpt
        fileName=str(RSE_Constants.FILENAME_FORMAT % (self.H,self.K,self.L))

        dataFileNames=[file for file in os.listdir(self.folder) if \
                                    file.startswith("H") and not file.endswith("pdf")]

        try:

            msg = 'using bin_e with padded upper AND lower bound'
            # use full padded window; returns False if too much bad data
            if (self.makeRawSlice(self.bin_h, self.bin_k, self.bin_l, self.bin_e, self.folder, \
                            fileName, self.params.MinPointsInDataBackgroundFile)==1):
                
                # try with reduced upper bound
                msg = 'using bin_e with padded lower bound ONLY'
                self.bin_e=[self.bin_e[0],self.bin_e[1],self.energy_end]
                if (self.makeRawSlice(self.bin_h, self.bin_k, self.bin_l, self.bin_e, self.folder, \
                            fileName, self.params.MinPointsInDataBackgroundFile)==1):

                    # try with reduced lower bound
                    msg = 'using bin_e with padded upper bound ONLY'
                    self.bin_e=[self.energy_start,self.bin_e[1],self.bin_e[2]]
                    if (self.makeRawSlice(self.bin_h, self.bin_k, self.bin_l, self.bin_e, self.folder, \
                                fileName, self.params.MinPointsInDataBackgroundFile)==1 ):

                        # reduce both bounds to the min as last resort
                        msg = 'using bin_e from file'
                        self.bin_e=[self.energy_start,self.bin_e[1],self.energy_end]

            # now go and get the data with whatever bin_e was selected
            print('\n '+msg+'\n\n')
            self.makeRawSlice(self.bin_h, self.bin_k, self.bin_l, self.bin_e, self.folder, \
                        fileName, self.params.MinPointsInDataBackgroundFile)

        except Exception as e:
            print(e)
            print("no slice 1")

        print(f'\n maxfiles  {self.params.maxFiles}\n')


    # ----------------------------------------------------------------------------------------------

    def prep_BG_cuts(self):

        """
        added by Tyler Sterling Apr 2022

        get Q-points array, binning params, filenames, etc. we do it in advance to 
        enable parallelism
        """

        # use whatever tools were specified in params file
        sys.path.insert(0,"Background Tools")
        B = __import__(self.params.BackgroundAlgorithm)

        # hold all of the data needed to get the Q-points
        self.num_Qpts = 2*self.params.maxFiles+1 # there is a bonus one for raw data Q-pt
        self.Qpts = np.zeros((self.num_Qpts,3))
        self.bg_bin_h = np.zeros((self.num_Qpts,2))
        self.bg_bin_k = np.zeros((self.num_Qpts,2))
        self.bg_bin_l = np.zeros((self.num_Qpts,2))
        self.bg_file_names = []
        self.bg_mult = np.zeros(self.num_Qpts)

        # put the raw data Q-pt in the array for one of the procs to do; always the 0th index
        self.Qpts[0,:] = self.Q[:]
        self.bg_bin_h[0,:] = [self.Q[0]-self.params.Deltah, self.Q[0]+self.params.Deltah]
        self.bg_bin_k[0,:] = [self.Q[1]-self.params.Deltak, self.Q[1]+self.params.Deltak]
        self.bg_bin_l[0,:] = [self.Q[2]-self.params.Deltal, self.Q[2]+self.params.Deltal]
        self.bg_mult[0] = 1
        self.bg_file_names.append(self.raw_data_file)

        for ii in range(1,self.num_Qpts):

            BkgQ = B.BackgroundQ(self.H,self.K,self.L,self.params,ii) # generates a phi/theta pair

            if BkgQ.flag==1:
                return 3

            self.Qpts[ii,:] = BkgQ.Qslash[:]

            #CalcQslash takes input in reciprocal angstoms, result returned in r.l.u
            self.bg_bin_h[ii,:] = [BkgQ.Qslash[0]-self.params.Deltah, BkgQ.Qslash[0]+self.params.Deltah]
            self.bg_bin_k[ii,:] = [BkgQ.Qslash[1]-self.params.Deltak, BkgQ.Qslash[1]+self.params.Deltak]
            self.bg_bin_l[ii,:] = [BkgQ.Qslash[2]-self.params.Deltal, BkgQ.Qslash[2]+self.params.Deltal]
            self.bg_file_names.append(BkgQ.fileName)
            self.bg_mult[ii] = BkgQ.mult

    # ----------------------------------------------------------------------------------------------

    def generate_BG_cuts(self):

        """
        added by Tyler Sterling Apr 2022

        loop over the assigned Q-pts. at each iteration, exit if numFiles==maxFiles, i.e. if enough
        BG files have been produced. this enables multiple processes to work on the same directory at
        once without wasting effort.
        """

        # startsup the matlab engine
        self.initialize()

        folder = self.folder
        
        # the inds for this proc to do
        inds = list(range(self.num_Qpts))

        # the raw data file MUST be cut into the dir, so do it no matter what
        if 0 in inds: # the raw data file is always the 0th in the list of inds

            ind = inds.index(0)
            inds.pop(ind)

            bin_h = self.bg_bin_h[ind,:].tolist()
            bin_k = self.bg_bin_k[ind,:].tolist()
            bin_l = self.bg_bin_l[ind,:].tolist()
            file_name = self.bg_file_names[ind]
            mult = self.bg_mult[ind]

            try: # try to cut data from file
                print(folder+file_name)
                self.makeRawSlice(bin_h, bin_k, bin_l, self.bin_e,folder,file_name, \
                    self.params.MinPointsInDataBackgroundFile,mult)

            except Exception as e: # if garbage, return 1 and it will try again
                print(e)
                print("no slice 2")

        # start doing the rest of the BG cuts
        for ii in range(len(inds)):

            # --- get the number of files already present ---

            dataFileNames=[file for file in os.listdir(folder) if file.startswith("H") and not file.endswith("pdf")]
            numFiles=len(dataFileNames)

            print(f'\n numFiles: {numFiles}\n')

            # break the loop if there are enough files already

            if numFiles==self.params.maxFiles:
                break

            # --- otherwise go get the cuts ---

            ind = inds[ii]

            bin_h = self.bg_bin_h[ind,:].tolist()
            bin_k = self.bg_bin_k[ind,:].tolist()
            bin_l = self.bg_bin_l[ind,:].tolist()
            file_name = self.bg_file_names[ind]
            mult = self.bg_mult[ind]

            try: # try to cut data from file
                print(folder+file_name)
                self.makeRawSlice(bin_h, bin_k, bin_l, self.bin_e,folder,file_name, \
                    self.params.MinPointsInDataBackgroundFile,mult)

            except Exception as e: # if garbage, return 1 and it will try again
                print(e)
                print("no slice 2")

# =======================================================================================================================
# -----------------------------------------------------------------------------------------------------------------------
# =======================================================================================================================


