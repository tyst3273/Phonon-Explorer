#________________________________________________________________________%
#              University of Colorado Boulder                            %
#                                                                        %
#              Dmitry Reznik, Irada Ahmadova, Aaron Sokolik                             %       
#                                                                        %
#    Work supported by the DOE, Office of Basic Energy Sciences,         %
#    Office of Science, under Contract No. DE-SC0006939                  %                  
#________________________________________________________________________%

"""

Modified by Tyler Sterling, Apr. 2023

Whats new: The binning args. and actual Q-point found in the file are attached as attributes
    to this class. The binning args are self.Delta* and self.e_step and the Q-point is
    self.Qpoint_rlu. These here since they are specific to the data that is cut.

"""

import math
import sys
import os
import io
import numpy as np  # !!!!
from RSE_Constants import *
#import matlab.engine

class RawData: 

    """
    Dataset can be either a single cut at one Q or several cuts put together for the purposes of multizone fitting,
    """
    
    # ----------------------------------------------------------------------------------------------

    def __init__(self, params):

        import matlab.engine # -- better to import this once when module is imported (T.S.)

        RSE_Constants.fileHandle=matlab.engine.start_matlab()
        RSE_Constants.fileHandle.addpath(os.getcwd()+'/Tools to access raw data/') 
        self.sqw_path = params.sqw_path 

        # the binning is tracked by the RawData() classes
        self.Deltah = params.Deltah
        self.Deltak = params.Deltak
        self.Deltal = params.Deltal
        self.e_step = params.e_step

        # this will only run if params has HoracePath attr, which is only true for linux
        if hasattr(params,'HoracePath'):
            HoracePath = params.HoracePath
            RSE_Constants.fileHandle.addpath(HoracePath)
            RSE_Constants.fileHandle.herbert_on()
            RSE_Constants.fileHandle.horace_on()

        """
        # old method to do above code
        try: #This is needed only for Linux. Windows and Mac will throw an exception, which is safely caught.
            dataFile = params.sqw_path
            HoracePath = params.HoracePath
            RSE_Constants.fileHandle.addpath(HoracePath)
            RSE_Constants.fileHandle.herbert_on()
            RSE_Constants.fileHandle.horace_on()
        except Exception as e:
            print("IGNORE THIS ERROR UNLESS RUNNING ON LINUX:", e)
        """

        return

    # ----------------------------------------------------------------------------------------------

    def GetSlice(self, bin_h, bin_k, bin_l, bin_e, Projection_u, Projection_v):

        """
        #Projection_u , v are passed as strings
        # RETURNS 1 if slice is not written, 0 if SUCCESS
        """

        #eng = matlab.engine.start_matlab();
        import matlab.engine # -- better to import this once when module is imported (T.S.

        # attach this as attribute
        self.Qpoint_rlu = [np.array(bin_h).mean(),np.array(bin_k).mean(),np.array(bin_l).mean()]
        self.Qpoint_rlu = np.array(self.Qpoint_rlu)
        self.Projection_u = Projection_u
        self.Projection_v = Projection_v

        proj = dict(
            u=matlab.double(list(Projection_u)),
            v=matlab.double(list(Projection_v)),
            type='rrr',
            uoffset=matlab.double([0,0,0,0])
            )
        try:
            import StringIO
            out = StringIO.StringIO()
            err = StringIO.StringIO()
        except:
            out = io.StringIO()
            err = io.StringIO()

        self.Energy=[]
        self.Intensity=[]
        self.Error=[]

        bin_h_m = matlab.double([bin_h])
        bin_k_m = matlab.double([bin_k])
        bin_l_m = matlab.double([bin_l])
        bin_e_m = matlab.double([bin_e])

        try:
            ourCut = RSE_Constants.fileHandle.Getslice( 
                self.sqw_path,  #### !!!!
                proj,
                bin_h_m,
                bin_k_m,
                bin_l_m,
                bin_e_m,
                '-nopix',
               stdout=out,
               stderr=err)
        except Exception as e:
            print("I am here")
            print(e)
            return 1
    
        print("ourCut done")
        print("cut done")
        print(out.getvalue())
        print(err.getvalue())

        NumPoints=len(ourCut['p'])
        self.Intensity=np.zeros(NumPoints) # !!!!
        self.Error=np.zeros(NumPoints) # !!!!
        self.Energy=np.zeros(NumPoints) # !!!!

        for i in range(0,NumPoints-1):
            self.Energy[i]=ourCut['p'][i][0]
            self.Intensity[i]=ourCut['s'][i][0]
            self.Error[i]=ourCut['e'][i][0]

        return 0

        # ----------------------------------------------------------------------------------------------

