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

from TextFile import *
import math
import sys
import numpy
from numpy import *
from RSE_Constants import *
from mantid.simpleapi import LoadMD, MDNorm, CutMD, mtd, ConvertMDHistoToMatrixWorkspace, SaveAscii
#, ConvertToMD, BinMD, ConvertUnits, Rebin
from mantid.api import Projection

# sys.path.append('C:/MantidInstall/bin') this is probably not needed ...

class RawData: 

    """
    Dataset can be either a single cut at one Q or several cuts put together for the purposes of multizone fitting
    """

    # ----------------------------------------------------------------------------------------------

    def __init__(self,params):

        self.params = params
        RSE_Constants.fileHandle = LoadMD(params.sqw_path,OutputWorkspace='proj_md',FileBackEnd=True)

        # the binning is tracked by the RawData() classes
        self.Deltah = self.params.Deltah
        self.Deltak = self.params.Deltak
        self.Deltal = self.params.Deltal
        self.e_step = self.params.e_step

        return

    # ----------------------------------------------------------------------------------------------

    def GetSlice(self, bin_h, bin_k, bin_l, bin_e, Projection_u, Projection_v):
        
        # attach this as attribute
        self.Qpoint_rlu = [np.array(bin_h).mean(),np.array(bin_k).mean(),np.array(bin_l).mean()]
        self.Qpoint_rlu = np.array(self.Qpoint_rlu)

        print("bin_h: [{:01.3f}, {:01.3f}]".format(bin_h[0], bin_h[1]))
        print("bin_k: [{:01.3f}, {:01.3f}]".format(bin_k[0], bin_k[1]))
        print("bin_l: [{:01.3f}, {:01.3f}]".format(bin_l[0], bin_l[1]))
        print("bin_e: [{:01.3f}, {:01.3f}, {:01.3f}]".format(bin_e[0], bin_e[1], bin_e[2]))
        ProjVec_u=list(Projection_u)
        ProjVec_v=list(Projection_v)
        ProjVec_p=list(numpy.cross(ProjVec_u, ProjVec_v))
        Projection_p=str(ProjVec_p[0])+','+str(ProjVec_p[1])+','+str(ProjVec_p[2])
        projection={'QDimension0':Projection_u,'QDimension1':Projection_v,'QDimension2':Projection_p}
#        MDNorm(RSE_Constants.fileHandle,# InputWorkspace='proj_md',#SolidAngleWorkspace='van138487', 
        MDNorm('proj_md',# InputWorkspace='proj_md',#SolidAngleWorkspace='van138487', 
        Dimension0Name='QDimension0',Dimension0Binning="{},{}".format(bin_h[0], bin_h[1]), 
        Dimension1Name='QDimension1',Dimension1Binning="{},{}".format(bin_k[0], bin_k[1]), 
        Dimension2Name='QDimension2',Dimension2Binning="{},{}".format(bin_l[0], bin_l[1]), 
        Dimension3Name='DeltaE', Dimension3Binning="{},{},{}".format(bin_e[0], bin_e[1], bin_e[2]), 
        OutputWorkspace='proj_ws', OutputDataWorkspace='jd', OutputNormalizationWorkspace='jn',**projection)
        proj_ws_2d_h = ConvertMDHistoToMatrixWorkspace(InputWorkspace='proj_ws', OutputWorkspace='proj_ws_2d')
        E_edges = proj_ws_2d_h.extractX() # Bin edges
        E_centers = (E_edges[0][:-1] + E_edges[0][1:])/2
        I = 1000 * proj_ws_2d_h.extractY()[0]
        e = 1000 * proj_ws_2d_h.extractE()[0]
        self.Error = e
        self.Intensity = I
        self.Energy = E_edges[0][:-1]

        return 1

    # ----------------------------------------------------------------------------------------------
    
