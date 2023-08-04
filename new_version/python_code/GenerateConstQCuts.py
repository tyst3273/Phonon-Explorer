#________________________________________________________________________%
#              University of Colorado Boulder                            %
#                                                                        %
#              Dmitry Reznik, Irada Ahmadova, Aaron Sokolik                             %       
#                                                                        %
#    Work supported by the DOE, Office of Basic Energy Sciences,         %
#    Office of Science, under Contract No. DE-SC0006939                  %                  
#________________________________________________________________________%

# ------------------------------------------------------------------------
# parallelism added to Data.py by Tyler Sterling Apr. 2022
# ------------------------------------------------------------------------

from TextFile import *
from Data import *
from RSE_Constants import *
from plotDataWithFit import *
from Display import *
import os
from Utils import timer
from SubtractPolyBackgr import *


# --------------------------------------------------------------------------------------------------
# initialize stuff 

params=Parameters()

# generate data on uniform grid of Q-points
if params.QMode==0: 
    testData=DataSmall_q(params, params.path_data)

# generate data on set of Q-points read from file
if params.QMode==1: 
    testData=CollectionOfQs(params)


# --------------------------------------------------------------------------------------------------
# generate the data

_t = timer('const. Q cuts')

testData.Generate()

# if not getting background cuts, plot data now
if params.BkgMode==0:
    plot=Plot(params)
    plot.Plot()
    Disp=Display()
    Disp.MakePlotSummary(params.path_data,params.ProcessedDataName)

_t.stop()

# --------------------------------------------------------------------------------------------------
# get the background cuts

if params.BkgMode==1:

    _t = timer('background cuts')
    start_time = timeit.default_timer()

    print('\n\n-------------------------------------------------------\n'+ \
        ' now getting BG cuts for each Q-pt.\n  this might take a while...\n\n')
    
    RSE_Constants.FLAG=0
    BackgroundFiles=DataBackgroundQs(params)
    BackgroundFiles.GenerateAllFiles()

    _t.stop()

# --------------------------------------------------------------------------------------------------

if params.BkgMode==2:
    SubtractPolyBackground()


