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
#import matlab.engine
from RSE_Constants import *
from plotDataWithFit import *
from Display import *
import os
from Utils import timer, arg_parser

# --------------------------------------------------------------------------------------------------
# initialize stuff 

_args = arg_parser()
_args.get_input_file()
print(_args.input_file)
print(_args.input_path)


params=Parameters(RSE_Constants.INPUTS_PATH, RSE_Constants.INPUTS_FILENAME)

if params.QMode==0: # generate data on uniform grid of Q-points
    testData=DataSmall_q(params, params.path_data)
if params.QMode==1: # generate data on set of Q-points read from file
    testData=CollectionOfQs(params)


# --------------------------------------------------------------------------------------------------
# generate the data
_t = timer('const. Q cuts')

print('\n now getting const. Q cuts from file.\n  this might take a while...\n\n')
testData.Generate()

# if not getting background cuts, plot data now
if params.BkgMode==0:
    plot=Plot(params)
    plot.Plot()
    Disp=Display()
    Disp.MakePlotSummary(params.path_data,params.ProcessedDataName)

_t.stop()


# --------------------------------------------------------------------------------------------------
# otherwise, get the background cuts

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


