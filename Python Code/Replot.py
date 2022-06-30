#________________________________________________________________________%
#              University of Colorado Boulder                            %
#                                                                        %
#              Dmitry Reznik, Irada Ahmadova, Aaron Sokolik                             %       
#                                                                        %
#    Work supported by the DOE, Office of Basic Energy Sciences,         %
#    Office of Science, under Contract No. DE-SC0006939                  %                  
#________________________________________________________________________%

# ------------------------------------------------------------------------
# remake plots without having to cut data again
# ------------------------------------------------------------------------

from TextFile import *
from RSE_Constants import *
from plotDataWithFit import *
from Display import *

# --------------------------------------------------------------------------------------------------

params=Parameters(RSE_Constants.INPUTS_PATH, RSE_Constants.INPUTS_FILENAME)

plot=Plot(params)
plot.Plot()
Disp=Display()
Disp.MakePlotSummary(params.path_data,params.ProcessedDataName)

# --------------------------------------------------------------------------------------------------


