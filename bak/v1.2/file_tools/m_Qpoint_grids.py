
"""
Author: Tyler C. Sterling
Email: ty.sterling@colorado.edu
Affil: University of Colorado Boulder, Raman Spectroscopy and Neutron Scattering Lab
Date: 03/01/2024
Description:
    wraper for m_MDE_tools.py and MDEventCreator.py (custom script to merge data from JPARC 
    into MDE objects using mantid algorithms). m_MDE_tools was designed to take a single set
    of binnings and Q-points and write a file, and the Q-points were set by the bin spacing. 
    if different offset Q-point gride were desired, a seperate file was needed and when 
    analyzing data w/ phonon-explorer, if the different grids were used, the data file arg
    had to change. THIS WAS A MESS. 
    this wrapper is designed to bin data on numerous offset grids and put them all in the same
    file, to make life easier!
"""

# --------------------------------------------------------------------------------------------------

class c_Qpoint_grids:

    # ----------------------------------------------------------------------------------------------

    def __init__(self,label,units='s'):
        """
        small tool for timing and printing timing info
        """

# --------------------------------------------------------------------------------------------------



