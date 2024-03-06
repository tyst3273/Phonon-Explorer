
"""
Author: Tyler C. Sterling
Email: ty.sterling@colorado.edu
Affil: University of Colorado Boulder, Raman Spectroscopy and Neutron Scattering Lab
Date: 03/01/2024
Description:
    turns out that having a seperate file for every binning is a nuisance. we want to have 
    multiple binnings/grid offsets in a single file to make our lives easier
"""

from timeit import default_timer
import os

class c_Q_point_grids:

    # ----------------------------------------------------------------------------------------------

    def __init__(self,H_bins=None,K_bins=None,L_bins=None,E_bins=None,
            H_offsets=[0.0],K_offsets=[0.0],L_offsets=[0.0],
            histo_file_name='MD_histo.hdf5',u=[1,0,0],v=[0,1,0],w=None):
        """
        bin the events in the MDE workspace into a histogram workspace; note that this doesnt
        really return anything or create new attributes. the produced data are stored in the
        histogram workspace

        H_bin = [start, step, stop]. note, these are now interpreted like the args to mantid
            MDNorm and to Horace cut_sqw: the binning will actually start at 'start' with spacing
            equal to the step size. i.e. the bin centers will be [start+1*step/2, start+2*step/2,
            start+3*step/2, ...]. similarly for K_bins, L_bins.

        H_offets = [offset_1, offset_2, ..., offset_N]. loop over the offets in H_offsets and bin
        using the args in H_bin shifted by the offets. e.g. for H_offset = [0.0, 0.01] and 
        H_bins = [-0.05, 0.1, 0.55], we will generate data on a grid intergrated around 
        H=0.0+-0.05, H=0.1+-0.05, ..., H=0.5+-0.05 and similarly H=0.01+-0.05, H=0.11+-0.05, ..., 
        H=0.51+-0.05 etc. similarly for K_offets, L_offsets.

        hist_file_name is the name of the output data. each offset grid is appended to the file

        u, v, w are the projections. w is optional and the cross product of u and v if not given.
        """

    
    # ----------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------



