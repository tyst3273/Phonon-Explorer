
"""
Author: Andrei T. Savici
Email: saviciat@ornl.gov

Updated by: Tyler C. Sterling
Email: ty.sterling@colorado.edu
Affil: University of Colorado Boulder, Raman Spectroscopy and Neutron Scattering Lab

Date: 03/05/2024
Description:
    - parse files NXSPE files from JPARC, merge into MDE dataset, use MDNorm to bin, 
        then write the histogrammed data to .hdf5 file that can be used by phonon 
        explorer. 
    - all the code to talk to the JPARC files was written by Andrei. Tyler wrote the 
        wrapper that loops over Q-pt grids and he wrote the code that writes histo 
        data to file.

    - note: to get mantid to correctly accept files from JPARC, you have to use a 
        modified "Facilities.xml" file. there is one provied in 'file_tools' dir. 
        allegedly we can just copy this to our $HOME/.mantid/instrument dir, but that 
        didnt work recently (it used to). I had to edit $HOME/.mantid/Mantid.user.properties 
        file and add instrumentDefinition.directory=/SNS/users/YOUR-USER-NAME/.mantid/instrument
        where you should replace YOUR-USER-NAME by your user name on the SNS server. 
        for other computers, the path should ultimately it should point to 
        $HOME/.mantid/instrument
"""

import numpy as np
from mantid.simpleapi import *

from file_tools.m_save_MDE_to_hdf5 import save_MDE_to_hdf5
from file_tools.m_file_utils import crash

# --------------------------------------------------------------------------------------------------

def PreprocessNXSPE(H_bins,K_bins,L_bins,event_files,goniometer_file,UB_params):
    """
    wrap the code to load and preprocess NXSPE files from JPARC measurement.

    the code below was written by Andrei T. Savici at ORNL (saviciat@ornl.gov). I have no idea
    what it does and don't know how to modify it. I am simply wrapping it in a function to 
    programmatically pass binning args in and get the data out.
    """

    msg = '\n*** WARNING ***\n'
    msg += 'PreprocessNXSPE is currently configured for the manganite experiment\n'
    msg += 'at JPARC! It is not generic and must be redone for other experiments!\n'
    print(msg)

    d = np.loadtxt(goniometer_file,skiprows=1,dtype=str)
    runs = d[:,0]
    gon = d[:,11].astype(float)-82.388
    gon_dict = dict()
    for r,g in zip(runs, gon):
        gon_dict[r] = g

    event_files = ','.join(event_files)
    wg = Load(event_files)
    for wsn in wg.getNames():
        ang = gon_dict[wsn.split('S0')[1]]
        SetGoniometer(Workspace=wsn, Axis0 = f'{ang},0,1,0,1')

    axis_deltaE = wg[0].readX(0)
    emin = axis_deltaE[0]
    emax = axis_deltaE[-1]
    de_step = axis_deltaE[1]-emin

    ConvertFromDistribution(wg)
    wge = ConvertToEventWorkspace(wg)
    wge = CropWorkspaceForMDNorm(wge,XMin=emin,XMax=emax)

    SetUB(wge,**UB_params)

    AddSampleLog(Workspace='wge', LogName='gd_prtn_chrg', LogText='1.0', LogType='Number',
        NumberType='Double')

    mdparts = ConvertToMD(InputWorkspace=wge, QDimensions='Q3D', dEAnalysisMode='Direct',
        Q3DFrames='Q_sample')

    mde = MergeMD(mdparts)
    #SaveMD(InputWorkspace="mde",Filename="/SNS/ARCS/IPTS-26347/shared/jpark/MDEventFile.nxs")

    #UB_params = {'a':3.93,'b':3.93,'c':3.93,'alpha':90,'beta':90,'gamma':90,
    #    'u':'1,0,-0.1','v':'0,1,-0.1'}
    #SetUB(mde,**UB_params)

    #print(axis_deltaE)
    # test MDNorm - elastic slice
    MDNorm(InputWorkspace=mde,
           QDimension0='1,0,0',
           QDimension1='0,1,0',
           QDimension2='0,0,1',
           Dimension0Name='QDimension0',
           Dimension0Binning="{},{},{}".format(H_bins[0],H_bins[1],H_bins[2]),
           Dimension1Name='QDimension1',
           Dimension1Binning="{},{},{}".format(K_bins[0],K_bins[1],K_bins[2]),
           Dimension2Name='QDimension2',
           Dimension2Binning="{},{},{}".format(L_bins[0],L_bins[1],L_bins[2]),
           Dimension3Name='DeltaE',
           Dimension3Binning="{}".format(de_step),
    #       SymmetryOperations='x,y,z;-x,-y,z;x,y,-z;-x,-y,-z;y,x,z;y,x,-z',
           OutputWorkspace='o',
           OutputDataWorkspace='d',
           OutputNormalizationWorkspace='n')

    #SaveMD(InputWorkspace="o",Filename="/SNS/ARCS/IPTS-26347/shared/jpark/MDHistoFile.nxs")

    MD_workspace = mtd['o']

    # to save data from MDNorm, we just need the output MD_workspace
    return MD_workspace

    # ----------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------

def bin_NXSPE_with_offsets(event_files,goniometer_file,output_file,UB_params,
    H_bins,K_bins,L_bins,H_offsets=None,K_offsets=None,L_offsets=None):
    """
    bin the events in the MDE workspace into a histogram workspace; note that this doesnt
    really return anything or create new attributes. the produced data are stored in the
    histogram workspace.

    turns out that having a seperate file for every grid offset is a nuisance. we want to 
    have multiple binnings/grid offsets in a single file to make our lives easier.

    H_bin = [start, step, stop]. note, these are now interpreted like the args to mantid
        MDNorm and to Horace cut_sqw: the binning will actually start at 'start' with spacing
        equal to the step size. i.e. the bin centers will be [start+1*step/2, start+2*step/2,
        start+3*step/2, ...]. similarly for K_bins, L_bins.

    H_offets = [offset_1, offset_2, ..., offset_N]. loop over the offets in H_offsets and bin
        using the args in H_bin shifted by the offets. e.g. for H_offset = [0.0, 0.01] and
        H_bins = [-0.05, 0.1, 0.55], we will generate data on a grid intergrated around
        H=0.0+-0.05, H=0.1+-0.05, ..., H=0.5+-0.05 and similarly H=0.01+-0.05, H=0.11+-0.05, ...,   
        H=0.51+-0.05 etc. similarly for K_offets, L_offsets.

    num_chunks is number of chunks to split binning in Q-space. E.g. if num_chunks = [2,1,1]
    and there are 10 bins along H, the data will be cut in 2 go's with the first 5 H-bins in
    the first go and the last 5 in the seconds. 

    merged_file_name is the name of the output data. each offset grid is appended to the file

    u, v, w are the projections. w is optional and the cross product of u and v if not given.
    """

    H_bins = np.array(H_bins,dtype=float)
    K_bins = np.array(K_bins,dtype=float)
    L_bins = np.array(L_bins,dtype=float)

    # check if any offsets are defined and initialize null offsets otherwise
    loop_over_offsets = False
    if H_offsets is None:
        H_offsets = np.array([0.0],dtype=float)
    else:
        H_offsets = np.array(H_offsets,dtype=float)
        loop_over_offsets = True
    if K_offsets is None:
        K_offsets = np.array([0.0],dtype=float)
    else:
        K_offsets = np.array(K_offsets,dtype=float)
        loop_over_offsets = True
    if L_offsets is None:
        L_offsets = np.array([0.0],dtype=float)
    else:
        L_offsets = np.array(L_offsets,dtype=float)
        loop_over_offsets = True

    print('H_bins:',H_bins)
    print('K_bins:',K_bins)
    print('L_bins:',L_bins,'\n')

    # bin the data on the unshifted grid first. this initiates file etc.
    MD_workspace = PreprocessNXSPE(H_bins,K_bins,L_bins,event_files,goniometer_file,UB_params)
    save_MDE_to_hdf5(MD_workspace,output_file)

    # only loop over offsets if one is defined
    if loop_over_offsets:

        num_offsets = H_offsets.size*K_offsets.size*L_offsets.size

        print('\n----------------------------------------------------------------\n')
        print('H_offsets:',H_offsets)
        print('K_offsets:',K_offsets)
        print('L_offsets:',L_offsets)
        print('num_offsets:',num_offsets)

        msg = '\nlooping over offsets'
        print(msg)

        _H_offset = np.array([0,0,0],dtype=float)
        _K_offset = np.array([0,0,0],dtype=float)
        _L_offset = np.array([0,0,0],dtype=float)

        _n = 0

        # loop over the offsets and bin the data at each offset
        for ii, _H in enumerate(H_offsets):
            for jj, _K in enumerate(K_offsets):
                for kk, _L in enumerate(L_offsets):

                    _H_offset[0] = _H; _H_offset[2] = _H
                    _K_offset[0] = _K; _K_offset[2] = _K
                    _L_offset[0] = _L; _L_offset[2] = _L

                    print(f'\noffset[{_n+1}]:',_H,_K,_L,'\n')
                    _n += 1

                    # go and bin the data 
                    MD_workspace = PreprocessNXSPE(H_bins+_H_offset,K_bins+_K_offset,
                            L_bins+_L_offset,event_files,goniometer_file,UB_params)
                    save_MDE_to_hdf5(MD_workspace,output_file,overwrite=False)

# --------------------------------------------------------------------------------------------------
