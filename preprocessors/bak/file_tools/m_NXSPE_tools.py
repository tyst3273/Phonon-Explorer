
"""
PreprocessNXSPE function:
author: Andrei T. Savici
email: saviciat@ornl.gov

Everything else:
author: Tyler C. Sterling
Email: ty.sterling@colorado.edu
Affil: University of Colorado Boulder, Raman Spectroscopy and Neutron Scattering Lab

Date: 08/22/2024
Description:
    - parse files NXSPE files from JPARC, merge into MDE dataset, use MDNorm to bin, 
        then write the histogrammed data to .hdf5 file that can be used by phonon 
        explorer. 
    - all the code to talk to the JPARC files was written by Andrei. Tyler wrote the 
        wrapper that loops over Q-pt grids and he wrote the code that writes histo 
        data to file.
    - now allows user to pass other mantid key-word arguments to MDnorm, e.g. to pass
        symmetry args into MDNorm to create symmetrized hisotgram file.

    - note: to get mantid to correctly accept files from JPARC, you have to use a 
        modified "Facilities.xml" file. there is one provied in 'file_tools' dir. 
        allegedly we can just copy this to our $HOME/.mantid/instrument dir, but that 
        didnt work recently (it used to). I had to edit $HOME/.mantid/Mantid.user.properties 
        file and add instrumentDefinition.directory=$HOME/.mantid/instrument where you should 
        replace $HOME by the actual path. this might break some stuff, however, as mantid 
        needs to find the correct *.xml configuration files for each instrument. i found these
        under /opt/mantid*/instrument where the * means the correct version. i copied them 
        to $HOME/.mantid/instrument and all is well! :)
    - update: this broke some other mantid functionality later and i had to undo it, but 
        this script is still work with the unmodified .xml file! 
"""

import numpy as np
from mantid.simpleapi import *

from file_tools.m_save_MDE_to_hdf5 import save_MDE_to_hdf5
from file_tools.m_file_utils import crash, c_timer

"""
### DEV
# reload the modules if running interactively while modifying modules
from importlib import reload
import file_tools.m_save_MDE_to_hdf5 as m_save_MDE_to_hdf5
reload(m_save_MDE_to_hdf5)
save_MDE_to_hdf5 = m_save_MDE_to_hdf5.save_MDE_to_hdf5
### DEV
"""

# --------------------------------------------------------------------------------------------------

def PreprocessNXSPE(H_bins,K_bins,L_bins,event_files,goniometer_file,UB_params,**MDNorm_kw_args):
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

    """
    # !!! DEV !!!
    # set _stride != 0 to use only every _stride^th files to test run
    _stride = 100
    if _stride:
        event_files = event_files[::_stride]  
    # !!! DEV !!!
    """

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

    print('MDNorm keyword-args:')
    print(MDNorm_kw_args,'\n')

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
           OutputWorkspace='o',
           OutputDataWorkspace='d',
           OutputNormalizationWorkspace='n',
           **MDNorm_kw_args)

    #SaveMD(InputWorkspace="o",Filename="/SNS/ARCS/IPTS-26347/shared/jpark/MDHistoFile.nxs")

    MD_workspace = mtd['o']

    # to save data from MDNorm, we just need the output MD_workspace
    return MD_workspace

    # ----------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------

def _setup_shifts(H_step,K_step,L_step,H_bin,K_bin,L_bin):
    """
    """
    _H_shifts = np.arange(H_bin)*H_step
    _K_shifts = np.arange(K_bin)*K_step
    _L_shifts = np.arange(L_bin)*L_step

    print('H_shifts:',_H_shifts)
    print('K_shifts:',_K_shifts)
    print('L_shifts:',_L_shifts,'\n')

    _H_shifts, _K_shifts, _L_shifts = np.meshgrid(_H_shifts,_K_shifts,_L_shifts,indexing='ij')
    _H_shifts = _H_shifts.flatten()
    _K_shifts = _K_shifts.flatten()
    _L_shifts = _L_shifts.flatten()
    shifts = np.c_[_H_shifts,_K_shifts,_L_shifts]

    ind = np.argwhere((shifts == 0.0).all(axis=1))
    shifts = np.delete(shifts,ind,axis=0)

    return np.atleast_2d(shifts)

# --------------------------------------------------------------------------------------------------

def _get_bins(lo,hi,step):
    """
    """
    lo = float(lo); hi = float(hi); step = float(step)
    _mod = np.modf((hi-lo)/step)
    num_bins = _mod[1].astype(int)
    decimal = _mod[0].round(6)
    if decimal != 0.0:
        msg = '\n*** WARNING ***\n'
        msg += f'upper bin center {hi: 6.3f} is not commensurate with\n'
        msg += f'lower bin center {lo: 6.3f} and bin size {step:6.3f}\n'
        hi = lo+step*(num_bins+1)
        msg += f'tweaking upper bin center to {hi: 6.3f}\n'
        print(msg)
    lo = np.round(lo-step/2,6)
    hi = np.round(hi+step/2,6)
    step = np.round(step,6)
    bins = [lo,step,hi]
    return bins

# --------------------------------------------------------------------------------------------------

def bin_NXSPE(event_files,goniometer_file,output_file,UB_params,H_lo,H_hi,H_step,K_lo,K_hi,K_step,
            L_lo,L_hi,L_step,H_bin=1,K_bin=1,L_bin=1,**MDNorm_kw_args):

    """
    wrapper to translate variable names for interface with PreprocessNXSPE. note, E-binning is 
    fixed by the binning already present in the *.nxspe files

    MDNorm_kw_args are any bonus arguments you want to pass to MDNorm. Note, you have to give 
    them as kw-args in function call:
        ex. bin_NXSPE(...,SymmetryOperations='x,y,z;-x-y-z')

    """

    timer = c_timer('bin_NXSPE',units='m')

    H_bin = int(H_bin); K_bin = int(K_bin); L_bin = int(L_bin)
    if H_bin < 1: H_bin = 1
    if K_bin < 1: K_bin = 1
    if L_bin < 1: L_bin = 1

    # convert binning args to ones used by c_MDE_tools
    H_bins = _get_bins(H_lo,H_hi,H_step*H_bin)
    K_bins = _get_bins(K_lo,K_hi,K_step*K_bin)
    L_bins = _get_bins(L_lo,L_hi,L_step*L_bin)
    print('H_bins:',H_bins)
    print('K_bins:',K_bins)
    print('L_bins:',L_bins,'\n')

    # convert shift args to useful form
    shifts = _setup_shifts(H_step,K_step,L_step,H_bin,K_bin,L_bin)
    if shifts.size == 0:
        loop_over_shifts = False
    else:
        loop_over_shifts = True

    # bin and save data with no offset
    MD_workspace = PreprocessNXSPE(H_bins,K_bins,L_bins,event_files,goniometer_file,UB_params,
                                   **MDNorm_kw_args)
    save_MDE_to_hdf5(MD_workspace,output_file)

    # only loop over offsets if atleast one is defined
    if loop_over_shifts:

        num_shifts = shifts.shape[0]

        print('\n----------------------------------------------------------------\n')
        print('num_shifts:',num_shifts)
        print('shifts:')
        print(shifts)

        msg = '\nlooping over shifts'
        print(msg)

        _H_shift = np.array([0,0,0],dtype=float)
        _K_shift = np.array([0,0,0],dtype=float)
        _L_shift = np.array([0,0,0],dtype=float)

        for ii in range(num_shifts):

            _H, _K, _L = shifts[ii,:]
            print(f'\nshifts[{ii+1}]:',_H,_K,_L,'\n')

            _H_shift[0] = _H; _H_shift[2] = _H
            _K_shift[0] = _K; _K_shift[2] = _K
            _L_shift[0] = _L; _L_shift[2] = _L

            # go and bin the data -- these MUST be appended.
            MD_workspace = PreprocessNXSPE(H_bins+_H_shift,K_bins+_K_shift,L_bins+_L_shift,
                                           event_files,goniometer_file,UB_params,**MDNorm_kw_args)
            save_MDE_to_hdf5(MD_workspace,output_file,append=True)

    timer.stop()

# --------------------------------------------------------------------------------------------------    


