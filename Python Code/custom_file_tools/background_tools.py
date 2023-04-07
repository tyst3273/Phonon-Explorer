
from timeit import default_timer
import numpy as np
import h5py
import shutil
import os
from scipy.signal import convolve
import matplotlib.pyplot as plt


# --------------------------------------------------------------------------------------------------

def crash(err_msg=None,exception=None):
    """
    stop execution in a safe way
    """
    msg = '\n*** error ***\n'
    if err_msg is not None:
        msg += err_msg+'\n'
    if exception is not None:
        msg += '\nException:\n'+str(exception)+'\n'
    print(msg)
    raise KeyboardInterrupt

# --------------------------------------------------------------------------------------------------

def check_file(file_name):
    """
    check if the specified file exists
    """
    if not os.path.exists(file_name):
        msg = f'the file:\n  \'{file_name}\' \nwas not found!'
        crash(msg)

# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------

class c_timer:

    # ----------------------------------------------------------------------------------------------

    def __init__(self,label,units='s'):
        """
        small tool for timing and printing timing info
        """

        self.label = label
        if units == 'm':
            self.units = 'm'
            self.scale = 1/60
        elif units == 'ms':
            self.units = 'ms'
            self.scale = 1000
        else:
            self.units = 's'
            self.scale = 1
        self.start_time = default_timer()

    # ----------------------------------------------------------------------------------------------

    def stop(self):
        """
        stop timer and print timing info
        """

        elapsed_time = default_timer()-self.start_time
        elapsed_time *= self.scale
        msg = f'timing:   {self.label} {elapsed_time:9.5f} [{self.units}]'
        print(msg)

# --------------------------------------------------------------------------------------------------



class c_background_tools:

    # ----------------------------------------------------------------------------------------------

    def __init__(self):

        """
        class that reads raw file, smooths the data in it, and then subtracts the background from
        Q-points in the raw file by doing "rocking scans" and taking the minimum of each. 
        """

        pass

    # ----------------------------------------------------------------------------------------------

    def make_smoothed_file(self,raw_file,smoothed_file,smoothing_fwhm=None,num_blocks=1):

        """
        split the Q-point set in raw file in num_blocks blocks. for each, block
        interpolate the data onto a fine energy grid and then Gaussian smooth the data. 
        """

        check_file(raw_file)

        with h5py.File(raw_file,'r') as raw_db, h5py.File(smoothed_file,'w') as smooth_db:
           
            self.num_Q = raw_db['H_rlu'].size
            self.energy = raw_db['DeltaE'][...]

            if smoothing_fwhm is None:
                self.smoothing_fwhm = (self.energy[1]-self.energy[0])*3.0
            else:
                self.smoothing_fwhm = float(smoothing_fwhm)

            self._create_smooth_datasets(raw_db,smooth_db)

            # go and do the smoothing 
            self._smooth_by_looping_over_Q_blocks(num_blocks,raw_db,smooth_db)

    # ----------------------------------------------------------------------------------------------

    def _smooth_by_looping_over_Q_blocks(self,num_blocks,raw_db,smooth_db):

        """
        split Q into blocks and loop over the blocks, smoothing the Q-point in each block
        """

        _t = c_timer('loop_over_Q_blocks')

        Q_index_blocks = np.arange(self.num_Q)
        Q_index_blocks = np.array_split(Q_index_blocks,num_blocks)

        print(f'\nsplitting data into {num_blocks} blocks\n')

        for block in range(num_blocks):

            _bt = c_timer(f'block[{block}]')
            
            Q_indices = Q_index_blocks[block]
            signal = raw_db['signal'][Q_indices]
            smooth_db['smoothed_signal'][Q_indices,:] = \
                    self._gaussian_smooth(self.energy,signal,self.smoothing_fwhm)

            _bt.stop()

        _t.stop()

    # ----------------------------------------------------------------------------------------------

    def _create_smooth_datasets(self,raw_db,smooth_db):

        """
        copy 'header' info from raw_db to smooth_db, e.g. Qpts, lattice_vectors, etc.
        """

        _t = c_timer('create_smooth_datasets')

        print('\ncreating datasets in smoothed file\n')

        keys = list(raw_db.keys())
       
        # don't want to copy raw signal. will write 
        smooth_db.create_dataset('smoothed_signal',shape=raw_db['signal'].shape)
        keys.pop(keys.index('signal'))

        for key in keys:
            smooth_db.create_dataset(key,data=raw_db[key])

        smooth_db.create_dataset('smoothing_fwhm',data=self.smoothing_fwhm)

        _t.stop()

    # ----------------------------------------------------------------------------------------------

    def _gaussian_smooth(self,x,y,fwhm):

        """
        smooth data using gaussian convolution. x is a 1d array, y is a 1d array or 2d array 
        with same shape as x along 2nd axis. calls another method to actually do the smoothing
        """

        # gaussian stddev
        sigma = fwhm/np.sqrt(8*np.log(2))
        
        # step size along x axis
        dx = x[1]-x[0]

        # x axis for gaussian function
        gx = np.arange(-6*sigma,6*sigma,dx)
        gaussian = np.exp(-0.5*(gx/sigma)**2)

        if y.ndim == 1:
            smooth = self._gaussian_smooth_1d(x,y,gaussian)
        else:
            smooth = np.zeros(y.shape)
            for ii in range(smooth.shape[0]):
                smooth[ii,:] = self._gaussian_smooth_1d(x,y[ii,:],gaussian)

        return smooth

    # ----------------------------------------------------------------------------------------------

    def _gaussian_smooth_1d(self,x,y,gaussian):

        """
        smooth data using gaussian convolution. if there are nans in input array, interpolate to 
        remove the nans. 
        """

        # mask infinites with nans
        y = np.nan_to_num(y,nan=np.nan,posinf=np.nan,neginf=np.nan)

        nans = np.isnan(y)
        mask_nans = np.any(nans)

        if mask_nans:
            y_interpolated = np.interp(x,x[~nans],y[~nans])
        else:
            y_interpolated = np.copy(y)

        # handle 1d or 2d arrays separately
        smooth = convolve(y_interpolated, gaussian, mode='same', method='auto') 

        # normalize to same norm as input data
        norm_y = np.sqrt(np.sum(y_interpolated**2)) 
        norm_smooth = np.sqrt(np.sum(smooth**2))

        if norm_smooth < 1e-6:
            smooth[:] = 0.0
        else:
            smooth *= norm_y/norm_smooth

        # replace indices with nans with nans
        smooth[nans] = np.nan

        return smooth

    # ----------------------------------------------------------------------------------------------

    def calculate_background(self,raw_file,smoothed_file,num_bg_cuts=10,num_max_tries=100,
            delta_Q_length=0.5,delta_polar_angle=10,delta_azimuthal_angle=10,num_Q_point_procs=1):

        """
        do a 'rocking scan' around each Q-pt in the raw file and get random neighboring points 
        that lie in a piece of spherical shell the same Q-pt. 
        
        num_bg_cuts determines the number of neighboring points to use. if a neighboring point
        is all NaNs, the code tries again up to maximum num_max_tries attempts. 

        the neighboring points are bounded by |Q'| in |Q| +- delta_Q_length (in 1/A) and
        angle' in angle +- delta_angle (in degrees)

        the neighboring points are read from the smoothed file and the background is taken as the
        point-by-point minimum of all of the Q-points (neighboring and original).

        the idea is that the background is spherical, while the single-crystal scattering is 
        not. so looking at data lying on the same sphere but offset by a small angle means 
        we will move away from the phonons and only look at the background. ofcourse it is 
        possible that we might find other phonons at another Q-pt... so we use multiple 
        background cuts that are randomly chosen as a fail safe.

        NOTE: will parallelize the loop over Qpts to find background from neighbors. 

        """

        check_file(raw_file)
        check_file(smoothed_file)

        with h5py.File(raw_file,'r') as raw_db:
            self.num_Q = raw_db['H_rlu'].size
            self.energy = raw_db['DeltaE'][...]

        # split Qpt indices onto multiple process to get and calculate background 
        self.Q_indices_on_procs = np.array_split(np.arange(self.num_Q),num_Q_point_procs)

        self.delta_polar_angle = delta_polar_angle*np.pi/180
        self.delta_azimuthal_angle = delta_azimuthal_angle*np.pi/180
        self.delta_Q_length = delta_Q_length

        proc = 0
        self._get_background_by_rocking_scans(proc,raw_file,smoothed_file,
                num_bg_cuts,num_max_tries)

    # ----------------------------------------------------------------------------------------------

    def _get_background_by_rocking_scans(self,proc,raw_file,smoothed_file,
                num_bg_cuts,num_max_tries):

        """
        do a 'rocking scan' around each Q-pt in the raw file and get random neighboring points
        that lie in a piece of spherical shell the same Q-pt.
        """

        if proc == 0:
            _t = c_timer(f'get_background_by_rocking_scans')

        Q_indices = self.Q_indices_on_procs[proc]
        num_Q_on_proc = Q_indices.size

        if proc == 0:
            print(f'\n{num_Q_on_proc} Q-points on proc 0\n')

        with h5py.File(raw_file,mode='r',libver='latest',swmr=True) as raw_db, \
        h5py.File(smoothed_file,mode='r',libver='latest',swmr=True) as smooth_db:

            self.Q_len = raw_db['Q_len'][...]
            self.polar_angle = raw_db['polar_angle'][...]
            self.azimuthal_angle = raw_db['azimuthal_angle'][...]

            for ii, Q_index in enumerate(Q_indices):

                if proc == 0:
                    if ii%1000 == 0:
                        print(f'  Qpt {ii} out of {num_Q_on_proc}')

                Q_ii = self.Q_len[Q_index]
                polar_ii = self.polar_angle[Q_index]
                azimuthal_ii = self.azimuthal_angle[Q_index]

                self._get_neighboring_Q_points_in_shell(Q_ii,polar_ii,azimuthal_ii)

        if proc == 0:
            _t.stop()

    # ----------------------------------------------------------------------------------------------

    def _get_neighboring_Q_points_in_shell(self,Q,polar,azimuthal):

        """
        find neighboring Q-pts that are bounded by |Q|+-dQ, polar+-d_polar, azi+-d_azi
        """

        dQ = self.delta_Q_length
        dpolar = self.delta_polar_angle
        dazi = self.delta_azimuthal_angle

        # Q_len bounds
        len_inds = np.flatnonzero(self.Q_len >= Q-dQ)
        len_inds = np.intersect1d(np.flatnonzero(self.Q_len <= Q+dQ),len_inds)
            
        # polar angle bounds
        polar_inds = np.flatnonzero(self.polar_angle >= polar-dpolar)
        polar_inds = np.intersect1d(np.flatnonzero(self.polar_angle <= polar+dpolar),polar_inds)

        # polar angle bounds
        azimuthal_inds = np.flatnonzero(self.azimuthal_angle >= azimuthal-dazi)
        azimuthal_inds = np.intersect1d(
            np.flatnonzero(self.azimuthal_angle <= azimuthal+dazi),azimuthal_inds)

    # ----------------------------------------------------------------------------------------------






# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    raw_file = 'LSNO25_300K_parallel_test.hdf5'
    smoothed_file = 'LSNO25_300K_parallel_teste_SMOOTH.hdf5'

    bg_tools = c_background_tools()

#    bg_tools.make_smoothed_file(raw_file,smoothed_file,smoothing_fwhm=1.5,num_blocks=10)

    bg_tools.calculate_background(raw_file,smoothed_file,num_Q_point_procs=16)




