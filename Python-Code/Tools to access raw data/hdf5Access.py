
import h5py 
import numpy as np


# --------------------------------------------------------------------------------------------------

class RawData:

    def __init__(self,params):
        """
        interface to Dmitry's phonon explorer code
        """
        
        self.params = params
        file_name = params.sqw_path

        # this class wraps another one that actually interacts with the file
        self._data_access_class = access_data_in_hdf5(file_name)
        self._get_meta_data_from_hdf5_file()

    # ----------------------------------------------------------------------------------------------

    def _get_meta_data_from_hdf5_file(self):
        """
        check binning vs what is in file
        """

        self.dh = self._data_access_class.dh
        self.dk = self._data_access_class.dk
        self.dl = self._data_access_class.dl
        self.dE = self._data_access_class.dE
        self.u = self._data_access_class.u
        self.v = self._data_access_class.v
        self.w = self._data_access_class.w

    # ----------------------------------------------------------------------------------------------

    def GetSlice(self, bin_h, bin_k, bin_l, bin_e, Projection_u, Projection_v):
        """
        get signal, error from file at requested Qpt
        projections and binning size are ignored. bin_* args are used to determine the bin
        center and then the nearest one is returned.
        """
        
        # requested bin center 
        Q = [np.array(bin_h).mean(),np.array(bin_k).mean(),np.array(bin_l).mean()]

        print("\nQ_h: {:01.3f}".format(Q[0]))
        print("Q_k: {:01.3f}".format(Q[1]))
        print("Q_l: {:01.3f}".format(Q[2]))
        print("bin_e: [{:01.3f}, {:01.3f}, {:01.3f}]".format(bin_e[0], bin_e[1], bin_e[2]),'\n')

        # get the data
        self.Energy, self.Intensity, self.Error = self._data_access_class.get_signal_and_error(Q)
        self.Intensity *= 1000
        self.Error *= 1000

        # if requested grid is smaller than whats in file, downsample the data onto that grid
        _E_inds = np.intersect1d(np.flatnonzero(self.Energy >= bin_e[0]),
                    np.flatnonzero(self.Energy <= bin_e[2]))
        self.Energy = self.Energy[_E_inds]
        self.Intensity = self.Intensity[_E_inds]
        self.Error = self.Error[_E_inds]

        return 1

    # ----------------------------------------------------------------------------------------------



# --------------------------------------------------------------------------------------------------

class access_data_in_hdf5:

    def __init__(self,file_name):
        """
        class that gets data from hdf5 file
        """

        self.file_name = file_name
        self._load_Q_and_E()

    # ----------------------------------------------------------------------------------------------

    def _load_Q_and_E(self):
        """
        get the Q_points and DeltaE array from hdf5 file. only want to have to do this once
        """

        with h5py.File(self.file_name,'r') as db:
            self.num_Q_in_file = db['H_rlu'].size
            self.Q_points_rlu = np.zeros((self.num_Q_in_file,3))
            self.Q_points_rlu[:,0] = db['H_rlu'][...]
            self.Q_points_rlu[:,1] = db['K_rlu'][...]
            self.Q_points_rlu[:,2] = db['L_rlu'][...]
            self.E = db['DeltaE'][...]
            self.dh = db['H_bin_args'][...][1]/2.0
            self.dk = db['K_bin_args'][...][1]/2.0
            self.dl = db['L_bin_args'][...][1]/2.0
            self.dE = db['E_bin_args'][...][1]
            self.u = db['u'][...]
            self.v = db['v'][...]
            self.w = db['w'][...]

        # allocate these now and just reassign later; saves time since we have to do it repeatedly

        # vector from Qpts in file to Qpt requested by user (in rlu)
        self.Q_file_to_Q_user_vector = np.zeros((self.num_Q_in_file,3))

        # distance from Qpts in file to Qpt requested by user (in rlu)
        self.Q_file_to_Q_user_distance = np.zeros(self.num_Q_in_file)

    # ----------------------------------------------------------------------------------------------

    def get_signal_and_error(self,Q,cutoff=0.5):
        """
        take Q in rlu, find nearest Qpt in file, and return signal and error arrays. NOTE:
        if Qpt in file is sufficiently far from what is requested, should not return any data
        """

        _prec = 4; _eps = 0.005
        
        # get distance from user Q to all Qpts in file
        _Qpts = self.Q_points_rlu
        _Qvec = self.Q_file_to_Q_user_vector
        _Qdist = self.Q_file_to_Q_user_distance
        _Qvec[:,0] = _Qpts[:,0]-Q[0]
        _Qvec[:,1] = _Qpts[:,1]-Q[1]
        _Qvec[:,2] = _Qvec[:,2]-Q[2]
        _Qdist[...] = np.round(np.sqrt(np.sum(_Qvec**2,axis=1)),_prec)

        # find the closest Q-point
        Q_ind = np.argsort(_Qdist)[0]

        # raise Exception if empty slice (i.e. no data at Q-pt within distance cutoff of requested Q)
        if _Qdist[Q_ind] >= cutoff:
            print('*** NOTE ***\nno slice at requested Q!\n')
            raise Exception

        # print a warning if this Q-point isnt exaclty the same as requested by user
        _Q = _Qpts[Q_ind]
        for ii in range(3):
            if np.abs(_Q[ii]-Q[ii]) >= _eps:
                msg = '\n*** WARNING ***\n'
                msg += 'nearest Q-point is not exactly what is requested!\n'
                msg += f'nearest Q = ({_Q[0]: .3f},{_Q[1]: .3f},{_Q[2]: .3f})\n'
                msg += 'if this is not what is expected, pick a Q-point commensurate with the\n' \
                       'data in the file and run again.\n'
                print(msg)

        # now get and return the data
        sig, err = self._get_cut_from_hdf5(Q_ind)

        return self.E, sig, err

    # ----------------------------------------------------------------------------------------------

    def _get_cut_from_hdf5(self,Q_ind):
        """
        attempt to get the data from hdf5 file. if fails, return nans
        """
        try:
            with h5py.File(self.file_name,'r') as db:
                sig = db['signal'][Q_ind,...]
                err = db['error'][Q_ind,...]
        except Exception as ex:
            msg = '\n*** WARNING ***\n'
            msg += f'couldnt get signal/error from hdf5 file \'{self.file_name}\'.\n' \
                  'see the exception below for what went wrong.\n\n' \
                  'continuing, but this is bad! I hope you know what you are doing ...\n' 
            msg += '\nException:\n'+str(ex)+'\n\n'
            print(msg)
            sig = np.full(self.num_E_in_file,np.nan); err = np.full(self.num_E_in_file,np.nan)

        return sig, err

    # ----------------------------------------------------------------------------------------------






# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    import matplotlib.pyplot as plt

    Q = [8,-2,0]

    lw = 1
    ms = 7

    hdf5_file_name = 'LSNO25_300K_parallel.hdf5'
    access_300K = access_data_in_hdf5(hdf5_file_name)
    E, sig, err = access_300K.get_signal_and_error(Q)
    sig[np.flatnonzero(sig == 0.0)] = np.nan
    sig = sig/(1+1/(np.exp(E/(0.08617*300))-1))
    #plt.errorbar(E,sig,yerr=err,barsabove=True,ls='-',lw=lw,marker='o',ms=ms,c='b',label='300K')
    plt.plot(E,sig,ls='-',lw=lw,marker='o',ms=ms,c='b',label='300K')

    hdf5_file_name = 'LSNO25_220K_parallel.hdf5'
    access_220K = access_data_in_hdf5(hdf5_file_name)
    E, sig, err = access_220K.get_signal_and_error(Q)
    sig[np.flatnonzero(sig == 0.0)] = np.nan
    sig = sig/(1+1/(np.exp(E/(0.08617*220))-1))
    #plt.errorbar(E,sig,yerr=err,barsabove=True,ls='-',lw=lw,marker='o',ms=ms,c='r',label='220K')
    plt.plot(E,sig,ls='-',lw=lw,marker='o',ms=ms,c='r',label='220K')

    hdf5_file_name = 'LSNO25_5K_parallel.hdf5'
    access_5K = access_data_in_hdf5(hdf5_file_name)
    E, sig, err = access_5K.get_signal_and_error(Q)
    sig[np.flatnonzero(sig == 0.0)] = np.nan
    sig = sig/(1+1/(np.exp(E/(0.08617*5))-1))
    #plt.errorbar(E,sig,yerr=err,barsabove=True,ls='-',lw=lw,marker='o',ms=ms,c='m',label='5K')
    plt.plot(E,sig,ls='-',lw=lw,marker='o',ms=ms,c='m',label='5K')

    plt.legend()
    plt.show()

