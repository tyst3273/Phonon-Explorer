
import h5py 
import numpy as np


"""

File added by Tyler Sterling, Apr. 2023

Description: data access class to interface with custom *.hdf5 files produced with 'MDE_tools.py' 
    (see Phonon-Explorer/Python-Code/custom_file_tools/MDE_tools.py). 

Explanation: RawData is the class that interfaces with phonon-explorer. These interfaces are
    the same for all file types. RawData takes the TextFile.Parameters class as an argument
    and gets the file path from it. the RawData.GetSlice method cuts data from the file. The 
    input options are the same as for the other data access classes. 

    RawData opens the file and gets the Q-points and energy arrays from it. It does all of 
    this inside of a 'plugin' called access_data_in_hdf5. The Q-points from the file are 
    compared to what the user requests. The nearest (or exact) Q-point is found by calculating the 
    Euclidean distance (in RLU) between the requested Q-point and the ones in the file. 
    The index of the nearest Q-point is used to slice the intensity and error from the file. 

    The energy, intensity, error, and Q-point found are all attached to RawData class as 
    attributes and are fetched upstream by phonon explorer. NOTE: if the user requests a 
    smaller energy range than what is in the file, the data are down sampled onto that energy
    grid. 

Whats new: The binning args and actual Q-point found in the file are attached as attributes
    to this class. The binning args are self.Delta* and self.e_step and the Q-point is
    self.Qpoint_rlu. These are here since they are specific to the data that is cut, so 
    are attached to the object holding the data.

"""



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

        # binning arguments in the hdf5 file
        self.Deltah = self._data_access_class.dh
        self.Deltak = self._data_access_class.dk
        self.Deltal = self._data_access_class.dl
        self.e_step = self._data_access_class.dE
        self.Projection_u = self._data_access_class.u
        self.Projection_v = self._data_access_class.v

    # ----------------------------------------------------------------------------------------------

    def GetSlice(self, bin_h, bin_k, bin_l, bin_e, Projection_u, Projection_v):
        """
        get signal, error from file at requested Qpt
        projections and binning size are ignored. bin_* args are used to determine the bin
        center and then the nearest one in the file is returned.
        """
        
        # requested bin center 
        Q = [np.array(bin_h).mean(),np.array(bin_k).mean(),np.array(bin_l).mean()]
        E_min = bin_e[0]; E_max = bin_e[2]
        
        # print the Q-point and E-range we are trying to get. 
        print("\nQ_h: {:01.3f}".format(Q[0]))
        print("Q_k: {:01.3f}".format(Q[1]))
        print("Q_l: {:01.3f}".format(Q[2]))
        print("bin_e: [{:01.3f}, {:01.3f}]".format(E_min, E_max),'\n')

        # try to get the data
        try:
            self.Energy, self.Intensity, self.Error, self.Qpoint_rlu = \
                self._data_access_class.get_intensity_and_error(Q)
        except Exception as e:
            print("No slice!")
            return 1

        self.Intensity *= 1000
        self.Error *= 1000

        # if requested grid is smaller than whats in file, downsample the data onto that grid
        _E_inds = np.intersect1d(np.flatnonzero(self.Energy >= bin_e[0]),
                    np.flatnonzero(self.Energy <= bin_e[2]))
        self.Energy = self.Energy[_E_inds]
        self.Intensity = self.Intensity[_E_inds]
        self.Error = self.Error[_E_inds]

        return 0

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
        get the Q_points and DeltaE array from hdf5 file. only want to have to do this once.
        also gets binning and projections args which are attached to RawData as attributes 
        for upstream code to deal with.
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

        print('Binning args in hdf5 file:')
        print('Deltah:',self.dh)
        print('Deltak:',self.dk)
        print('Deltal:',self.dl)
        print('e_step:',self.dE)
        print('Projection_u:',self.u)
        print('Projection_v:',self.v)
        print('\n-------------------------------------------\n')

        # allocate these now and just overwrite later; 
        # saves time since we have to do it repeatedly

        # vector from Qpts in file to Qpt requested by user (in rlu)
        self.Q_file_to_Q_user_vector = np.zeros((self.num_Q_in_file,3))

        # distance from Qpts in file to Qpt requested by user (in rlu)
        self.Q_file_to_Q_user_distance = np.zeros(self.num_Q_in_file)

    # ----------------------------------------------------------------------------------------------

    def get_intensity_and_error(self,Q,cutoff=0.1):
        """
        take Q in rlu, find nearest Qpt in file, and return signal and error arrays. 

        NOTE: if Qpt in file is sufficiently far from any Q-point in file, this not should not 
            return any data. raises an Exception which is handles upstream by phonon explorer.

            'cutoff' controls the minimum distance to a point in the file before
            an exception is raised. 
        """
        
        # get distance from user Q to all Qpts in file
        _Qpts = self.Q_points_rlu
        _Qvec = self.Q_file_to_Q_user_vector
        _Qdist = self.Q_file_to_Q_user_distance
        _Qvec[:,0] = _Qpts[:,0]-Q[0]
        _Qvec[:,1] = _Qpts[:,1]-Q[1]
        _Qvec[:,2] = _Qpts[:,2]-Q[2]
        _Qdist[...] = np.round(np.sqrt(np.sum(_Qvec**2,axis=1)),4)

        # find the closest Q-point
        Q_index = np.argsort(_Qdist)[0]

        # raise Exception if empty slice 
        # (i.e. no data at Q-pt within distance cutoff of requested Q)
        if _Qdist[Q_index] >= cutoff:
            raise Exception

        Qpoint_from_file = _Qpts[Q_index]

        # now get and return the data
        intensity, error = self._get_cut_from_hdf5(Q_index)

        return self.E, intensity, error, Qpoint_from_file

    # ----------------------------------------------------------------------------------------------

    def _get_cut_from_hdf5(self,Q_index):
        """
        attempt to get the data from hdf5 file. if fails, return nans
        """
        try:
            with h5py.File(self.file_name,'r') as db:
                intensity = db['signal'][Q_index,...]
                error = db['error'][Q_index,...]
        except Exception as ex:
            intensity = np.full(self.num_E_in_file,np.nan)
            error = np.full(self.num_E_in_file,np.nan)

        return intensity, error

    # ----------------------------------------------------------------------------------------------





