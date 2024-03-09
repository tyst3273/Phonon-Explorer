
import h5py 
import numpy as np
import time

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

    def GetMultiDimData(self, H_min, H_max, K_min, K_max, L_min, L_max, E_min, E_max, Projection_u=0, Projection_v=0, hbin=0, kbin=0, lbin=0, ebin=0): #last 6 params are for compatibility of the API with other file formats. Not used here.
    
        if H_max-H_min<2*self.Deltah:
            mid=(H_max+H_min)/2
            H_min=mid-self.Deltah
            H_max=mid+self.Deltah
        if K_max-K_min<2*self.Deltak:
            mid=(K_max+H_min)/2
            K_min=mid-self.Deltak
            K_max=mid+self.Deltak
        if L_max-L_min<2*self.Deltak:
            mid=(L_max+L_min)/2
            L_min=mid-self.Deltal
            L_max=mid+self.Deltal
        if E_max-E_min<self.e_step:
            mid=(E_max+E_min)/2
            E_min=mid-self.e_step/2
            E_max=mid+self.e_step/2

        min=[H_min,K_min,L_min]
        max=[H_max,K_max,L_max]

        try:
            self.Energy, self.Intensity, self.Error, self.Qpoint_rlu, code = \
                self._data_access_class.get_intensity_and_error_in_Q_range(min_vals, max_vals)
        except Exception as e:
            print("No slice!")
            return -1

        mask = np.logical_and(self.Energy >= E_min, self.Energy <= E_max)
        _E_inds = np.where(mask)[0]
        
        self.Energy = self.Energy[_E_inds]
        self.Intensity = self.Intensity[_E_inds]
        self.Error = self.Error[_E_inds]

        return 0
        
    def GetRockingScan(self, Phi_min, Phi_max, Theta_min, Theta_max):
        return

    def GetCollectionOfQs(self,Q):
#        if (1==1):
        try:
            self.Energy, intensityArray, errorArray, Qpoints, code = \
            self._data_access_class.get_intensity_and_error_in_Array_of_Qs(Q)
            return self.Energy, intensityArray, errorArray, Qpoints, code
        except Exception as e:
            print("No slice!")
            return -1

        return 0


    def GetSlice(self, bin_h, bin_k, bin_l, bin_e, Projection_u, Projection_v):
        """
        get signal, error from file at requested Qpt
        projections and binning size are ignored. bin_* args are used to determine the bin
        center and then the nearest one in the file is returned.
        """
        
        # requested bin center 
        Q = [np.array(bin_h).mean(),np.array(bin_k).mean(),np.array(bin_l).mean()]
        E_min = bin_e[0]; E_max = bin_e[2]

#        if (1==1):
        try:
            self.Energy, self.Intensity, self.Error, self.Qpoint_rlu, code = \
                self._data_access_class.get_intensity_and_error_Single_Q(Q)
        except Exception as e:
            print("No slice!")
            return -1
        self.Intensity *= 1000
        self.Error *= 1000

        # if requested grid is smaller than whats in file, downsample the data onto that grid
#        _E_inds = np.intersect1d(np.flatnonzero(self.Energy >= bin_e[0]), np.flatnonzero(self.Energy <= bin_e[2]))
        mask = np.logical_and(self.Energy >= bin_e[0], self.Energy <= bin_e[2])
        _E_inds = np.where(mask)[0]
        
        self.Energy = self.Energy[_E_inds]
        self.Intensity = self.Intensity[_E_inds]
        self.Error = self.Error[_E_inds]

        return code

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
#            self.Q_len = db['Q_len'][...]
#            self.phi=db['polar_angle'][...]
#            self.psi=db['azimuthal_angle'][...]
#            for i in range(0,10):
#                print("{:.2f}, {:.2f}, {:.2f}".format(self.Q_points_rlu[i,0],self.Q_points_rlu[i,1],self.Q_points_rlu[i,2]))
#            time.sleep(20)
            self.energy = db['DeltaE'][...]
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

    def filter_array(self, arr, min_vals, max_vals):
    
        mask = np.all((arr >= min_vals) & (arr <= max_vals), axis=1)
        filtered_indices = np.where(mask)[0]
        filtered_arr = arr[mask]
        
        return filtered_arr, filtered_indices
        
    def get_intensity_and_error_in_Q_range(self, min_vals, max_vals):

        Q_Points=[]
        Q_Point_Indices=[]
        
        for i in range(0,len(min_vals)):
            Q_Point, Q_Point_Index=self.filter_array(self.Q_points_rlu, min_vals[i], max_vals[i])
            try:
                for i in range(0,len(Q_Point)):
                    Q_Points.append(Q_Point[i])
                    Q_Point_Indices.append(Q_Point_Index[i])
            except:
                pass
        code=0

        Num_Q_Points = len(Q_Point_Indices)
        sorted_indices = np.argsort(Q_Point_Indices)
        Q_Points_sorted = []
        sorted_Q_Point_Indices = []
        for i in range (0,len(sorted_indices)):
            positions = np.where(sorted_Q_Point_Indices == Q_Point_Indices[sorted_indices[i]])[0]
            if len(positions)==0:
                sorted_Q_Point_Indices.append(Q_Point_Indices[sorted_indices[i]])
                Q_Points_sorted.append(Q_Points[sorted_indices[i]])
        
        intensity = np.full(Num_Q_Points,np.NaN)
        error = np.full(Num_Q_Points,np.NaN)

        intensityArray, errorArray = self._get_cut_from_hdf5(sorted_Q_Point_Indices)
        
        return Q_Points_sorted, intensityArray, errorArray, code


    def get_intensity_and_error_Single_Q(self,Q,cutoff=0.5):
            
        min=[[Q[0]-self.dh,Q[1]-self.dk,Q[2]-self.dl]]
        max=[[Q[0]+self.dh,Q[1]+self.dk,Q[2]+self.dl]]
        
        Qpoints, intensity, error, code = self.get_intensity_and_error_in_Q_range(min, max)
        Qpoint_from_file=Qpoints[0]

        return self.energy, intensity[0], error[0], Qpoint_from_file, code

    def get_intensity_and_error_in_Array_of_Qs(self,Q):
        
        min=Q-[self.dh,self.dk,self.dl]
        max=Q+[self.dh,self.dk,self.dl]


        Qpoints, intensity, error, code = self.get_intensity_and_error_in_Q_range(min, max)

        return self.energy, intensity, error, Qpoints, code




