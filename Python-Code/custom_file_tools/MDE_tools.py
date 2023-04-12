
"""
Author: Tyler C. Sterling
Email: ty.sterling@colorado.edu
Affil: University of Colorado Boulder, Raman Spectroscopy and Neutron Scattering Lab
Date: 01/24/2022
Description:
    tools to:
        - programatically get data from mantid MDE files (in nexus format) using mantid
            and write to custom hdf5 for quick access from phonon explorer
        - calculate and subtract background from raw data in custom hdf5 file
"""

from timeit import default_timer
import numpy as np
import h5py 
import os
from scipy.signal import convolve
import multiprocessing as mp

# --------------------------------------------------------------------------------------------------
# auxilliary utilities for classes below
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


# --------------------------------------------------------------------------------------------------
# class to interact with mantid
# --------------------------------------------------------------------------------------------------

class c_MDE_tools:

    # ----------------------------------------------------------------------------------------------

    def __init__(self,MDE_file_name=None):
        """
        class to handle everything to reduce MDE => binned histo => textfile
        the point is using mantid is slow as hell with phonon expolorer so we want to 
        convert MDE to a 'database' that can be parsed quickly with SQL...
        """

        self.import_mantid()
        
        if MDE_file_name.endswith('.nxs'):
            self.MDE_ws_name = MDE_file_name[:-len('.nxs')]
            self.MDE_ws_name = os.path.basename(self.MDE_ws_name)
        else:
            _type = MDE_file_name.split('.')[-1].strip()
            msg = f'unknown filetype \'{_type}\'. use a *.nxs file!\n'
            crash(msg) 
        
        self.histo_ws_name = 'histo_ws'

        self.load_MDE(MDE_file_name)
        self.print_dimensions(self.MDE_ws_name)
        self.get_lattice(self.MDE_ws_name)

    # ----------------------------------------------------------------------------------------------

    def import_mantid(self):
        """
        try to import mantid. put in this method since mantid is a PITA and want control 
        over debugging errors
        """
        try:
            import mantid.simpleapi as msi
        except Exception as _ex:
            err_msg = 'couldnt import mantid! see exception below...'
            crash(err_msg,_ex)
        self.msi = msi

    # ----------------------------------------------------------------------------------------------

    def load_MDE(self,file_name):
        """
        load the MDE workspace
        """
        _t = c_timer('load_MDE',units='m')
        check_file(file_name)
        self.load_ws(file_name=file_name,ws_name=self.MDE_ws_name)
        _t.stop()

    # ----------------------------------------------------------------------------------------------

    def check_ws(self,ws_name,crash=True,bonus_msg=None):
        """
        check if a workspace exists in memory on mantid backend
        if crash=True, it crashes if ws not loaded. 
        if crash=False, it just returns (bool) wheter or not its loaded 
        """
        loaded = self.msi.mtd.doesExist(ws_name)
        if not crash:
            return loaded
        if not loaded:
            msg = f'the workspace \'{ws_name}\' isnt loaded!\n'
            if bonus_msg is not None:
                msg += bonus_msg+'\n'
            crash(msg)

    # ----------------------------------------------------------------------------------------------

    def get_ws(self,ws_name):
        """
        get object handle for workspace that exists in mantid backend
        """
        self.check_ws(ws_name)
        return self.msi.mtd[ws_name]

    # ----------------------------------------------------------------------------------------------
    
    def load_ws(self,file_name,ws_name,force_load=False):
        """
        check if the MD file is already loaded. if not, load it
        """
        load = not(self.check_ws(ws_name,crash=False))
        if force_load or load:
            msg = f'\nloading workspace \'{ws_name}\' from file \'{file_name}\'\n'
            print(msg)
            self.msi.Load(Filename=file_name,OutputWorkspace=ws_name)       
        else:
            msg = f'\nworkspace \'{ws_name}\' is already loaded! continuing...\n'
            print(msg)

    # ----------------------------------------------------------------------------------------------

    def print_dimensions(self,ws_name=None):
        """
        print the dimensions in the workspace
        """
        if ws_name is None:
            ws_name = self.MDE_ws_name

        ws = self.get_ws(ws_name)
        ndims = ws.getNumDims()

        dim_names = np.zeros(ndims,dtype=object)
        dim_mins = np.zeros(ndims)
        dim_maxs = np.zeros(ndims)
        dim_units = np.zeros(ndims,dtype=object)
        dim_num_bins = np.zeros(ndims,dtype=int)

        msg = f'dimensions in workspace \'{ws_name}\':\n'
        msg += f'num_dimensions: {ndims}\n'
        msg += '--------------- (name) -------- (units) ------- (min) ---- (max) -- (bins)'
        for ii in range(ndims):
            _dim = ws.getDimension(ii)
            _min = _dim.getMinimum()
            _max = _dim.getMaximum()
            _num_bins = _dim.getNBins()
            _name = _dim.name
            _units = _dim.getUnits()

            dim_names[ii] = _name
            dim_mins[ii] = _min
            dim_maxs[ii] = _max
            dim_units[ii] = _units
            dim_num_bins[ii] = _num_bins

            msg += f'\ndim[{ii}]: {_name:>16s} {_units:>16s}   ' \
                   f'{_min: >10.4f} {_max: >10.4f} {_num_bins:5g}'
        print(msg+'\n')

    # ----------------------------------------------------------------------------------------------

    def get_lattice(self,ws_name):
        """
        get the lattice vectors for the data in the MDE file
        """
        ws = self.get_ws(ws_name)

        _ = ws.getExperimentInfo(0)
        _ = _.sample()

        if not _.hasOrientedLattice():
            msg = f'couldnt find OrientedLattice object in workspace \'{ws_name}\'\n'
            crash(msg)

        _uc = _.getOrientedLattice()

        # lattice vectors lenghts in Angstrom
        a = _uc.a1(); b = _uc.a2(); c = _uc.a3()
        self.a = a; self.b = b; self.c = c

        # unitcell angles in Radians
        alpha = _uc.alpha1(); beta = _uc.alpha2(); gamma = _uc.alpha3()
        self.alpha = alpha; self.beta = beta; self.gamma = gamma

        self.get_lattice_vectors_from_params()
        self.get_reciprocal_lattice_vectors()

        # print lattice info to screen
        alpha *= 180/np.pi; beta *= 180/np.pi; gamma *= 180/np.pi
        msg = f'lattice parameters from workspace \'{ws_name}\'\n'
        msg += '---- a ------- b ------- c --- (Angstrom)\n'
        msg += f'{a:9.5f} {b:9.5f} {c:9.5f}\n'
        msg += '-- alpha --- beta ---- gamma - (degrees)\n'
        msg += f'{alpha:9.5f} {beta:9.5f} {gamma:9.5f}\n'
        _x = self.lattice_vectors[0,:]
        _y = self.lattice_vectors[1,:]
        _z = self.lattice_vectors[2,:]
        msg += '\n lattice vectors (Angstrom)\n'
        msg += f'{_x[0]: 9.5f} {_x[1]: 9.5f} {_x[2]: 9.5f}\n'
        msg += f'{_y[0]: 9.5f} {_y[1]: 9.5f} {_y[2]: 9.5f}\n'
        msg += f'{_z[0]: 9.5f} {_z[1]: 9.5f} {_z[2]: 9.5f}\n'
        _x = self.recip_lattice_vectors[0,:]
        _y = self.recip_lattice_vectors[1,:]
        _z = self.recip_lattice_vectors[2,:]
        msg += '\n reciprocal lattice vectors (1/Angstrom)\n'
        msg += f'{_x[0]: 9.5f} {_x[1]: 9.5f} {_x[2]: 9.5f}\n'
        msg += f'{_y[0]: 9.5f} {_y[1]: 9.5f} {_y[2]: 9.5f}\n'
        msg += f'{_z[0]: 9.5f} {_z[1]: 9.5f} {_z[2]: 9.5f}\n'
        print(msg)

    # ----------------------------------------------------------------------------------------------

    def get_lattice_vectors_from_params(self):
        """
        take lattice params (cell lenghts/angles) and return lattice vector array
        note, lengths should be in Angstrom, angles in radians

        note, lattice_vectors are ROW vectors
        """
        lattice_vectors = np.zeros((3,3))
        lattice_vectors[0,:] = [self.a,0,0]
        lattice_vectors[1,:] = [self.b*np.cos(self.gamma),self.b*np.sin(self.gamma),0]
        lattice_vectors[2,:] = [self.c*np.cos(self.beta),
                                    -self.c*np.sin(self.beta)*np.cos(self.alpha),self.c]
        self.lattice_vectors = lattice_vectors

    # ----------------------------------------------------------------------------------------------

    def get_reciprocal_lattice_vectors(self):
        """
        get reciprocal lattice vectors from lattice vectors. we need them to go from HKL to 1/A
        note: the UnitCell class of mantid seems to offer methods for this (e.g. qFromHKL()) but
        I cant figure out how to use them. I hate mantid.
        """
        lattice_vectors = self.lattice_vectors
        _2pi = np.array([[1,0,0],[0,1,0],[0,0,1]],dtype=float)*2*np.pi
        self.recip_lattice_vectors = np.linalg.solve(lattice_vectors,_2pi).T

    # ----------------------------------------------------------------------------------------------

    def bin_MDE_chunks(self,H_bin_args,K_bin_args,L_bin_args,E_bin_args,num_chunks=[1,1,1],
        merged_file_name='merged_sparse_histo.hdf5',u=[1,0,0],v=[0,1,0],w=None):

        """
        split requested binning into chunks and bin over small chunks separately. merge the 
        results of them all into a single file
        """
        
        _t = c_timer('bin_MDE_chunks',units='m')

        self.merged_file_name = merged_file_name
        if os.path.exists(merged_file_name):
            msg = 'merged file already exists. removing it ...\n'
            print(msg)
            os.remove(merged_file_name)
        
        self.u = u; self.v = v; self.w = w
    
        # get bin edges; needed to for my algorithm to partion array into chunks
        self.H_bin_args = H_bin_args
        self.H_bin_edges, self.H_range, self.dH = self._get_bin_edges(H_bin_args)
        self.K_bin_args = K_bin_args
        self.K_bin_edges, self.K_range, self.dK = self._get_bin_edges(K_bin_args)
        self.L_bin_args = L_bin_args
        self.L_bin_edges, self.L_range, self.dL = self._get_bin_edges(L_bin_args)
        self.E_bin_args = E_bin_args
        self.E_bin_edges, self.E_range, self.dE = self._get_bin_edges(E_bin_args)

        H_bin_centers = np.round((self.H_bin_edges[:-1]+self.H_bin_edges[1:])/2,2)
        K_bin_centers = np.round((self.K_bin_edges[:-1]+self.K_bin_edges[1:])/2,2)
        L_bin_centers = np.round((self.L_bin_edges[:-1]+self.L_bin_edges[1:])/2,2)
        E_bin_centers = np.round((self.E_bin_edges[:-1]+self.E_bin_edges[1:])/2,2)
        print('H_bin_centers\n',H_bin_centers)
        print('K_bin_centers\n',K_bin_centers)
        print('L_bin_centers\n',L_bin_centers)
        print('E_bin_centers\n',E_bin_centers)

        # split bin edges (well actually bin list args) into chunks
        self.num_chunks = num_chunks        
        self._split_bins_on_chunks()

        # get the data for each voxel
        self._loop_over_chunks()

        _t.stop()

    # ----------------------------------------------------------------------------------------------

    def _loop_over_chunks(self):
        """
        loop over the bin edges and cut the data
        """

        msg = f'\nsplitting binning into {self.num_grid} voxels\n'
        print(msg)

        count = 0
        # loop over grid voxels
        for ii, H_bins in enumerate(self.H_chunk_bin_args):
            for ii, K_bins in enumerate(self.K_chunk_bin_args):
                for ii, L_bins in enumerate(self.L_chunk_bin_args):

                    _vt = c_timer(f'chunk[{count}]')

                    # bining per chunk 
                    msg = '-----------------------------------------------------------------'
                    msg += f'\nnow on grid voxel {count} spanning '
                    msg += f'\n H: {H_bins[0]:.3f} => {H_bins[2]:.3f},'
                    msg += f'  K: {K_bins[0]:.3f} => {K_bins[2]:.3f},'
                    msg += f'  L: {L_bins[0]:.3f} => {L_bins[2]:.3f}\n'
                    print(msg)

                    # get the data for this voxel
                    self.bin_MDE(H_bins,K_bins,L_bins,self.E_bin_args,
                            self.u,self.v,self.w)

                    # append all the non-empty bins to the file
                    self.append_sparse_to_hdf5(self.merged_file_name)

                    count += 1

                    _vt.stop()
                    print('')

    # ----------------------------------------------------------------------------------------------

    def get_dim_array(self,dim):
        """
        convert mantid dimension object to numpy array containing bin centers
        """
        _min = dim.getMinimum()
        _max = dim.getMaximum()
        _nbins = dim.getNBins()
        _d = (_max-_min)/_nbins
        _bins = np.arange(_min+_d/2,_max,_d)
        return _bins

    # ----------------------------------------------------------------------------------------------

    def get_dim_mesh(self,dim_array_list):
        """
        take list of arrays containing bin centers and convert to an Ndim coordinate grid
        """
        return np.meshgrid(*dim_array_list,indexing='ij')
    
    # ----------------------------------------------------------------------------------------------

    def get_cartesian_coords(self):
        """
        go from rlu to cartesian coords; vectorized method below is ~ 10x faster than looping,
        but neither are slow
        """

        _t = c_timer('get_cartesian_coords',units='s')

        _r_lat_vecs = self.recip_lattice_vectors
        _nQ = self.Q_mesh.shape[0]

        self.cartesian_Q_mesh = np.zeros((_nQ,3),dtype=float)
        
        for ii in range(3):
            _b = _r_lat_vecs[ii,:].reshape(1,3) # row vectors
            _b = np.tile(_b,reps=(_nQ,1))
            _Q = self.Q_mesh[:,ii].reshape(_nQ,1)
            _Q = np.tile(_Q,reps=(1,3))
            self.cartesian_Q_mesh += _Q*_b

        _t.stop()

    # ----------------------------------------------------------------------------------------------

    def get_polar_coords(self):
        """
        go from cartesian coords to polar. 0th column is magnitude of |Q|, 1st is polar angle, 2nd
        is azimuthal angle. the "z-axis" is the polar axis. 
        """

        _t = c_timer('get_polar_coords',units='s')

        _cart_Q = self.cartesian_Q_mesh
        _nQ = _cart_Q.shape[0]

        self.polar_Q_mesh = np.zeros((_nQ,3),dtype=float)

        # magnitude of Q
        self.polar_Q_mesh[:,0] = np.sqrt(np.sum(_cart_Q**2,axis=1))

        # mask where |Q| = 0; 
        # note, this will make the angle ~0 degrees but since |Q| ~ 0, it doesnt matter
        _mag_Q = np.copy(self.polar_Q_mesh[:,0])
        _mag_Q[np.flatnonzero(_mag_Q < 1e-6)] = 1.0

        # polar angle
        self.polar_Q_mesh[:,1] = np.arccos(_cart_Q[:,2]/_mag_Q)

        # azimuthal angle
        self.polar_Q_mesh[:,2] = np.arctan2(_cart_Q[:,1],_cart_Q[:,0])

        _t.stop()

    # ----------------------------------------------------------------------------------------------
    
    def append_sparse_to_hdf5(self,output_file_name):
        """
        write sparse data to hdf5 file; needs to be resizable... gonna be hard
        """

        _t = c_timer('append_sparse_to_hdf5')

        # get the sparse data 
        self.get_sparse_arrays_from_histo_ws()

        # get Qpts in cartesian coords
        self.get_cartesian_coords()

        # get Qpts in polar coords of data
        self.get_polar_coords()
    
        msg = f'\nwriting sparse data to hdf5 file:\n  \'{output_file_name}\'\n'
        print(msg)

        _exists = os.path.exists(output_file_name)
        with h5py.File(output_file_name,mode='a',libver='latest') as db:

            # create datasets if file doesnt exist
            if not _exists:
                self._create_datasets(db)

            # if file already exists, resize and add data to datasets
            else:
                self._append_datasets(db)

        _t.stop()

    # ----------------------------------------------------------------------------------------------

    def _create_datasets(self,db):
        """
        create datasets in file if it doesnt already exits
        """

        _nE = self.dim_arrays[3].size

        # see below... i used 'chunks=True' to let h5py estimate these sizes
        _chunks = (125,12) # chunk size nQxnE arrays (resizable along 1st dim)
        _flat_chunks = (1000,) # chunk size for nQ arrays (resizable)

        db.create_dataset('u',data=self.u)
        db.create_dataset('v',data=self.v)
        db.create_dataset('w',data=self.w)
        db.create_dataset('lattice_vectors',data=self.lattice_vectors)
        db.create_dataset('recip_lattice_vectors',data=self.recip_lattice_vectors)
        db.create_dataset('a',data=self.a)
        db.create_dataset('b',data=self.b)
        db.create_dataset('c',data=self.c)
        db.create_dataset('alpha',data=self.alpha)
        db.create_dataset('beta',data=self.beta)
        db.create_dataset('gamma',data=self.gamma)

        db.create_dataset('Dim_0_name',data=self.dim_names[0])
        db.create_dataset('Dim_1_name',data=self.dim_names[1])
        db.create_dataset('Dim_2_name',data=self.dim_names[2])
        db.create_dataset('Dim_3_name',data=self.dim_names[3])

        db.create_dataset('H_bin_args',data=np.array([self.H_range[0],self.dH,self.H_range[1]]))
        db.create_dataset('K_bin_args',data=np.array([self.K_range[0],self.dK,self.K_range[1]]))
        db.create_dataset('L_bin_args',data=np.array([self.L_range[0],self.dL,self.L_range[1]]))
        db.create_dataset('E_bin_args',data=np.array([self.E_range[0],self.dE,self.E_range[1]]))
        db.create_dataset('DeltaE',data=self.dim_arrays[3],dtype=np.float32)

        # resizable arrays
    
        # Q rlu
        db.create_dataset('H_rlu',data=self.Q_mesh[:,0],
                                maxshape=(None,),dtype=np.float32,chunks=_flat_chunks)
        db.create_dataset('K_rlu',data=self.Q_mesh[:,1],
                                maxshape=(None,),dtype=np.float32,chunks=_flat_chunks)
        db.create_dataset('L_rlu',data=self.Q_mesh[:,2],
                                maxshape=(None,),dtype=np.float32,chunks=_flat_chunks)

        # Q cartesian
        db.create_dataset('H_cartesian',data=self.cartesian_Q_mesh[:,0],
                                maxshape=(None,),dtype=np.float32,chunks=_flat_chunks)
        db.create_dataset('K_cartesian',data=self.cartesian_Q_mesh[:,1],
                                maxshape=(None,),dtype=np.float32,chunks=_flat_chunks)
        db.create_dataset('L_cartesian',data=self.cartesian_Q_mesh[:,2],
                                maxshape=(None,),dtype=np.float32,chunks=_flat_chunks)

        # Q polar
        db.create_dataset('Q_len',data=self.polar_Q_mesh[:,0],
                                maxshape=(None,),dtype=np.float32,chunks=_flat_chunks)
        db.create_dataset('polar_angle',data=self.polar_Q_mesh[:,1],
                                maxshape=(None,),dtype=np.float32,chunks=_flat_chunks)
        db.create_dataset('azimuthal_angle',data=self.polar_Q_mesh[:,2],
                                maxshape=(None,),dtype=np.float32,chunks=_flat_chunks)

        # the data
        db.create_dataset('signal',data=self.signal,
                                maxshape=(None,_nE),dtype=np.float64,chunks=_chunks)
        db.create_dataset('error',data=self.error,
                                maxshape=(None,_nE),dtype=np.float64,chunks=_chunks)
        #db.create_dataset('num_events',data=self.num_events,
        #                        maxshape=(None,_nE),dtype=int,chunks=_chunks)

    # ----------------------------------------------------------------------------------------------

    def _append_datasets(self,db):
        """
        append datasets to file if it already exists
        """

        _nE = self.dim_arrays[3].size
        _nQ = self.Q_mesh[:,0].size
        _nfile = db['signal'].shape[0]
        msg = f'resizing arrays from ({_nfile},{_nE}) to ({_nQ+_nfile},{_nE})\n'
        print(msg)

        # rlu
        db['H_rlu'].resize(_nQ+_nfile,axis=0)
        db['H_rlu'][_nfile:] = self.Q_mesh[:,0]
        db['K_rlu'].resize(_nQ+_nfile,axis=0)
        db['K_rlu'][_nfile:] = self.Q_mesh[:,1]
        db['L_rlu'].resize(_nQ+_nfile,axis=0)
        db['L_rlu'][_nfile:] = self.Q_mesh[:,2]

        # cartesian
        db['H_cartesian'].resize(_nQ+_nfile,axis=0)
        db['H_cartesian'][_nfile:] = self.cartesian_Q_mesh[:,0]
        db['K_cartesian'].resize(_nQ+_nfile,axis=0)
        db['K_cartesian'][_nfile:] = self.cartesian_Q_mesh[:,1]
        db['L_cartesian'].resize(_nQ+_nfile,axis=0)
        db['L_cartesian'][_nfile:] = self.cartesian_Q_mesh[:,2]
        
        # polar
        db['Q_len'].resize(_nQ+_nfile,axis=0)
        db['Q_len'][_nfile:] = self.polar_Q_mesh[:,0]
        db['polar_angle'].resize(_nQ+_nfile,axis=0)
        db['polar_angle'][_nfile:] = self.polar_Q_mesh[:,1]
        db['azimuthal_angle'].resize(_nQ+_nfile,axis=0)
        db['azimuthal_angle'][_nfile:] = self.polar_Q_mesh[:,2]

        # data
        db['signal'].resize(_nQ+_nfile,axis=0)
        db['signal'][_nfile:,:] = self.signal[...]
        db['error'].resize(_nQ+_nfile,axis=0)
        db['error'][_nfile:,:] = self.error[...]
        #db['num_events'].resize(_nQ+_nfile,axis=0)
        #db['num_events'][_nfile:,:] = self.num_events[...]

    # ----------------------------------------------------------------------------------------------

    def get_sparse_arrays_from_histo_ws(self,zero_nans=False):
        """
        remove all (Q-point) bins with 0 events
        """

        _t = c_timer('sparsify')

        msg = '\nsparsifying...'
        print(msg)

        self.get_arrays_from_histo_ws()

        _shape = self.signal.shape
        _nQ = np.prod(_shape[:-1])
        _nE = _shape[-1]

        self.signal = np.reshape(self.signal,(_nQ,_nE))
        self.error = np.reshape(self.error,(_nQ,_nE))
        self.num_events = np.reshape(self.num_events,(_nQ,_nE))

        # find Q-pts where whole array (along E) is nans
        _nans = np.isnan(self.signal) # where nans are
        _infs = np.isinf(self.signal) # where +/- infs are
        _inds = ~(~_nans * ~_infs) # True for Q,E that are either nan or inf
        _inds = ~np.all(_inds,axis=1) # True for Q-pts that are all nans or infs
    
        # replace infs/nans with 0s ?
        if zero_nans:
            self.signal[(_nans)] = 0.0
            self.signal[(_infs)] = 0.0
            self.error[(_nans)] = 0.0
            self.error[(_infs)] = 0.0
            self.num_events[(_nans)] = 0.0
            self.num_events[(_infs)] = 0.0

        # strip empties
        self.signal = self.signal[_inds,:]
        self.error = self.error[_inds,:]
        self.num_events = self.num_events[_inds,:]

        # flatten Q mesh after stripping empties        
        _Q_mesh = [_Q.flatten()[_inds] for _Q in self.Q_mesh]
        self.Q_mesh = np.zeros((_Q_mesh[0].size,3),dtype=float)
        self.Q_mesh[:,0] = _Q_mesh[0][...]
        self.Q_mesh[:,1] = _Q_mesh[1][...]
        self.Q_mesh[:,2] = _Q_mesh[2][...]

        # print info 
        msg = 'original shape: '+'x'.join([str(_) for _ in _shape])
        msg += f'\nnumber of Q-points: {_nQ:d}'
        msg += f'\nnumber of non-empty Q-points: {_inds.size:d}\n'
        print(msg)

        _t.stop()

    # ----------------------------------------------------------------------------------------------

    def get_arrays_from_histo_ws(self):
        """
        strip the data in the histo workspace in numpy arrays
        NOTE: see input arg NumEvNorm and effects in 'SaveMDToAscii' provided by A. Savici
        """

        self.check_ws(self.histo_ws_name)
        ws = self.get_ws(self.histo_ws_name)

        self.signal = ws.getSignalArray()
        self.error = np.sqrt(ws.getErrorSquaredArray())
        self.num_events = ws.getNumEventsArray()

        # get bin center arrays for all dimensions
        _dims = [ws.getDimension(ii) for ii in range(ws.getNumDims())]
        _dim_arrays = [self.get_dim_array(d) for d in _dims]
        _dim_mesh = self.get_dim_mesh(_dim_arrays)
        _dim_names = [d.name for d in _dims]
        _dim_num_bins = [d.getNBins() for d in _dims]
        _Q_mesh = self.get_dim_mesh(_dim_arrays[:-1])
        _Q_names = [d.name for d in _dims[:-1]]

        self.dim_arrays = _dim_arrays
        self.dim_mesh = _dim_mesh
        self.dim_names = _dim_names
        self.Q_names = _Q_names
        self.Q_mesh = _Q_mesh

    # ----------------------------------------------------------------------------------------------

    def _split_bins_on_chunks(self):
        """
        new method explicity enumerates bin edges and then grabs the 'chunks' of edges to loop
        over; old method of splitting range and then coming up with ranges didnt necessarily
        result in bins commensurate with the grid user expects. 
        """

        # number of chunks to split binning along each axis
        nH = self.num_chunks[0]; nK = self.num_chunks[1]; nL = self.num_chunks[2]
        self.num_grid = nH*nK*nL
        self.nH = nH; self.nK = nK; self.nL = nL

        self.H_chunk_bin_args = self._get_chunk_bin_args(self.H_bin_edges,self.dH,self.nH)
        self.K_chunk_bin_args = self._get_chunk_bin_args(self.K_bin_edges,self.dK,self.nK)
        self.L_chunk_bin_args = self._get_chunk_bin_args(self.L_bin_edges,self.dL,self.nL)

    # ----------------------------------------------------------------------------------------------

    def _get_chunk_bin_args(self,edges,d,n):
        """
        get binning range for each 'chunk' along axis
        """
        if n >= edges.size:
            msg = 'you are trying to split binning into more chunks than there are bins.\n' \
                  'pick a bigger range to bin or use fewer chunks.\n'
            crash(msg)
        _split = np.array_split(edges,n)
        chunk_bin_args = []
        for ii in range(n):             
            _s = _split[ii]
            # pad by d; upper bin is too big, but will never be empty (will crash if empty)
            chunk_bin_args.append([_s[0],d,_s[-1]+d])  
        return chunk_bin_args

    # ----------------------------------------------------------------------------------------------
    
    def _get_bin_edges(self,bins):
        """     
        return array of bin edges
        """
        if len(bins) != 3:
            msg = 'intergrating out dimensions is not supported.\n' \
                  'all binning args must be given as [lo, d, hi]\n'
            crash(msg)

        bin_range = [bins[0],bins[2]]
        d = bins[1]
        edges = np.arange(bin_range[0],bin_range[1]+d,d)

        return edges, bin_range, d

    # ----------------------------------------------------------------------------------------------

    def bin_MDE(self,H_bins=None,K_bins=None,L_bins=None,E_bins=None,u=[1,0,0],v=[0,1,0],w=None):
        """
        bin the events in the MDE workspace into a histogram workspace; note that this doesnt 
        really return anything or create new attributes. the produced data are stored in the 
        histogram workspace
        """

        MDE_ws = self.get_ws(self.MDE_ws_name)

        # get u, v vectors to determin axes along which to bin
        u = np.array(u); v = np.array(v)
        if w is None:
            w = np.cross(u,v)
        else:
            w = np.array(w)

        # copy projection to attribute the convert to str to pass to mantid
        self.u = u
        self.v = v
        self.w = w
        u = self.get_bin_str(u)
        v = self.get_bin_str(v)
        w = self.get_bin_str(w)

        # copy bin args as attributes then convert to str to pass to mantid
        H_bins = self.get_bin_str(H_bins)
        K_bins = self.get_bin_str(K_bins)
        L_bins = self.get_bin_str(L_bins)
        E_bins = self.get_bin_str(E_bins)

        msg = 'binning MDE workspace\n\n'
        msg += f'u: {u:>9}\n'
        msg += f'v: {v:>9}\n'
        msg += f'w: {w:>9}\n\n'
        msg += f'H_bins: {H_bins:>9}\n'
        msg += f'K_bins: {K_bins:>9}\n'
        msg += f'L_bins: {L_bins:>9}\n'
        msg += f'E_bins: {E_bins:>9}\n\n'
        print(msg)

        # call MDNorm to bin MDE workspace
        self.msi.MDNorm(InputWorkspace=MDE_ws,
                        QDimension0=u,
                        QDimension1=v,
                        QDimension2=w,
                        Dimension0Name='QDimension0',
                        Dimension1Name='QDimension1',
                        Dimension2Name='QDimension2',
                        Dimension3Name='DeltaE',
                        Dimension0Binning=H_bins,
                        Dimension1Binning=K_bins,
                        Dimension2Binning=L_bins,
                        Dimension3Binning=E_bins,
                        OutputWorkspace=self.histo_ws_name,
                        OutputDataWorkspace='_d',
                        OutputNormalizationWorkspace='_n')
    
    # ----------------------------------------------------------------------------------------------

    def get_bin_str(self,bins):
        """
        get bin string arg for MDNorm
        """
        if bins is not None:
            _bins = ''
            for _ in bins:
                _bins += f'{_:.5f},'
            bins = _bins.strip(',')
        else:
            bins = ''
        return bins
    
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
# class to calculate and subtract about background
# --------------------------------------------------------------------------------------------------
    
class c_background_tools:

    # ----------------------------------------------------------------------------------------------

    def __init__(self,raw_file):

        """
        class that reads raw file, smooths the data in it, and then subtracts the background from
        Q-points in the raw file by doing "rocking scans" and taking the minimum of each. 
        """

        check_file(raw_file)
        self.raw_file = raw_file

        self.smoothed_file = self._append_suffix_to_file_name(raw_file,suffix='SMOOTHED')
        self.background_file = self._append_suffix_to_file_name(raw_file,suffix='BACKGROUND')
        self.background_subtracted_file = self._append_suffix_to_file_name(
                    raw_file,suffix='BACKGROUND_SUBTRACTED')

    # ----------------------------------------------------------------------------------------------

    def _append_suffix_to_file_name(self,file_name,suffix):

        """
        remove file type from file name and append a suffix to the leading text. put filetype 
        back on
        """

        file_type = file_name.split('.')[-1]
        file_type = '.'+file_type
        file_prefix = file_name[:-len(file_type)]
        new_file_name = file_prefix+'_'+suffix+file_type

        return new_file_name

    # ----------------------------------------------------------------------------------------------

    def make_smoothed_file(self,smoothing_fwhm=None,num_blocks=1):

        """
        split the Q-point set in raw file in num_blocks blocks. for each, block
        interpolate the data onto a fine energy grid and then Gaussian smooth the data. 
        """

        print('\n*** smoothing data ***\n')
        
        check_file(self.raw_file)

        with h5py.File(self.raw_file,'r') as raw_db, \
                h5py.File(self.smoothed_file,'w') as smooth_db:
           
            self.num_Q = raw_db['H_rlu'].size
            self.energy = raw_db['DeltaE'][...]

            if smoothing_fwhm is None:
                self.smoothing_fwhm = (self.energy[1]-self.energy[0])*3.0
            else:
                self.smoothing_fwhm = float(smoothing_fwhm)

            self._copy_header_datasets(raw_db,smooth_db)

            # add new header info
            smooth_db.create_dataset('smoothing_fwhm',data=self.smoothing_fwhm)

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
            signal = raw_db['signal'][Q_indices,:]
            smooth_db['signal'][Q_indices,:] = \
                    self._gaussian_smooth(self.energy,signal,self.smoothing_fwhm)

            _bt.stop()

        _t.stop()

    # ----------------------------------------------------------------------------------------------

    def _copy_header_datasets(self,source_db,dest_db):

        """
        copy 'header' info from database_1 to database_2
        """

        _t = c_timer('copy_datasets')

        print('\ncopying_datasets ...\n')

        keys = list(source_db.keys())
       
        # don't want to copy raw signal. will write new data later 
        keys.pop(keys.index('signal'))
        dest_db.create_dataset('signal',shape=source_db['signal'].shape)

        for key in keys:
            dest_db.create_dataset(key,data=source_db[key])

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

    def calculate_background(self,num_bg_cuts=10,delta_Q_length=0.5,delta_polar_angle=10,
                                     delta_azimuthal_angle=10,num_Q_point_procs=1):

        """
        do a 'rocking scan' around each Q-pt in the raw file and get random neighboring points 
        that lie in a piece of spherical shell the same Q-pt. 
        
        num_bg_cuts determines the number of neighboring points to use.

        the neighboring points are bounded by |Q'| in |Q| +- delta_Q_length (in 1/A) and
        angle' in angle +- delta_angle (in degrees)

        the neighboring points are read from the smoothed file and the background is taken as the
        point-by-point minimum of all of the Q-points (neighboring and original).

        the idea is that the background is spherical, while the single-crystal scattering is 
        not. so looking at data lying on the same sphere but offset by a small angle means 
        we will move away from the phonons and only look at the background. ofcourse it is 
        possible that we might find other phonons at another Q-pt... so we use multiple 
        background cuts that are randomly chosen as a fail safe.
        """
        
        print('\n*** calculating background ***\n')
        
        _t = c_timer('calculate_background')

        check_file(self.raw_file)
        check_file(self.smoothed_file)

        with h5py.File(self.smoothed_file,'r') as _db:
            self.num_Q = _db['H_rlu'].size
            self.energy = _db['DeltaE'][...]

        self.num_E = self.energy.size
        self.background_signal = np.zeros((self.num_Q,self.num_E),dtype=float)

        # split Qpt indices onto multiple process to get and calculate background 
        self.Q_indices_on_procs = np.array_split(np.arange(self.num_Q),num_Q_point_procs)
        print('\nnumber of Q-pts on each proc:')
        msg = ''
        for ii in range(num_Q_point_procs):
            _ = f'proc[{ii}]'
            msg += f'{_:10} {self.Q_indices_on_procs[ii].size}\n'
        print(msg)

        # +- bounds for the polar coordinates
        self.delta_polar_angle = delta_polar_angle*np.pi/180 # degrees => radians
        self.delta_azimuthal_angle = delta_azimuthal_angle*np.pi/180 # degrees => radians
        self.delta_Q_length = delta_Q_length

        # Queue for passing data between procs
        self.queue = mp.Queue()

        proc_list = []
        for pp in range(num_Q_point_procs):
            proc_list.append(mp.Process(target=self._get_background_by_rocking_scans,
                                args=[pp,num_bg_cuts]))

        # start execution
        for pp in proc_list:
            pp.start()

        # get the results from queue
        for pp in range(num_Q_point_procs):

            proc, background_signal_on_proc = self.queue.get()
            Q_indices = self.Q_indices_on_procs[proc]
            self.background_signal[Q_indices,:] = background_signal_on_proc[...]

        # blocking; wait for all procs to finish before moving on
        for pp in proc_list:
            pp.join()
            pp.close()

        self._write_background_file()

        _t.stop()

    # ----------------------------------------------------------------------------------------------

    def _write_background_file(self):

        """
        create background file and copy datasets from raw file into it
        """

        _t = c_timer('write_background_file')

        with h5py.File(self.smoothed_file,'r') as smooth_db, \
                h5py.File(self.background_file,'w') as bg_db:

            self._copy_header_datasets(smooth_db,bg_db)
            bg_db['signal'][...] = self.background_signal[...]

        _t.stop()

    # ----------------------------------------------------------------------------------------------

    def _get_background_by_rocking_scans(self,proc,num_bg_cuts):

        """
        do a 'rocking scan' around each Q-pt in the raw file and get random neighboring points
        that lie in a piece of spherical shell centered Q-pt.

        note, since Qpts are split over procs, only arg we need is which proc this is

        outline:    loop over all Qpts on this proc. 
                    for each one, randomly choose upto num_bg_cuts number of neighboring Qpts. 
                    get the smooth data for all of them 
                    calculate BG as point-by-point minimum in the smooth data
                    fill a BG array for all of the Qpts on this proc.
                    return to main proc using multiprocessing Queue 
                    main proc writes BG arrays to file
        """

        if proc == 0:
            _t = c_timer('get_background_by_rocking_scans')
            print('\nonly printing progress for proc 0\n')

        Q_indices = self.Q_indices_on_procs[proc]
        num_Q_on_proc = Q_indices.size
        background_signal = np.zeros((num_Q_on_proc,self.num_E),dtype=float)

        with h5py.File(self.raw_file,mode='r',libver='latest',swmr=True) as raw_db, \
            h5py.File(self.smoothed_file,mode='r',libver='latest',swmr=True) as smooth_db:

            # need these to find Qpt and neighbors in spherical coords
            self.Q_len = raw_db['Q_len'][...]
            self.polar_angle = raw_db['polar_angle'][...]
            self.azimuthal_angle = raw_db['azimuthal_angle'][...]

            for ii, Q_index in enumerate(Q_indices):

                num_cuts = num_bg_cuts

                if proc == 0:
                    if ii%1000 == 0:
                        print(f'  Qpt {ii} out of {num_Q_on_proc}')
                
                # get this Qpt
                Q_ii = self.Q_len[Q_index]
                polar_ii = self.polar_angle[Q_index]
                azimuthal_ii = self.azimuthal_angle[Q_index]
                
                # get neighbors
                neighbor_inds = self._get_neighboring_Q_points_in_shell(Q_ii,polar_ii,azimuthal_ii)
                num_neighbors = neighbor_inds.size

                # if no neighbors, background is 0
                if num_neighbors == 0:
                    msg = '\n*** WARNING ***\n'
                    msg += f'no bg cuts for Qpt num. {Q_index}!\n'
                    msg += 'setting bg to 0 for this Qpt\n'
                    print(msg)

                    background_signal[ii,:] = np.zeros(self.num_E,dtype=float)

                    continue

                # warn user if not enough neighbors
                if num_neighbors < num_bg_cuts:

                    # get all neighbors if num_neighbors < num_bg_cuts
                    num_cuts = num_neighbors 

                    msg = '\n*** WARNING ***\n'
                    msg += f'found {num_neighbors} bg cuts for Qpt num. {Q_index}\n'
                    print(msg)

                # randomly pick neighbors
                np.random.shuffle(neighbor_inds) # shuffles in place ...
                neighbor_inds = neighbor_inds[:num_cuts]

                background_signal[ii,:] = self._calculate_background_from_pt_by_pt_min(
                        Q_index,neighbor_inds,smooth_db)

        # put in queue to return to main proc
        self.queue.put([proc,background_signal])

        if proc == 0:
            _t.stop()

    # ----------------------------------------------------------------------------------------------

    def _calculate_background_from_pt_by_pt_min(self,Q_index,neighbor_inds,smooth_db):

        """
        cut Q_pt and neighbors from smoothed data and calculate background as the 
        point-by-point minimum of the data
        """

        _mask = 1e8
        
        neighbor_inds = np.sort(neighbor_inds)
        num_neighbors = neighbor_inds.size
        smooth_signal = np.zeros((num_neighbors+1,self.num_E),dtype=float)
        smooth_signal[0,:] = smooth_db['signal'][Q_index,:]
        smooth_signal[1:,:] = smooth_db['signal'][neighbor_inds,:]

        smooth_signal = np.nan_to_num(smooth_signal,nan=_mask,posinf=_mask,neginf=_mask)
        background_signal = np.min(smooth_signal,axis=0)

        background_signal[np.flatnonzero(~(background_signal < _mask))] = 0.0

        return background_signal

    # ----------------------------------------------------------------------------------------------

    def _get_neighboring_Q_points_in_shell(self,Q,polar,azimuthal):

        """
        find neighboring Q-pts that are bounded by |Q|+-dQ, polar+-d_polar, azi+-d_azi,
        i.e. all atoms in the shell centered on the Q-pt defined earlier
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

        # intersection of all of these is the shell
        neighbor_inds = np.intersect1d(len_inds,polar_inds)
        neighbor_inds = np.intersect1d(neighbor_inds,azimuthal_inds)

        return neighbor_inds
    
    # ----------------------------------------------------------------------------------------------

    def subtract_background(self,num_blocks=1):

        """
        open raw and bg files and subtract bg from raw data
        """
        
        print('\n*** subtracting data ***\n')
        
        _t = c_timer('subtract_background')

        with h5py.File(self.raw_file,'r') as raw_db, \
            h5py.File(self.background_file,'r') as bg_db, \
            h5py.File(self.background_subtracted_file,'w') as sub_db:

            self._copy_header_datasets(bg_db,sub_db)
            
            Q_index_blocks = np.arange(self.num_Q)
            Q_index_blocks = np.array_split(Q_index_blocks,num_blocks)

            print(f'\nsplitting data into {num_blocks} blocks\n')

            for block in range(num_blocks):
                
                Q_indices = Q_index_blocks[block]
                raw_signal = raw_db['signal'][Q_indices,:]
                bg_signal = bg_db['signal'][Q_indices,:]
                sub_db['signal'][Q_indices,:] = raw_signal-bg_signal

        _t.stop()


    # ----------------------------------------------------------------------------------------------



# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    _t = c_timer('MDE_tools',units='m')

    # temp and projection
    T = 300

    MDE_file_name = f'../merged_mde/LSNO25_Ei_120meV_300K.nxs'
    out_file_name = f'LSNO25_300K_parallel_test.hdf5'
    
    u = [ 1, 0, 0]
    v = [ 0, 1, 0]
    w = [ 0, 0, 1]

    #H_bins = [  -5.025,  0.05,  15.025]
    #K_bins = [ -12.025,  0.05,   7.525]
    #L_bins = [ -11.000,   2.0,  11.000]
    #E_bins = [  -20.25,   0.5,  100.25]

    H_bins = [  -0.050,  0.10,   8.050]
    K_bins = [  -0.050,  0.10,   8.050]
    L_bins = [  -2.500,   5.0,   2.500]
    E_bins = [  -20.25,   0.5,  100.25]

    num_chunks = [4,4,1]

    # class to do the stuff
    MDE_tools = c_MDE_tools(MDE_file_name)
    MDE_tools.bin_MDE_chunks(H_bins,K_bins,L_bins,E_bins,num_chunks,out_file_name,u,v,w)

    _t.stop()









