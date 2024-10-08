
"""
this file is part of the phonon explorer package!
author: tyler c. sterling & dmitry reznik.
email: ty.sterling@colorado.edu
affil: University of Colorado Boulder, Raman Spectroscopy and Neutron Scattering Lab
date: 08/29/2024
description:
    take binning and projection args from the user and call mantid MDNorm (on the backend)
    to "prebin" data for later analysis with phonon explorer. can be split into "chunks" to 
    save memory. see the description of the variables below.
"""

from file_tools.m_file_utils import *
from file_tools.m_symmetry import check_symmetry_args

import numpy as np
import h5py 
import os

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
        merged_file_name='merged_histo.hdf5',u=[1,0,0],v=[0,1,0],w=None,append=False,
        SymmetryOperations=None):

        """
        split requested binning into chunks and bin over small chunks separately. merge the 
        results of them all into a single file

        H_bin = [start, step, stop]. note, these are now interpreted like the args to mantid 
            MDNorm and to Horace cut_sqw: the binning will actually start at 'start' with spacing
            equal to the step size. i.e. the bin centers will be [start+1*step/2, start+2*step/2, 
            start+3*step/2, ...]

        etc. for K_bin, L_bin, E_bin

        num_chunks is number of chunks to split binning in Q-space. E.g. if num_chunks = [2,1,1]
        and there are 10 bins along H, the data will be cut in 2 go's with the first 5 H-bins in
        the first go and the last 5 in the seconds. 

        merged_file_name is the *.hdf5 where all bins will be written.

        u, v, w are the projections. w is optional and the cross product of u and v if not given.
        """
        
        _t = c_timer('bin_MDE_chunks',units='m')

        self.SymmetryOperations = SymmetryOperations
        self.merged_file_name = merged_file_name

        # overwrite by default, but dont overwrite in some cases, e.g. assembling multiple offsets
        _exists = os.path.exists(merged_file_name)
        if append:
            if not _exists:
                msg = f'creating file \'{self.merged_file_name}\'\n'
                print(msg)
            else:
                msg = f'appending to file \'{self.merged_file_name}\'\n'
                print(msg)
        else:
            if _exists:
                msg = f'file \'{self.merged_file_name}\' already exists. removing it ...\n'
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

        H_bin_centers = np.round((self.H_bin_edges[:-1]+self.H_bin_edges[1:])/2,3)
        K_bin_centers = np.round((self.K_bin_edges[:-1]+self.K_bin_edges[1:])/2,3)
        L_bin_centers = np.round((self.L_bin_edges[:-1]+self.L_bin_edges[1:])/2,3)
        E_bin_centers = np.round((self.E_bin_edges[:-1]+self.E_bin_edges[1:])/2,3)
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
                            self.u,self.v,self.w,SymmetryOperations=self.SymmetryOperations)

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

    def get_sparse_arrays_from_histo_ws(self):
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

    def bin_MDE(self,H_bins=None,K_bins=None,L_bins=None,E_bins=None,u=[1,0,0],v=[0,1,0],w=None,
                SymmetryOperations=None):
        """
        bin the events in the MDE workspace into a histogram workspace; note that this doesnt 
        really return anything or create new attributes. the produced data are stored in the 
        histogram workspace

        H_bin = [start, step, stop]. note, these are now interpreted like the args to mantid
            MDNorm and to Horace cut_sqw: the binning will actually start at 'start' with spacing
            equal to the step size. i.e. the bin centers will be [start+1*step/2, start+2*step/2,
            start+3*step/2, ...]

        etc. for K_bin, L_bin, E_bin

        u, v, w are the projections. w is optional and the cross product of u and v if not given.
        """

        if not hasattr(self,'H_range'):
            self.H_bin_edges, self.H_range, self.dH = self._get_bin_edges(H_bins)
            #self.K_bin_args = K_bin_args
            self.K_bin_edges, self.K_range, self.dK = self._get_bin_edges(K_bins)
            #self.L_bin_args = L_bin_args
            self.L_bin_edges, self.L_range, self.dL = self._get_bin_edges(L_bins)
            #self.E_bin_args = E_bin_args
            self.E_bin_edges, self.E_range, self.dE = self._get_bin_edges(E_bins)

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

        print('\nSymmetryOperations:')
        print(SymmetryOperations,'\n')

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
                        OutputNormalizationWorkspace='_n',
                        SymmetryOperations=SymmetryOperations)
    
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

def bin_MDE(MDE_file_name,H_lo,H_hi,H_step,K_lo,K_hi,K_step,L_lo,L_hi,L_step,E_lo,E_hi,E_step,
            H_bin=1,K_bin=1,L_bin=1,merged_file_name='merged_file.hdf5',
            u=[1,0,0],v=[0,1,0],w=None,num_chunks=[1,1,1],SymmetryOperations=None):

    """
    wrapper to translate variable names for interface with c_MDE_tools.bin_MDE_chunks
    to more sensible names for the user.
    
    MDNorm_kw_args are any bonus arguments you want to pass to MDNorm. Note, you have to give 
    them as kw-args in function call:
        ex. bin_NXSPE(...,SymmetryOperations='x,y,z;-x-y-z')
    ...

    """

    timer = c_timer('bin_MDE',units='m')

    check_symmetry_args(SymmetryOperations)

    # class to do all of the stuff
    MDE_tools = c_MDE_tools(MDE_file_name)

    H_bin = int(H_bin); K_bin = int(K_bin); L_bin = int(L_bin)
    if H_bin < 1: H_bin = 1
    if K_bin < 1: K_bin = 1
    if L_bin < 1: L_bin = 1

    # convert binning args to ones used by c_MDE_tools
    H_bins = _get_bins(H_lo,H_hi,H_step*H_bin)
    K_bins = _get_bins(K_lo,K_hi,K_step*K_bin)
    L_bins = _get_bins(L_lo,L_hi,L_step*L_bin)
    E_bins = _get_bins(E_lo,E_hi,E_step)
    print('H_bins:',H_bins)
    print('K_bins:',K_bins)
    print('L_bins:',L_bins)
    print('E_bins:',E_bins,'\n')

    # convert shift args to useful form
    shifts = _setup_shifts(H_step,K_step,L_step,H_bin,K_bin,L_bin)
    if shifts.size == 0:
        loop_over_shifts = False
    else:
        loop_over_shifts = True

    # bin the data on the unshifted grid first. this initiates file etc.
    MDE_tools.bin_MDE_chunks(H_bins,K_bins,L_bins,E_bins,num_chunks,merged_file_name,u,v,w,
                             SymmetryOperations=SymmetryOperations)

    # only loop over offsets if atleast one is defined
    if loop_over_shifts:

        num_shifts = shifts.shape[0]

        print('\n----------------------------------------------------------------\n')
        print('num_shifts:',num_shifts)
        print('shifts:')
        print(shifts,'\n')

        msg = 'looping over shifts'
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
            MDE_tools.bin_MDE_chunks(H_bins+_H_shift,K_bins+_K_shift,L_bins+_L_shift,E_bins,
                num_chunks,merged_file_name,u,v,w,append=True,
                SymmetryOperations=SymmetryOperations)

    timer.stop()

# --------------------------------------------------------------------------------------------------    

    







