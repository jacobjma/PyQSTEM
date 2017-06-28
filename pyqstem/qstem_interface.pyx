import warnings
from libcpp.vector cimport vector
import numpy as np
from collections import OrderedDict
import matplotlib.pyplot as plt
from libcpp.unordered_map cimport unordered_map
from libcpp.pair cimport pair
from libcpp.string cimport string
from .wave import Wave, Potential
from .util import atoms_plot

TOL=1e-6

# c++ interface to cython
cdef extern from "QSTEM.h" namespace "shapes":
  cdef cppclass QSTEM:
        QSTEM(string) except +
        string mode

        int trans_array_state
        int wave_state
        int atoms_state
        int box_state
        int scan_range_state

        void set_atoms(int,int,vector[vector[double]] &,vector[double] &,
                        vector[double] &,vector[int] &)
        void set_positions(int,vector[vector[double]] &)
        void set_box(vector[double] &,int,int,float)

        void get_resolution(float*, float*)
        void get_potential_samples(int*, int*, int*)
        void get_probe_samples(int*, int*)
        void get_scan_range(float*, float*, int*, float*, float*, int*)
        void get_energy(float*)
        void get_potential_extent(float*,float*,float*,float*)

        void build_potential(int)
        void set_potential(vector[vector[vector[vector[double]]]] &,vector[int] &,vector[double] &)
        vector[vector[vector[vector[double]]]] get_potential_or_transfunc(float*,float*,float*)
        void calculate_transfunc()

        void set_scan_range(float,float,int,float,float,int)
        void build_wave(int,float,int,int,float,float)
        void build_probe(float,float,int,int,float,float,unordered_map[string,float])
        void set_wave(vector[vector[vector[double]]],float,int,int,float,float)
        vector[vector[vector[double]]] get_wave(float*,float*,float*)

        void create_detectors(vector[string],vector[vector[double]],int)
        void run(int)
        vector[vector[double]] read_detector(int)

# cython wrapper class
cdef class PyQSTEM:

    cdef object _atoms
    cdef char* c_mode

    cdef QSTEM *thisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self, str mode):
        self._atoms=None

        if isinstance(mode, unicode):
            c_mode = mode.encode('UTF-8')
        else:
            c_mode = mode

        self.thisptr = new QSTEM(c_mode)

        self.thisptr.mode = c_mode;

        # 0: potential not created
        # 1: potential created
        # 2: potential converted to transmission function
        # -1: potential created, but have to be recreated (atoms or simulation box changed)
        # -2: potential created, but have to be recreated (energy changed)
        self.thisptr.trans_array_state=0

        # 0: wave not created
        # 1: wave created
        # -1: wave created, but have to be recreated
        self.thisptr.wave_state=0

        # 0: atoms set
        # 1: atoms not set
        self.thisptr.atoms_state=0

        # 0: simulation box set
        # 1: simulation box not set
        self.thisptr.box_state=0

    _detectors=OrderedDict()

    def __dealloc__(self):
        del self.thisptr

    def set_atoms(self,atoms,set_box=True,TDS=False,B_factors=None):
        self._atoms=atoms

        positions=atoms.get_positions()

        occ=[1.]*len(atoms) # occupancy not implemeted
        q=[0.]*len(atoms) # charge not implemeted

        Znum=atoms.get_atomic_numbers()

        if set_box:
            box=np.diag(atoms.get_cell())
            self.set_box(box)

        self.thisptr.set_atoms(len(atoms),len(set(Znum)),positions,occ,q,Znum)
        self.thisptr.atoms_state = 1
        if self.thisptr.trans_array_state >= 1:
            self.thisptr.trans_array_state = -1

    def _set_positions(self,positions):
        self.thisptr.set_positions(len(self._atoms),positions)

    def set_box(self,box,periodic_xy=None,periodic_z=False,cell_div=1):
        if ((self.thisptr.mode == 'STEM')&(periodic_xy==True)):
            raise RuntimeError('Periodic boundary conditions not implemented for mode {0}'.format(self.thisptr.mode))

        if periodic_xy is None:
            if self.thisptr.mode == 'STEM':
                periodic_xy = False
            elif self.thisptr.mode == 'TEM':
                periodic_xy = True

        nonperiodic_xy = not periodic_xy
        nonperiodic_z = not periodic_z

        self.thisptr.set_box(box,nonperiodic_xy,nonperiodic_z,cell_div)
        self.thisptr.box_state = 1
        if self.thisptr.trans_array_state >= 1:
            self.thisptr.trans_array_state = -1

    def get_resolution(self):
        cdef float resolutionX, resolutionY
        self.thisptr.get_resolution(&resolutionX, &resolutionY)
        return resolutionX, resolutionY

    def get_potential_samples(self):
        cdef int potNx, potNy, slices
        self.thisptr.get_potential_samples(&potNx, &potNy, &slices)
        return potNx, potNy, slices

    def get_probe_samples(self):
        cdef int nx, ny
        self.thisptr.get_probe_samples(&nx, &ny)
        return nx, ny

    def get_scan_range(self):
        cdef float scanXStart, scanXStop, scanYStart, scanYStop
        cdef int scanXN, scanYN

        self.thisptr.get_scan_range(&scanXStart, &scanXStop, &scanXN, &scanYStart, &scanYStop, &scanYN)
        return ((scanXStart, scanXStop, scanXN), (scanYStart, scanYStop, scanYN))

    def get_energy(self):
        cdef float v0
        self.thisptr.get_energy(&v0)
        return v0

    def get_probe_extent(self):
        probe_samples = self.get_probe_samples()
        resolution = self.get_resolution()

    def get_potential_extent(self):
        cdef float potOffsetX, potOffsetY, potSizeX, potSizeY
        self.thisptr.get_potential_extent(&potOffsetX, &potOffsetY, &potSizeX, &potSizeY)
        return potOffsetX, potOffsetX+potSizeX, potOffsetY, potOffsetY+potSizeY

    def get_minimum_potential_extent(self):
        scan_range = self.get_scan_range()
        probe_samples = self.get_probe_samples()
        resolution = self.get_resolution()

        probe_size=[resolution[0]*probe_samples[0],resolution[1]*probe_samples[1]]
        return [scan_range[0][0]-.5*probe_size[0],scan_range[0][1]+.5*probe_size[0],
                scan_range[1][0]-.5*probe_size[0],scan_range[1][1]+.5*probe_size[0]]

    def view(self,ax=None):
        scan_range=self.get_scan_range()
        potential_extent=self.get_potential_extent()
        potential_samples=self.get_potential_samples()
        probe_extent=self.get_minimum_potential_extent()
        atoms_plot(self._atoms,scan_range=scan_range,potential_extent=potential_extent,probe_extent=probe_extent,ax=ax,legend=True)

    def build_potential(self,num_slices=None,slice_thickness=.5,scan_range=None,probe_position=None):

        if self.thisptr.mode == 'STEM':
            if scan_range is None:
                raise RuntimeError('Provide scan window for mode STEM')
            else:
                self.set_scan_range(scan_range)
        elif self.thisptr.mode == 'CBED':
            if probe_position is None:
                raise RuntimeError('Provide probe postion for mode CBED')
            else:
                scan_range = [[probe_position[0],probe_position[0],1],[probe_position[1],probe_position[1],1]]
                self.set_scan_range(scan_range)

        if ((self.thisptr.atoms_state==0)|(self.thisptr.box_state==0)):
            raise RuntimeError('Please set atoms and simulation box')

        if num_slices is None:
            num_slices = int(self._atoms.get_cell()[2,2]/slice_thickness)
        #nx_old,ny_old,slices_old = self.get_numsamples()

        self.thisptr.build_potential(num_slices)

        self.thisptr.trans_array_state = 1

        #if ((self.thisptr.wave_state == 1)&((nx_old!=shape[0])|(ny_old!=shape[1]))):
        #    self.thisptr.wave_state = -1

    def set_potential(self,potential,scan_range=None):
        if self.thisptr.mode == 'STEM':
            if scan_range is None:
                raise RuntimeError('Provide scan window for mode {0}'.format(self.thisptr.mode))
            else:
                self.set_scan_range(scan_range)

        array = potential.array
        nonperiodic_xy = not potential.periodic_xy
        nonperiodic_z = not potential.periodic_z
        size = array.shape
        extent = potential.get_extent()

        if self.thisptr.wave_state > 0:
            if not np.any(np.isclose(potential.sampling[:2],self.get_resolution(),rtol=1e-6,atol=1e-6)):
                raise RuntimeError('Potential resolution does not match wavefunction ({0},{1})!=({2},{3})'
                          .format(potential.sampling[0],potential.sampling[1],self.get_resolution()[0],self.get_resolution()[1]))

        array = np.concatenate((np.real(array)[:,:,:,np.newaxis],np.imag(array)[:,:,:,np.newaxis]),axis=3)

        self.thisptr.set_potential(array,size,extent)

        self.thisptr.trans_array_state = 1
        #if self.thisptr.wave_state == 1:
        #    self.thisptr.wave_state = -1

    def get_potential_or_transfunc(self):
        cdef float resolutionX, resolutionY, sliceThickness
        if self.thisptr.trans_array_state==0:
            raise RuntimeError('A potential have not been build')
        else:
            trans_array = self.thisptr.get_potential_or_transfunc(&resolutionX,&resolutionY,&sliceThickness)
            #trans_array = np.apply_along_axis(lambda args: [complex(*args)], 3, trans_array)[:,:,:,0]
            trans_array_noncplx = np.array(trans_array)
            trans_array = np.empty(trans_array_noncplx.shape[:-1], dtype=np.complex)
            trans_array.real = trans_array_noncplx[...,0]
            trans_array.imag = trans_array_noncplx[...,1]
            potential_extent = self.get_potential_extent()
            return Potential(trans_array,(resolutionX,resolutionY,sliceThickness))

    def calculate_transfunc(self):

        if self.thisptr.trans_array_state==0:
            raise RuntimeError('Please set or build a potential')
        elif self.thisptr.trans_array_state==1:
            self.thisptr.calculate_transfunc()
            self.thisptr.trans_array_state=2
        elif self.thisptr.trans_array_state==2:
            pass

    def set_scan_range(self,scan_range):

        self.thisptr.set_scan_range(scan_range[0][0],scan_range[0][1],scan_range[0][2],
                      scan_range[1][0],scan_range[1][1],scan_range[1][2])

    def build_probe(self,v0,alpha,num_samples,resolution=None,window_size=None,
                    defocus=0,Cs=0,C5=0,astig_mag=0,astig_angle=0,aberrations={}):

        if ((resolution is None)&(window_size is None)):
            raise RuntimeError("Please specify resolution or window size")

        if resolution is None:
            resolution = np.array(window_size)/np.array(num_samples)

        cdef unordered_map[string, float] aberrations_map

        symbols=[b'a22',b'phi22',b'a20',
                 b'a33',b'phi33',b'a31',b'phi31',
                 b'a44',b'phi44',b'a42',b'phi42',b'a40',
                 b'a55',b'phi55',b'a53',b'phi53',b'a51',b'phi51',
                 b'a66',b'phi66',b'a64',b'phi64',b'a62',b'phi62',b'a60']

        keys = [key.encode('utf-8') for key in aberrations.keys()]

        for symbol in symbols:
            if symbol in keys:
                aberrations_map[symbol] = aberrations[symbol.decode('utf-8')]
            elif symbol == b'a20':
                aberrations_map[symbol] = defocus
            elif symbol == b'a40':
                aberrations_map[symbol] = Cs
            elif symbol == b'a60':
                aberrations_map[symbol] = C5
            elif symbol == b'a22':
                aberrations_map[symbol] = astig_mag
            elif symbol == b'phi22':
                aberrations_map[symbol] = astig_angle
            else:
                aberrations_map[symbol] = 0.

        self.thisptr.build_probe(v0,alpha,num_samples[0],num_samples[1],resolution[0],resolution[1],aberrations_map)
        self.thisptr.wave_state = 1

    def build_wave(self,str wave_type,v0,num_samples,resolution=None):

        if wave_type=='plane':
            wave_type_int=0
        else:
            raise RuntimeError('Wave type {0} not recognized'.format(wave_type))

        old_v0 = self.get_energy()

        if resolution is None:
            resolution = [-1,-1]

        self.thisptr.build_wave(wave_type_int,v0,num_samples[0],num_samples[1],resolution[0],resolution[1])

        self.thisptr.wave_state = 1
        if self.thisptr.trans_array_state == 2:
            if not np.isclose(old_v0,v0,rtol=0,atol=1e-6):
                self.thisptr.trans_array_state = -2

    def set_wave(self,wave):
        nx=wave.array.shape[0]
        ny=wave.array.shape[1]

        nx_old,ny_old = self.get_probe_samples()
        resolutionX,resolutionY = self.get_resolution()

        if self.thisptr.trans_array_state>0:
            if ((nx_old!=nx)|(ny_old!=ny)):
                raise RuntimeError('Wave function shape does not match QSTEM: ({0},{1}) != ({2},{3})'.format(nx,ny,nx_old,ny_old))

        #if wave.sampling is not None:
        #    if not np.any(np.isclose(wave.sampling,(resolutionX,resolutionY),rtol=0,atol=1e-06)):
        #        warnings.warn('Wavefunction resolution will be changed to match the simulation box')

        array = np.concatenate((np.real(wave.array)[:,:,np.newaxis],np.imag(wave.array)[:,:,np.newaxis]),axis=2)

        if wave.sampling is None:
            resolution = [-1,-1]
        else:
            resolution = wave.sampling

        old_v0 = self.get_energy()
        self.thisptr.set_wave(array,wave.energy,nx,ny,resolution[0],resolution[1])

        self.thisptr.wave_state = 1
        if self.thisptr.trans_array_state == 2:
            if not np.isclose(old_v0,wave.energy,rtol=0,atol=1e-6):
                self.thisptr.trans_array_state = -2

    def get_wave(self):
        if self.thisptr.wave_state == 0:
            raise RuntimeError('A wavefunction have not been created')

        cdef float resolutionX, resolutionY, v0
        wave_array = self.thisptr.get_wave(&resolutionX,&resolutionY,&v0)
        wave_array = np.apply_along_axis(lambda args: [complex(*args)], 2, wave_array)[:,:,0]

        return Wave(wave_array,energy=v0,sampling=(resolutionX,resolutionY))

    def add_detector(self,str name,radii,shift=(0,0)):
        self._detectors[name]=[False,radii[0],radii[1],shift[0],shift[1]]

    def remove_detector(self,str name):
      if name not in self._detectors.keys():
          raise RuntimeError('Detector name {0} not recognized'.format(name))

      del self._detectors[name]

    def read_detector(self,str name,dwell_time=1,current=1):

        if name not in self._detectors.keys():
            raise RuntimeError('Detector name {0} not recognized'.format(name))

        if not self._detectors[name][0]:
            raise RuntimeError('Detector {0} is empty, add the detector before running'.format(name))

        img = np.array(self.thisptr.read_detector(list(self._detectors.keys()).index(name)))
        img*=dwell_time*current/1.6021773e-4

        return img

    def run(self,slices_per_division=None,slice_thickness=.5,cell_divisions=1,scan_range=None,probe_position=None,display_progress_interval=None):

        if self.thisptr.wave_state==0:
            raise RuntimeError('A wavefunction have not been created')
        elif self.thisptr.wave_state==-1:
            raise RuntimeError('The wavefunction have to be recreated after changing the potential sampling')
        elif self.thisptr.trans_array_state==-1:
            raise RuntimeError('The potential have to be recreated after changing the atoms or simulation box')
        elif self.thisptr.trans_array_state==-2:
            raise RuntimeError('The transmission function energy does not match wavefunction')

        cdef unordered_map[string, vector[double]] detector_map

        if self.thisptr.mode == 'STEM':
            pot_ext = self.get_potential_extent()
            min_pot_ext = self.get_minimum_potential_extent()

            if ((pot_ext[0]-min_pot_ext[0]>TOL)|(min_pot_ext[1]-pot_ext[1]>TOL)|
                (pot_ext[2]-min_pot_ext[2]>TOL)|(min_pot_ext[3]-pot_ext[3]>TOL)):
                raise RuntimeError('The potential is too small to accomodate the probe size for the chosen scan range')

            properties=[]
            for name,values in self._detectors.items():
                values[0]=True
                properties.append(values[1:])

            names = [name.encode('utf-8') for name in self._detectors.keys()]
            self.thisptr.create_detectors(names,properties,len(self._detectors))

        if display_progress_interval is None:
            display_progress_interval = -1

        box=np.diag(self._atoms.get_cell())
        h=box[2]/cell_divisions
        if slices_per_division is None:
            slices_per_division = int(h/slice_thickness)

        if cell_divisions == 1:
            if self.thisptr.trans_array_state==0:
                self.build_potential(slices_per_division,scan_range=scan_range,probe_position=probe_position)

            self.calculate_transfunc()
            self.thisptr.run(display_progress_interval)
        else:
            box_div=[box[0],box[1],h]
            self.set_box(box_div)
            positions=self._atoms.get_positions()

            for i in range(cell_divisions):
                self._set_positions(positions-[0,0,i*h])

                self.build_potential(slices_per_division,scan_range=scan_range,probe_position=probe_position)
                self.calculate_transfunc()
                self.thisptr.run(display_progress_interval)
