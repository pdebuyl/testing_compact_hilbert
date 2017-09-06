Testing the compact Hilbert curve
=================================

This code compares the implementation of the compact Hilbert index across Chris
Hamilton's original code and my Python version.

Authors and licenses
--------------------

1. For the submodule `libhilbert`: Chris Hamilton, LGPLv2, see
   https://github.com/pdebuyl/libhilbert
   This submodule is not *distributed* here, you have to fetch the submodule.
2. For the file `hilbert.cxx`: Peter Colberg, MIT, see
   https://colberg.org/nano-dimer/
3. For the rest: Pierre de Buyl, BSD.

Usage
-----

1. Requirements:
    - A C++ compiler
	- Python
    - [HDF5](https://support.hdfgroup.org/HDF5/)
	- [NumPy](http://www.numpy.org/)
	- [h5py](http://www.h5py.org/)
	
2. Setup:
    1. Build libhilbert

            cd libhilbert
            ./configure
            make
            cd ..

    2. Know the location to HDF5, if is not known yet.

2. Run the default comparisons.

    make HDF5_HOME=/home/hdf5-1.8.18

(adjust the location of HDF5 accordingly).

This will compare the Python version and the C++ version of the code after build the
executable `hilbert`.
