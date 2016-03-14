UP2D
====

Two dimensional openMP code (spectral, volume penalization)

This is a slender code for medium sized 2D simulations (typically 1024^2 -- 4096^2)
It uses the same mathematics as the 3D flusi code, but it is drastically simplified.

It has the following features
* primitive variables (u,p) penalized Navier--Stokes
* time stepping is RungeKutta 2 with integrating factor
* shared memory parallelism
* FFTW (fftw3_threads) used for FFTs
* hdf5 output, compatible with FLUSI code
* you can use mpi2vis to load the hdf5 data to paraview


Missing features 
* drag/lift/moment/power computation
* dynamic mean flow
* fluid--structure interaction

We also have a MATLAB version of this code, for even faster prototyping.
