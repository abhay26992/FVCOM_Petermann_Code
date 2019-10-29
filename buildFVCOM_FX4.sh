
module purge

module load intel/19.0.4 intel-mpi/18.0.5 intel-mkl/2019.0.4
module load esmf
module load elmer/fisoc
module load fvcom/FX4
module load cmake
module load netcdf/4.7.0
module load netcdf-fortran/4.4.4
module list

python buildFVCOM_FX4.py

#export FC=mpiifort
#export CFLAGS=" -O3 -fPIC"
#export OPT=-O3
#export CPP= /usr/bin/cpp
#export COMPILER= -DMPIIFORT  


