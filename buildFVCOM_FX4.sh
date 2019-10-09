cd /projappl/project_2000339/source/esmf-fvcom-fisoc/FVCOM_source
cp makefile.api makefile
cp make.inc.puhti make.inc
#export FC=mpiifort
#export CFLAGS=" -O3 -fPIC"
#export OPT=-O3
#export CPP= /usr/bin/cpp
#export COMPILER= -DMPIIFORT  
cp /projappl/project_2000339/source/esmf-fvcom-fisoc/FVCOM_source/libfvcom_api.so /projappl/project_2000339/installs/FVCOM_FX4/
cp /projappl/project_2000339/source/esmf-fvcom-fisoc/FVCOM_source/fvcom_driver /projappl/project_2000339/installs/FVCOM_FX4/

make clean
make

