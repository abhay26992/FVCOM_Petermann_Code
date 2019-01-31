# ESMF-FVCOM-FISOC

Developing source code on coupling FVCOM with the ice sheet model Elmer/ICE through a recently developed Framework for Ice Shelf Ocean Coupling , build on Earth System Modeling Framework. 

To build on virtual machine / elemeruser prepare code by:
- copying FVCOM_source/make.inc.elmeruser to FVCOM_source/make.inc
- copying FVCOM_source/libs/makefile.vm_elemeruser to FVCOM_source/libs/makefile
- and execute it to build local libraries

Then build FVCOM via cmake following the instructions given in the readme.elmeruser of the uk-fvcom-iceshelf-vm-elmeruser branch of the fvcom-cmake repository