
import os
from shutil import copyfile
import subprocess


clean = True
#clean = False

#BuildPreReq = True
BuildPreReq = False

BuildFVCOM = True
#BuildFVCOM = False

COLD = "/projappl/project_2000339/"

prereq_install_dir =  "/projappl/project_2000339/installs/FVCOM_external/"
install_dir =  "/projappl/project_2000339/installs/FVCOM_FX4/"

prereq_dir =  "/projappl/project_2000339/source/fvcom4_fisoc/FVCOM_source/libs/"
source_dir =  "/projappl/project_2000339/source/fvcom4_fisoc/FVCOM_source/"

prereq_makefile_name = "makefile.puhti"
make_inc_name = "make.inc.puhti"
makefile_name = "makefile.api"
makedep_name = "makedepends.api"

# environment variables set here will be used over vars set with ?= in the makefile
os.environ["IOLIBS"] = "-L/appl/spack/install-tree/intel-19.0.4/netcdf-fortran-4.4.4-nbpz5p/lib -lnetcdff -lnetcdf"
os.environ["IOINCS"] = "-I/appl/spack/install-tree/intel-19.0.4/netcdf-fortran-4.4.4-nbpz5p/include"

os.environ["INSTALLDIR"] = prereq_install_dir
os.environ["TOPDIR"] = source_dir
           
if clean:
    make_args = "clean"
else:
    make_args = ""


print("\nBuilding FVCOM pre-requisites\n")

try:
    os.chdir(source_dir)
    copyfile(make_inc_name, "make.inc")
    copyfile(makefile_name, "makefile")
    copyfile(makedep_name, "makedepends")
    os.chdir(prereq_dir)
    copyfile(prereq_makefile_name, "makefile")
except:
    raise NameError("ERROR: Failed to replace makefiles\n")

if BuildPreReq:
    if clean:
        ret = subprocess.call(["make","clean"])
    else:
        ret = subprocess.call(["make"])

    if (ret != 0):
        raise NameError("ERROR: Failed to build FVCOM pre-requisites\n")
    

print("\nBuilding FVCOM \n")

os.chdir(source_dir)

if clean:
    print "cleaning"
    ret = subprocess.call(["make", "clean"])
    if (ret != 0):
        raise NameError("ERROR: Failed to build FVCOM\n")
    ret = subprocess.call(["make"])
else:
    ret = subprocess.call(["make"])
if (ret != 0):
    raise NameError("ERROR: Failed to build FVCOM\n")

print "copy to ",install_dir

try:
    copyfile("fvcom_driver", install_dir+"/fvcom_driver")
    copyfile("libfvcom_api.so", install_dir+"/libfvcom_api.so")
except:
    raise NameError("ERROR: Failed to copy lib and exe to install_dir\n")


print("\nFVCOM build complete (maybe it even works...) \n")

