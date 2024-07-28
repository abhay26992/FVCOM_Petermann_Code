<div align="justify"> 
  
# FVCOM_Petermann_Code

The Finite Volume Community Ocean Model (FVCOM) ([Chen et al., 2007](https://doi.org/10.1029/2006JC003485)) has recently been amended by an ice shelf module ([Zhou and Hattermann, 2020](https://doi.org/10.1016/j.ocemod.2019.101536)) that allows modelling of oceanic processes in ice shelf cavities bounded by complex coastal geometries and fjord bathymetry. A new sea ice module (<em>Ice Nudge</em>) has been implemented into FVCOM that allows the user to prescribe sea ice variables (sea ice concentration and thickness, bulk ice salinity, and sea ice velocities) as external surface boundary conditions. The <em>Ice Nudge</em> module reads in the sea ice concentration to determine where to modify the ocean momentum and thermodynamical surface fluxes. A technical description of the module has been detailed in [Prakash et al. (2022)](https://doi.org/10.1016/j.mex.2022.101668). This repository contains version 4.0 (v4.0) of the FVCOM code, augmented by both the ice shelf and sea ice modules, which was used to conduct the [sea ice arches](https://doi.org/10.5194/tc-17-5255-2023) and the subglacial discharge experiments for the Petermann ice shelf and fjord. 

# Running the subglacial discharge experiments

We have developed a high-resolution (200 m) unstructured grid 3-D regional FVCOM setup centered over the Petermann ice shelf and fjord, which utilizes a realistic seafloor and ice shelf basal topography ([Prakash et al., 2022](https://doi.org/10.1016/j.mex.2022.101668), [2023](https://doi.org/10.5194/tc-17-5255-2023)). Moreover, the poorly constrained sub-ice shelf bathymetry in the [BedMachine v3](https://doi.org/10.1002/2017GL074954) dataset has been remedied, together with the inclusion of an [inner sill](https://doi.org/10.1016/j.epsl.2015.04.009), to provide a [more accurate representation of the water column thickness beneath Petermann's floating tongue](https://doi.org/10.1016/j.mex.2022.101668). A novel study investigating the mechanisms behind subglacial discharge-driven increase in melt and their evolution with increasing discharge in a warmer future climate is presently under consideration for peer-review and publication in a suitable scientific journal. Using the model code provided in this repository together with the datasets provided [here](10.5281/zenodo.12803094), it is possible to reproduce the results of our numerical experiments. 

# Compilation

Please checkout the <em>betzyfv4fisoc</em> branch for instructions regarding model compilation. Note that you may need to contact your HPC support.  

</div>
