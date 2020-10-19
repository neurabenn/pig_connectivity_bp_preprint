**This is the git repository contains the files of the [Pig WM atlas and Connectivity Blueprint Preprint](https://doi.org/10.1101/2020.10.13.337436)**

Files in this repository provide the code and template files to recreate the figures in the paper. 

The raw and perprocessd diffision data is available on  [Zenodo](https://zenodo.org/record/4103368#.X43aD0IzZFQ)

The dense connectomes for the cross-species browser will be made available through a release on Xenodo with the raw data files. 
However, users can also recreate the browser by running the associated notebooks.

The contents of each folder are as follows 
* IDW_scripts
    - scripts to perform cross-species spatial prediction with Inverse Distance Weighted Interpolation as demonstrated in figures 9 and 10.
* PNI50_surfaces_FSICO5
    - The Freesurfer compatible average surfaces in Icosahedron 5 space. 
    - *this is not the space of the final blueprint* but can be used for spherical registration 
* T1_templatePNI50
    - The T1 weighted full and brain extracted template as well and segmentation priors. 
* Atlases 
    - The WM and adapted Saikali et al. Atlases. 
        * for use in FSLeyes add to $FSLDIR/data/atlases
* blueprints
    - .mat files corresponding to the connectivity blueprints of each species. 
    - All human blueprints were sourced from Mars et al, 2018 
        * https://github.com/elifesciences-publications/rmars-comparing-connectivity-blueprints   
    - Pig blueprints were generated as shown in figure 1. 
* data_driven_tracts 
    - Left and Right results of the data-driven tractography. 
    - These components were used to guide the definition of tractography protocols
    - Raw files refer to minimal processing of the ICAs where as gmm file have had a threshold applied via gaussian mixture modeling 
    - see: https://www.sciencedirect.com/science/article/pii/S1053811917304044
* Protocols
    - Contains the tractography protocols as defined in PNI50 space for use with the [Xtract](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/XTRACT) package. 
* Tract Projections 
    - Tractograms showing each tracts termiantino points on the cortical surface. 
    - see spec files
* Tracts
    - Individual Nifti files of the tracts in the WM atlas
* Surfaces
    - Surfaces used in the generation of the connectivity blueprint. 
    - Human surfaces are from the original Mars et al. [blueprint paper](https://github.com/elifesciences-publications/rmars-comparing-connectivity-blueprints)
    - Native pig surfaces to the PNI50 are in the pig_surfaces/Surf_Volume_geometry folder
    - The pig surfaces called in the notebooks and spec files have been oriented with a corrected anterior-posterior and inferior-superior labelling.
* Cross_species_browser
    - The cross species browser can be accessed through the cross species spec files
    - The actual files dense connectome images are too big for git ans are hosted on [Zenodo](https://zenodo.org/record/4103368#.X43aD0IzZFQ)
        - Simply place the L_hum_pig_double_click.dconn.nii and R_hum_pig_double_click.dconn.nii files in this folder once downloaded and open the spec file. 
