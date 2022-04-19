**This is the git repository contains the files of the [Pig WM atlas and Connectivity Blueprint Preprint](https://www.biorxiv.org/content/10.1101/2020.10.13.337436v4.full)**

Files in this repository provide the code and template files to recreate the figures in the paper. 

The raw and perprocessd diffusion data is available on  [Zenodo](https://zenodo.org/record/4103368#.X43aD0IzZFQ)

To recreate the contents of figures 4 and 5 use the scripts found in the IDW_scripts folder. <br />
To recreate the contents of figures 6-9 run the associated notebook after creating the environment using the following commands<br />
- conda create --name pigBP --file requirements.txt
- conda activate pigBP
To create a Cross Species Browser run the notebook. <br /> 
Prior to running any of the associated notebooks please ensure your environment is configured with the packages isted in the requirements.txt



The contents of each folder are as follows 
* IDW_scripts
    - scripts to perform cross-species spatial prediction with Inverse Distance Weighted Interpolation as demonstrated in figures 4 and 5.
    - For the manuscript regions were defined from the pig's [Saikali atlas](https://www.sciencedirect.com/science/article/abs/pii/S0165027010004164?via%3Dihub) 
    - For the human the Harvard/Oxford atlas was used. ROIs we're generated/saved in standard space in FSLeye's 
    - ROIs were projected to the surface using the fsl command surf_proj before being passed to the scripts included in this folder. 
    - The Calculate KL matlab script comes from the original connectivity blueprint paper and it's associated github repo which can be found [here](https://github.com/elifesciences-publications/rmars-comparing-connectivity-blueprints) 
* PNI50_surfaces_FSICO5
    - The Freesurfer compatible average surfaces in Icosahedron 5 space. 
    - *this is not the space of the final blueprint* but can be used for spherical registration. The figures use the AP corrected surfaces.
* T1_templatePNI50
    - The pig T1 weighted full and brain extracted template as well and segmentation priors. 
* Atlases 
    - The WM and adapted Saikali et al. Atlases. 
        * for use in FSLeyes add this folder as is to $FSLDIR/data/atlases
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
    - Tractograms showing each tracts termination points on the cortical surface. 
    - To visualize load the human and pig tract projection spec files
* Tracts
    - Individual Nifti files of the tracts in the WM atlas
* Surfaces
    - Surfaces used in the generation of the connectivity blueprint. 
    - Human surfaces are from the original Mars et al. [blueprint paper](https://github.com/elifesciences-publications/rmars-comparing-connectivity-blueprints)
    - Native pig surfaces to the PNI50 are in the pig_surfaces/Surf_Volume_geometry folder
    - The pig surfaces called in the notebooks and spec files have been oriented with a corrected anterior-posterior and inferior-superior labelling.
* Cross_species_browser
    - To use the cross species browser you need to run the Create_XspeciesBrowser.ipynb notebook. 
        - Doing so will create the dense connectome files using HCP workbench command and add them to the empty ?.CrossSpeciesBrowser.spec files
    - Alternatively the cifti files can be downloaded from [Zenodo](https://zenodo.org/record/4103368#.X43aD0IzZFQ) 
        - Once Downloaded open them from wb_view after loading the associated empty spec file. 
        - You will be able to click on a vertex in one species and find the vertex in the other species which is the closest in the connecotpic space