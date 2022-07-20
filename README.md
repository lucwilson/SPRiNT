[![DOI](https://zenodo.org/badge/449755672.svg)](https://zenodo.org/badge/latestdoi/449755672)

# SPRiNT: Time-resolved parameterization of aperiodic and periodic brain activity  

Code for simulating dynamic neural time series, spectrogram parameterization (using SPRiNT), and generating figures from the SPRiNT preprint.

**SPRiNT** is a MATLAB toolbox. This package has an associated [article](https://www.biorxiv.org/content/10.1101/2022.01.21.477243v1.abstract) where we validate the method on both simulated and empirical data.

# Information

Code for runnning the SPRiNT algorithm can be found in SPRiNT.m. The SPRiNT algorithm is also available from the [Brainstorm distribution](https://neuroimage.usc.edu/brainstorm/Introduction) (Tadel et al., 2011). See the [Brainstorm website](https://neuroimage.usc.edu/brainstorm/Tutorials/SPRiNT?highlight=%28SPRiNT%29) for a tutorial on how to run SPRiNT.

Code for simulating time series (see [Study 1](https://osf.io/c3gn4/)) can be found in timeseries_sim.m. Code for generating the figure for Study 1 & 3 can be found in the Figures folder. Code for the EEG dataset (see [Study 2](https://www.nature.com/articles/sdata2018308)) can be found in the LEMON_code folder.

Please contact the authors of this package for any questions and suggestions.

# Citation

If you use this software, please cite it.

To cite SPRiNT in publications, please use:

Wilson, L. E., da Silva Castanheira, J., & Baillet, S. (2022). Time-resolved parameterization of aperiodic and periodic brain activity. bioRxiv, 2022.2001.2021.477243. 
