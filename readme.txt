
Folders 'Fig1', 'Fig2', 'Fig3', 'Fig3_S1', 'Fig3_S2', etc.
contain MATLAB scripts and data for generating the corresponding figures.
A figure generator script always has the same name as its folder.
For example, 'Fig3_S2.m' MATLAB script in the folder 'Fig_S2'
generates Figure 3 figure supplement 2 in pdf and png formats.
Both outputs (pdf and png) can be found in the same folder.
Fig. 1 is an exception, this figure is constructed in Keynote,
therefore there is a Keynote file (.key) in folder 'Fig1'.

If you want to regenerate not only a figure but simulations too,
you have to run other scripts in the corresponding folder.
In this case some data files will be overwritten in that folder.
Some common predefined functions can be found in the 'function' folder.
Do not change the folder structure, because functions are loaded from this location.
Functions and scripts (.m files) have descriptions and annotations,
just open them with MATLAB or a simple text viewer.
At the beginning of each script used functions and data files are listed.

Most of the scripts run simulations and have very similar structure,
for understanding purpose script_2d in folder 'Fig3' is mostly recommended to study,
because it has the most detailed annotation.
'_2d' in its name has the meaning that the stimulus bank is 2-dimensional:
next to orientation we have a nuisance parameter, in this case, the phase.
For example 'script_1d' belongs to the no nuisance case,
and script_1d2d3d4d in folder 'Fig4' simulates all nuisance combinations:
stimulus orientation, phase, spatial period, and contrast.

The following syntax for variable names is used in most of the scripts:

N: number of cells (Gabor filters) in the population
K: number of discrete orientation categories
MK: number of stimulus orientation bins within a category (see M_theta in Table 1 in the paper)
J: number of bins for stimulus phase (see M_phi in Table 1 in the paper)
M1: repetition of the stimulus bank for training (see M_rep in Table 1 in the paper)
M2: repetition of the stimulus bank for independent test (see M_rep' in Table 1 in the paper)
disk_radius: aperture radius in degrees (see R in Table 1 in the paper)
V0: mean membrane potential [mV] (see u_DC in Table 1 in the paper)
V1: peak MP amplitude [mV] (see u_AC in Table 1 in the paper)
noise: std of Gaussian membrane potential noise [mV] (see sigma in Table 1 in the paper)
lambda: spatial period [deg] of the stimulus or the Gabor filter
threshold: FRNL threshold [mV]

If you have any questions just email to gasparmerse@gmail.com