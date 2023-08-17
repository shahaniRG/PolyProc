# PolyProc
PolyProc is to to import, process, and analyze 3DXRD datasets of varying dimensions, from 2D to higher dimensions.
Through effective computational routines, the toolbox allows us to align and track features in a highly accurate, efficient, and robust manner.

PolyProc is an open source MATLAB toolbox for the processing and analysis of grain maps across space and time. Its main features are ...

* registration of multiple 3D data sets
* grain map clean up routines (filtering low qualitiy data, removal of small grains, etc)
* statistical analysis - texture(IPF), grain volume histogram, topology(number of neighbors), grain boundary misorientation, compactness(sphericity)
* visualization - 3D volume, 2D slice, grain boundary, multimodal, PMDF
* grain tracking
    
# MATLAB requirements
1. MATLAB version – 2019a or newer version

2. Required MATLAB toolboxes
- Statistical and Machine learning Toolbox
- Global Optimization Toolbox
- Image processing toolbox
- Optimization Toolbox
- Text Analytics Toolbox

# Dependency
PolyProc requires the following dependencies:

1. MTEX is a required Dependency for the code (https://mtex-toolbox.github.io)
(All PolyProc codes are tested compatible with MTEX v. 5.1)

2. Other required Dependencies are already included in `PolyProc/utilities`.

# Getting Started
From MatLab:
Home tab > New > Project > from Git > (clone htpps above)
Download MTEX and save to newly created Project folder, add subfolders to MatLab path

...

# Trial data
Trial data can be downloaded from the following link.

https://bit.ly/2R8SZgC 
    
# Reference Publications
J. Kang, N. Lu, I. Loo, N. Senabulya, and A.J. Shahani, “PolyProc: A Modular Processing Pipeline for X-ray Diffraction Tomography,” Integr. Mater. Manuf. Innov. (2019). https://doi.org/10.1007/s40192-019-00147-2
