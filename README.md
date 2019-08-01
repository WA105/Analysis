# LArTPC Analysis
Set of useful libraries and examples of macros for the dual-phase LArTPC analysis. The folder provides a list of simple, but useful functions and examples for the data analysis of the 311, which can easly adapted to other detector geometry. For questions:<a href="mailto:andrea.scarpelli@cern.ch" target="_blank"> andrea.scarpelli@cern.ch </a>
## Prerequisites
The project requires ROOT 6.10: <a href="https://root.cern.ch/downloading-root" target="_blank">https://root.cern.ch/downloading-root</a> or a more recent version.
## Organization of the repository
The repository contains the basic classess, functions and common values for the analysis of the 311 charge data. 
 * _CommonValues_: this repository describes the common cuts and corrections to apply and which are used as default for the analysis
 * _Event-track-selection_: this repository contains all the routines and instructions necessary to operate a track selection using the [Highway algorithm] (https://github.com/ascarpel/Analysis/blob/master/Event-track-selection/HighwayAlgorithm/README.md). 
 * _dqdx_: example of a routine to perform a track selection and retrurn a ROOT file with relevant event quantities
 * _efficiency_: example of a routine to assess the reconstruction efficiency using the Monte Carlo sample
 * _noise_: example of routines for the study of the noise
 * _macro_: location of the ROOT macros
 * _source_: ocation of the custom functions necessary for the analysis
 * _include_: location of the headers of the files in source
 * _database_: .csv file with all the runs metadata
## Getting Started
### Run as ROOT macro
Analysis can be performed using standard ROOT macro. The custom functions defined in include and source are loaded into ROOT with the routine _loadLib.cc_ which must be called before your macro as follows:
```
root -l loadLib.cc myMacro.cc
```
or
```
root -l loadLib.cc
root[0] .x myMacro.cc
```
Do not forget to add the appropriate includes in your macro. _dQds.cc_ is an example macro to plit the dQ/ds disribution of a given ROOT file and fit with a Langaus distribution.
### Run as compiled executable
It is possible to complile a routine using CMake to link all the necessary libraries. The scripts in folder _dqdx_ are an example of this. The paths configured in _CMakeLists.txt_ should be modified by the user. First, a build repository is created: 
```
cd dqdx
mkdir build
cd build
```
Then, It is possible to generate the CMake files and build the project 
```
cmake ../
cmake --build .
```
The exectutable is now located in the build directory. 
