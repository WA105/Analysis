# LArTPC Analysis
Set of useful libraries and examples of macros for the dual-phase LArTPC analysis. The folder provides a list of simple, but useful functions and examples for the data analysis of the 311, which can easly adapted to other detector geometry. For questions:<a href="mailto:andrea.scarpelli@cern.ch" target="_blank"> andrea.scarpelli@cern.ch </a>
## Prerequisites
The project requires ROOT 6.10: <a href="https://root.cern.ch/downloading-root" target="_blank">https://root.cern.ch/downloading-root</a> or a more recent version.
## Organization of the repository
The repository contains the basic classess, functions and common values for the analysis of the 311 charge data. 
..* CommonValues: this repository describes the common cuts and corrections to apply and which are used as default for the analysis
..* Event-track-selection: this repository contains all the routines and instructions necessary to operate a track selection using the [Highway algorithm] (https://github.com/ascarpel/Analysis/blob/master/Event-track-selection/HighwayAlgorithm/README.md). 
..* dqdx: example of a routine to perform a track selection and retrurn a ROOT file with relevant event quantities
..* efficiency: example of a routine to assess the reconstruction efficiency using the Monte Carlo sample
..* noise: example of routines for the study of the noise
..* macro: location of the ROOT macros
..* source: ocation of the custom functions necessary for the analysis
..* include: location of the headers of the files in source
..* database: .csv file with all the runs metadata
## Getting Started
### Run as ROOT macro
Analysis can be performed using standard ROOT macro. The custom functions defined in include and source are loaded into ROOT with the routine loadLib.cc which must be called before your macro as follows:
```
root -l loadLib.cc myMacro.cc
```
or
```
root -l loadLib.cc
root[0] .x myMacro.cc
```
Do not forget to add the appropriate includes in your macro. dQds.cc is an example macro to plit the dQ/ds disribution of a given ROOT file and fit with a Langaus distribution.
### Run as compiled executable
It is possible to complile a routine using CMake to link all the necessary libraries. The scripts in dqdx are an example of this. The paths configured in CMakeLists.txt should be modified by the user. First, a build repository is created: 
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
