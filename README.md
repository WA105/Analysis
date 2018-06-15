Created by Andrea Scarpelli (andrea.scarpelli@cern.ch) on dec 5th

# 311Analysis

Add a description here

# Prerequisite

The project requires ROOT 6.10: <a href="https://root.cern.ch/downloading-root" target="_blank">`https://root.cern.ch/downloading-root`</a>

# Getting Started

> Download the folder

> Run a ROOT macro

To run a macro:
```
root -l loadLib.cc myMacro.cc
```
or
```
root -l loadLib.cc
root[0] .x myMacro.cc
```

> How to write a macro

# Folder organization

# Contribute

> Add your Library

> Suggest a modification


WARNING: This repository was at first created for author's use only. Please
report every error, comment or question to andrea.scarpelli@cern.ch

Set of root macro for the analysis of the 311 data. Macros are intended to be used
after loading the libraries inside lib folders containing useful function for the
analysis.

Some of the macro requires simple root and no previous knowledge of LArSoft, some
other requires the use of the ups package "Gallery" coming with LArSoft packages
How to use Gallery:

loadLib.cc script will load and compiled the libraries listed in a new folder called obj
Please do not push obj folder.

to run root loading all the libraries: root -l loadLib.cc to be rerun every time
root is launched. If libraries contents is not touched, the compilation won't be
done a second time.

Organization of data products:
  Some function requires as input track and hits classes
  the function read_tree will parse the input root file into that class, filling
  the relative important quantities. Use is rather simple and macros already in
  the repository will provide nice example on how to use them.

The others macro are divided into folders specific for the type of analysis for
which they are conceived:

utils: miscellaneous macro for different purposes
  1.) track_cut.cc read input from C. Alt AnaRootParserFiles and selects mips.
    Everything is stored in a new root file in a branch called with the type of
    cut performed.
  2.) evd_311.cc and evd_raw_311.C are gallery macro to display a single event,
  the latter with raw data (with pedestal) and the latter after data preparation
  (eg pedestal and noise removal, depending on the configuration files used)
  3.) ViewToDAQChan.C convert view channels into daq channels


Purity: macro purity.cc for LAr purity estimation
