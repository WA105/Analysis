////////////////////////////////////////////////////////////////////////////////
// Utilities functions and classes
//
//mailto:andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////

#ifndef  __UTILS_H
#define __UTILS_H

#include <glob.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h> //to use struct stat

#include "TFile.h"
#include "TTree.h"

using namespace std;

// Spare functions /////////////////////////////////////////////////////////////

inline bool ExistTest (const std::string& name);

TTree *getTTree( string filename );

#endif // __UTILS_H
