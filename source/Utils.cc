////////////////////////////////////////////////////////////////////////////////
// Utilities functions and classes method implementation
//
//mailto:andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////

#include <glob.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h> //to use struct stat

#include "TFile.h"
#include "Utils.hh"

using namespace std;

// Spare functions /////////////////////////////////////////////////////////////

inline bool ExistTest (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

inline vector<string> glob(const string& pat){
  glob_t glob_result;
  glob(pat.c_str(),GLOB_TILDE,NULL,&glob_result);
  vector<string> ret;
  for(unsigned int i=0;i<glob_result.gl_pathc;++i){
      string full_path = string(glob_result.gl_pathv[i]);
      string basename = full_path.substr(full_path.find_last_of("/")+1);
      ret.push_back(basename.data());
  }
  globfree(&glob_result);
  return ret;
}

////////////////////////////////////////////////////////////////////////////////

TTree *getTTree( string filename ){
  //find the ttree associated to that specific filename: assuming fixed tree name

  TTree *ttree;

  if( ExistTest(filename) ){

      cout << "Processing file: " << filename << endl;
      TFile *file = new TFile(filename.c_str(), "READ");

      if( file->IsOpen() )
        ttree = (TTree*)file->Get("analysistree/anatree");
      else
        cout << "getTTree::Error: Invalid TTree name" << endl;

    }else{
      cout << "getTTree::Error: Invalid file name " << filename << endl;
    }

    return ttree;
}

////////////////////////////////////////////////////////////////////////////////
