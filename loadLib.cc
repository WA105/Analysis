#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

//NB: Library cuts not yet completed

void loadLib(){

  char cwd[1024];
  getcwd(cwd, sizeof(cwd));

  gStyle->SetPalette(kRainBow);
  gStyle->SetOptFit(11111);

  std::string currentdir( cwd );
  std::string includepath="-I"+currentdir+"/include/";

  gSystem->SetBuildDir("obj",true);
  gSystem->AddIncludePath(includepath.c_str());
  gROOT->LoadMacro((currentdir+"/source/Geometry.cc+").c_str());
  gROOT->LoadMacro((currentdir+"/source/Run.cc+").c_str());
  gROOT->LoadMacro((currentdir+"/source/DataStructure.cc+").c_str());
  //gROOT->LoadMacro((currentdir+"/source/Cuts.cc+").c_str());

  //gROOT->LoadMacro((currentdir+"/source/311Lib.cc+").c_str());
  //gROOT->LoadMacro((currentdir+"/source/311style.cc+").c_str());

  #define __INITIALIZED__

}//end loadLib
