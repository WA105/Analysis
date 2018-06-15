#ifndef  __311STYLE_H
#define __311STYLE_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile2D.h>
#include <TF1.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>

//project libraries
#include <311Lib.h>

using namespace std;

Color_t getColor(int view);

#endif // __311STYLE_H
