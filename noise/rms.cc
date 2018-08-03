////////////////////////////////////////////////////////////////////////////////
// ROOT macro to calculate the chan rms for one event
//
//mailto:andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////

//c++ includes
#include <glob.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h> //to use struct stat

//root includes
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

//projects include
#include "Geometry.hh"
#include "Run.hh"
#include "DataStructure.hh"
#include "Cuts.hh"
#include "Utils.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////////////

void GetMeanAndRMS(vector<double> array, double & mean, double & rms){
  //use welford method..consisted with qScan

  mean = 0;
  rms  = 0;

  if( array.size() == 0 )
    return;

  double A = 0;
  double Q = 0;

  for(size_t i=0;i<array.size();i++){

    double d  = (double)array[i];
    double Ak = A + (d - A)/(i+1);
    double Qk = Q + (d - Ak)*(d-A);
    A = Ak;
    Q = Qk;
  }

  mean = A;
  rms  = sqrt( Q/(array.size()-1) );

  return;
}

////////////////////////////////////////////////////////////////////////////////

int getFileNumber( std::string filename ){
  //assuming the filename encorded in a format /path/to/file/run-subrun-Parser.root

  //isolate filenumber from path
  std::string s = filename;
  std::string delimiter = "/";

  size_t pos = 0;
  std::string token;

  while ((pos = s.find(delimiter)) != std::string::npos) {
    token = s.substr(0, pos);
    s.erase(0, pos + delimiter.length());
  }

  pos = s.find("-");
  std::string number = s.substr(0, pos);
  return stoi(number);
}


////////////////////////////////////////////////////////////////////////////////

void rms(string filename, int evtnumber ){

  //----------------------------------------------------------------------------

    int runNum = getFileNumber( filename );

    string name = "_"+to_string(runNum);
    string title = "Run: " to_string(runNum);

    TH1D *hChMean = new TH1D( ("hChMean"+name).c_str(), ("").c_str(),
                                                      Ch_0+Ch_1, 0, Ch_0+Ch_1 );
    TH1D *hChRMS = new TH1D( ("hChRMS"+name).c_str(), ("").c_str(),
                                                      Ch_0+Ch_1, 0, Ch_0+Ch_1 );

  //----------------------------------------------------------------------------

  LArParser *rawParser = new LArParser();
  TTree *rawTree = getTTree( filename );

  Run *run = new Run(runNum, "metadata/test.db");

  rawParser->setTTree(rawTree);
  rawParser->setRun(run);

  //check if the tree exists and has been correctly set
  if( !rawParser->isTreeGood() ){
     cout << "Invalid ttree" << endl;
     return 1;
  }

  vector<Channel> rawChannels;

  if( evtnumber < rawTree->GetEntries() ){
      rawParser->getRawChannelsEvent( rawChannels, evtnumber );
  }
  else{
      cout << "ERROR! Invalid event number" << endl;
      break;
  }

  for(auto rawChannel : rawChannels){

      double mean, rms;
      GetMeanAndRMS( noiseAdc, mean, rms );

      hChMean->Fill( rawChannel.channel, mean );
      hChRMS->Fill( rawChannel.channel, rms );

  } //end ch loop

  //define output filename and write the histograms in it ----------------------
  string ofilename = "chMeanAndRMS_"+to_string(runNum)+".root";
  TFile *ofile = TFile(ofilename.c_str(), "RECREATE");

  ofile->cd();
  hChMean->Write();
  hChRMS->Write();
  ofile->Close();

}//end macro
