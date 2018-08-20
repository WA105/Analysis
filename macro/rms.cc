////////////////////////////////////////////////////////////////////////////////
// ROOT macro to calculate the chan rms for one fullreco file
//
// Usage: root -l loadLib.cc rms.cc
// [0] .x rms("path/to/fullreco.file")
//
// mailto:andrea.scarpelli@cern.ch
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

void GetMeanAndRMS(vector<float> array, float & mean, float & rms){
  //use welford method..consisted with qScan

  mean = 0;
  rms  = 0;

  if( array.size() == 0 )
    return;

  float A = 0;
  float Q = 0;

  for(size_t i=2;i<array.size();i++){

    float d  = (float)array[i];
    float Ak = A + (d - A)/(i+1);
    float Qk = Q + (d - Ak)*(d-A);
    A = Ak;
    Q = Qk;
  }

  mean = A;
  rms  = sqrt( Q/(array.size()-1) );

  return;
}

////////////////////////////////////////////////////////////////////////////////

float getMeanVector( vector<float> v ){
  //quick evaluate the mean of one vector
  float mean = 0;

  if( v.size() == 0 ){
    return mean;
  }
  else{
    float sum = accumulate(std::begin(v), std::end(v), 0.0);
    mean = sum/(float)v.size();
    return mean;
  }
}

////////////////////////////////////////////////////////////////////////////////

int getFileNumber( std::string filename, subtra ){
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

void rms(string filename, string ofilename){

  //----------------------------------------------------------------------------

    int runNum = getFileNumber( filename );

    string name = "_"+to_string(runNum);
    string title = "Run: " + to_string(runNum);

    TH1D *hChMean = new TH1D( ("hChMean"+name).c_str(), (""+title).c_str(),
                                                      Ch_0+Ch_1, 0, Ch_0+Ch_1 );
    TH1D *hChRMS = new TH1D( ("hChRMS"+name).c_str(), (""+title).c_str(),
                                                      Ch_0+Ch_1, 0, Ch_0+Ch_1 );

  //----------------------------------------------------------------------------

  LArParser *rawParser = new LArParser();
  TTree *rawTree = getTTree( filename );

  map<int, vector<float>> fCh2Mean;
  map<int, vector<float>> fCh2RMS;

  int mod = 6;
  if( rawTree->GetEntries() < 300 ){
    mod = rawTree->GetEntries()/50; //it will process event by event only max 100 events
  }

  for(int evt=0; evt<rawTree->GetEntries(); evt++){

      if( evt % mod != 0 ){ continue; } //process one event every 6 ( about 50 events per subrun w 335 events )
      cout << " Processing event " << evt << endl;

      vector<Channel> rawChannels;
      rawParser->getRecoChannelsEvent( rawTree, rawChannels, evt );

      for(auto rawChannel : rawChannels){

        float mean, rms;

        rawChannel.subtractPedestal(true); //allow pedestal subtraction for that channel

        GetMeanAndRMS( rawChannel.signal, mean, rms );

        fCh2Mean[rawChannel.channel].push_back(mean);
        fCh2RMS[rawChannel.channel].push_back(rms);

      } //end ch loop


    map<int, vector<float>>::iterator mapIt;
    for( mapIt = fCh2Mean.begin(); mapIt != fCh2Mean.end(); ++mapIt ){

      float mean = getMeanVector( fCh2Mean[ mapIt->first ] );
      float rms = getMeanVector( fCh2RMS[ mapIt->first ] );

      int channel = ViewToDAQChan( mapIt->first );

      hChMean->SetBinContent( channel, mean );
      hChRMS->SetBinContent( channel, rms );

    }

  }//end event loop

  //define output filename and write the histograms in it ----------------------
  //string ofilename = "chMeanAndRMS_"+to_string(runNum)+".root";
  TFile *ofile = new TFile(ofilename.c_str(), "RECREATE");

  ofile->cd();
  hChMean->Write();
  hChRMS->Write();
  ofile->Close();

}//end macro
