////////////////////////////////////////////////////////////////////////////////
// Macro to calculate the rms for the given run and subrun
// TODO loop over subrubs to process
//
//mailto:andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////

//TODO dynamic path to input and output file

//c++ includes
#include <glob.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h> //to use struct stat

//root includes
#include "TFile.h"
#include "TEfficiency.h"

//projects include
#include "Run.hh"
#include "DataStructure.hh"
#include "Cuts.hh"
#include "Utils.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Main macro

int main( int argc, char* argv[] ){

  //parse the input argument ===================================================
  int runNum, subrun;

  for ( int i=1; i<argc; i=i+2 ) {
   if      ( string(argv[i]) == "-run" ) runNum= atoi(argv[i+1]);
   else if ( string(argv[i]) == "-subrun" ) subrun = atoi(argv[i+1]);
   else {
     //PrintUsage();
     return 1;
   }
  }

  string commonPath = "/eos/experiment/wa105/offline/LArSoft/Data/";
  string rawFile = commonPath+"/Raw/ROOT/"+to_string(runNum)+"/"+to_string(runNum)+"-"+to_string(subrun)+"-RawWaveform-Parser.root" ;
  string recoFile =  commonPath+"Reco/2018_June_24/ROOT/recofast/"+to_string(runNum)+"/"+to_string(runNum)+"-"+to_string(subrun)+"-RecoFast-Parser.root";

  //define output ==============================================================

  string outputFilename = "noiseRms-"+to_string(runNum)+".root";
  TFile *ofile = new TFile(outputFilename.c_str(), "RECREATE");
  if(!ofile->IsOpen()){
    cout << "File: " << outputFilename << " cannot be open!" << endl;
    return 1;
  }

  //here I initialize the parser object ========================================

  LArParser *rawParser = new LArParser();
  LArParser *recoParser = new LArParser();

  Run *run = new Run(runNum, "metadata/test.db");

  TTree *rawTree = getTTree( rawFile );
  TTree *recoTree = getTTree( recoFile );

  rawParser->setTTree(rawTree);
  rawParser->setRun(run);

  recoParser->setTTree(recoTree);
  recoParser->setRun(run);

  //check if the tree exists and has been correctly set
  if( !recoParser->isTreeGood() || !recoParser->isTreeGood() ){
   cout << "Invalid ttree" << endl;
   return 1;
  }

  //Event looper ===============================================================

  for(int evt=0; evt<recoTree->GetEntries(); evt++){

      //reco objects
      vector<Hit> recoHits; //hits not associated to a track
      recoParser->getRecoHitsEvent( recoHits, evt );

      //raw objects
      vector<Channel> rawChannels;
      rawParser->getChannelsEvent( rawChannels, evt );

      //loop over channels
      for(auto rawChannel : rawChannels){

        //Skip the channels with reconstructed hits

        //save the waveforms of the channels without hits

        //cout the number of different channels are found across the event

      }

      //calculate mean and rms for each view

      //fill the TTree and save the output

  }//end event loop

}//end macro
