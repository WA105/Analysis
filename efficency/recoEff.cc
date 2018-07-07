////////////////////////////////////////////////////////////////////////////////
// Macro studying the reconstruction efficiency of one MC sample
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
#include "Efficiency.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Functions

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

TTree *getTTree( string filename ){

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
// Main macro

int main( int argc, char* argv[] ){

  //parse the input argument
  std::string simFile;
  std::string recoFile;
  int mockRun;

  for ( int i=1; i<argc; i=i+2 ) {
   if      ( string(argv[i]) == "-s" ) simFile = argv[i+1];
   else if ( string(argv[i]) == "-r" ) recoFile = argv[i+1];
   else if ( string(argv[i]) == "-n" ) mockRun = atoi(argv[i+1]);
   else {
     //PrintUsage();
     return 1;
   }
  }

  //define here the output file
  Track recoTrack;
  MCTrack mcTrack;

  //define output file
  string outputFilename = "recoEfficiency.root";
  TFile *ofile = new TFile(outputFilename.c_str(), "RECREATE");
  if(!ofile->IsOpen()){
    cout << "File: " << outputFilename << " cannot be open!" << endl;
    return 1;
  }

  //here I define the parser object
  LArParser *mcParser = new LArParser();
  LArParser *recoParser = new LArParser();

  //define the Run object using the mockRun flag (always the same in this case)
  Run *run = new Run(mockRun, "metadata/test.db");

  //and here i define the class efficiency
  Efficiency *recoEfficiency = new Efficiency();

  TTree *mcTree = getTTree( simFile );
  TTree *recoTree = getTTree( recoFile );

  mcParser->setTTree(mcTree);
  mcParser->setRun(run);

  recoParser->setTTree(recoTree);
  recoParser->setRun(run);

  //check if the tree exists and has been correctly set
  if( !mcParser->isTreeGood() || !recoParser->isTreeGood() ){
   cout << "Invalid ttree" << endl;
   return 1;
  }

  //check if the two trees have the same number of entries. I ideally want to loop over one of them
  if( mcTree->GetEntries() != recoTree->GetEntries() ){
    cout << "Not the same number of events" << endl;
    return 1;
  }

  //loop over the events in mcTree. Should be the same for also recoTree
  for(int evt=0; evt<mcTree->GetEntries(); evt++){

      //data structures array
      vector<MCTrack> mcTracks;
      vector<Track> recoTracks;
      vector<Hit> recoHits; //hits not associated to a track

      mcParser->getMCTracksEvent(mcTracks, evt);
      recoParser->getRecoTracksEvent(recoTracks, evt);
      recoParser->getRecoHitsEvent( recoHits, evt );

      for( auto track : mcTracks ){
        //order truth track into a map sorted by their id so is easier to make
        //a match between true and reco
        recoEfficiency->setMapEntry( track.particleID, track );
      }

      //insert the reconstruced hits inside the event
      recoEfficiency->setRecoHits( recoHits );

      for( auto track : recoTracks ){
        //here we do the real efficiency analysis trying to match reconstructed
        //calculate efficiency for the track and fill the historams
        recoEfficiency->setRecoTrack( track );
        recoEfficiency->fill();
      }

      recoEfficiency->clean();

  }//end event loop

  //calculate the efficency plot

  ofile->cd();
  recoEfficiency->write();

  ofile->Close();

}//end macro
