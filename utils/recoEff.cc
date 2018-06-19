////////////////////////////////////////////////////////////////////////////////
// Example of a macro studying the reconstruction efficiency of one sample
//
//mailto:andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////

#include <glob.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h> //to use struct stat

#include "TChain.h"
#include "Run.hh"
#include "DataStructure.hh"
#include "Cuts.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//initalize some variables as global variables (easy to edit)

vector<int> fileList = {0}; //runs to process
int mockRun = 840; //query the metadata of this run from db

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

void fillFileList( ){
  //fill up the runList with all the file if empty
  if (fileList.size() == 0){
    #if verbose
    cout << "Processing all runs in " << path << "*..." << endl;
    #endif
    string wildcard_path = path + "*";
    for( auto irun : glob(wildcard_path) ){
      fileList.push_back(atoi(irun.data()));
    }
  }
  return;
}

TTree *getTTree( string path, string prefix, string suffix, int runNumber ){

  TTree *ttree;

  //no need to run over subruns (are mc). Just jump to filename definition
  string file=prefix+to_string(runNumber)+suffix;
  string filename=path+file;

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

void recoEff(){
  //Example of a macro to study the reconsturction eff.

  //define here the output file
  Track recoTrack;
  MCTrack mcTrack;

  TFile *ofile = new TFile("simOutput.root", "RECREATE");

  //mc tree
  TTree *mcOutputTree = new TTree("mcTracks", "contains geant tracks");
  mcOutputTree->Branch("mcTrack", &mcTrack );
  //reco tree
  TTree *recoOutputTree = new TTree("recoTracks", "contains reco tracks");
  recoOutputTree->Branch("recoTrack", &recoTrack );


  //here I define the parser object
  LArParser *mcParser = new LArParser();
  LArParser *recoParser = new LArParser();

  //define the Run object using the mockRun flag (always the same in this case)
  Run *run = new Run(mockRun, "metadata/test.db");


  //fill the filelist if empty
  fileList();

  //loop over all the runs
  for(int fileNumber : fileList){

     TTree *mcTree = getTTree("/Users/scarpell/cernbox/311/simulation/g4detsim/", "", "-G4Detsim-Parser.root", fileNumber );
     TTree *recoTree = getTTree("/Users/scarpell/cernbox/311/simulation/ana/", "MC5_", "_Ana.root", fileNumber );

     mcParser->setTTree(mcTree);
     mcParser->setRun(run);

     recoParser->setTTree(recoTree);
     recoParser->setRun(run);

     //check if the tree exists and has been correctly set
     if( !mcParser->isTreeGood() || !recoParser->isTreeGood() ){
       cout << "Invalid ttree" << endl;
       continue;
     }

    //check if the two trees have the same number of entries. I ideally want to loop over one of them
    if( mcTree->GetEntries() != recoTree->GetEntries() ){
      cout << "Not the same number of events" << endl;
      continue;
    }

    //loop over the events in mcTree. Should be the same for also recoTree
    for(int evt=0; evt<mcTree->GetEntries(); evt++){

      //data structures array
      vector<MCTrack> mcTracks;
      vector<Track> recoTracks;

      mcParser->getMCTracksEvent(mcTracks, evt);
      recoParser->getRecoTracksEvent(recoTracks, evt);

      for( auto track : mcTracks ){
        mcTrack = track;
        mcOutputTree->Fill();
      }

      for( auto track : recoTracks ){
        recoTrack = track;
        recoOutputTree->Fill();
      }

    }//end event loop
  }//end filelist run


  ofile->cd();
  mcOutputTree->Write();
  recoOutputTree->Write();
  ofile->Close();

}//end macro
