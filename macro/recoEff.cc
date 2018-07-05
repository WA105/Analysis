////////////////////////////////////////////////////////////////////////////////
// Macro studying the reconstruction efficiency of one MC sample
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
#include "TChain.h"
#include "TEfficiency.h"

//projects include
#include "Run.hh"
#include "DataStructure.hh"
#include "Cuts.hh"
#include "Efficiency.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//initalize some variables as global variables (easy to edit)

vector<int> fileList = {}; //runs to process
int startFile = 0;
int endFile = 10;

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

void fillFileList(string path){
  //fill up the runList with all the file if empty

  for(int i=startFile; i<endFile; i++){
    fileList.push_back(i);
  }

  /*
  if (fileList.size() == 0){
    #if verbose
    cout << "Processing all runs in " << path << "*..." << endl;
    #endif
    string wildcard_path = path + "*";
    for( auto irun : glob(wildcard_path) ){
      fileList.push_back(atoi(irun.data()));
    }
  }
  */

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
  //TTree *mcOutputTree = new TTree("mcTracks", "contains geant tracks");
  //mcOutputTree->Branch("mcTrack", &mcTrack );
  //reco tree
  //TTree *recoOutputTree = new TTree("recoTracks", "contains reco tracks");
  //recoOutputTree->Branch("recoTrack", &recoTrack );

  //here I define the parser object
  LArParser *mcParser = new LArParser();
  LArParser *recoParser = new LArParser();

  //define the Run object using the mockRun flag (always the same in this case)
  Run *run = new Run(mockRun, "metadata/test.db");

  //fill the filelist if empty
  fillFileList("/Users/scarpell/cernbox/311/simulation/ana/");

  //and here i define the class efficiency
  Efficiency *recoEfficiency = new Efficiency();

  //loop over all the runs
  for(int fileNumber : fileList){

     TTree *mcTree = getTTree("/Users/scarpell/cernbox/311/simulation/g4detsim/", "", "-G4Detsim-Parser.root", fileNumber ); //TODO dynamic path
     TTree *recoTree = getTTree("/Users/scarpell/cernbox/311/simulation/ana/", "", "-RecoFull-Parser.root", fileNumber ); //TODO: dyynamic path

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
  }//end filelist run


  //calculate the efficency plot
  recoEfficiency->makeEfficiencyPlot();

  ofile->cd();
  recoEfficiency->write();

  ofile->Close();

}//end macro
