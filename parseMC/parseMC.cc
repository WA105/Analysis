////////////////////////////////////////////////////////////////////////////////
// Macro to parse the MC ROOT files into an easier format using Traack and Hit
// data classes
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
<<<<<<< HEAD:utils/parseMC.cc
#include "TEfficiency.h"
=======
>>>>>>> cmake:parseMC/parseMC.cc

//projects include
#include "Run.hh"
#include "DataStructure.hh"
#include "Efficiency.hh"

////////////////////////////////////////////////////////////////////////////////
// Functions

inline bool ExistTest (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

TTree *getTTree( string path, string prefix, string suffix, int runNumber ){

  TTree *ttree;

  //no need to run over subruns (are mc). Just jump to filename definition
  string file=prefix+to_string(runNumber)+suffix;
  string filename=path+file;

  if( ExistTest(filename) ){

      cout << "Reading file: " << filename << endl;
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
     int fileNumber;
     std::string mcFolder;
     int mockRun;
     for ( int i=1; i<argc; i=i+2 ) {
       if      ( string(argv[i]) == "-n" ) fileNumber = atoi(argv[i+1]);
       else if ( string(argv[i]) == "-f" ) mcFolder = argv[i+1];
       else if ( string(argv[i]) == "-r" ) mockRun = atoi(argv[i+1]);
       else {
         //PrintUsage();
         return 1;
       }
     }

     //define here the output
     Hit recoHit;
     Track recoTrack;
     MCTrack mcTrack;
     //define output file
     string outputFilename = "/eos/user/a/ascarpel/311/simulation/"+mcFolder+"/Simulation_"+to_string(mockRun)+"_"+to_string(fileNumber)+".root";
     TFile *ofile = new TFile(outputFilename.c_str(), "RECREATE");
     if(!ofile->IsOpen()){
       cout << "File: " << outputFilename << " cannot be open!" << endl;
       return 1;
     }
     //mc tree
     TTree *mcOutputTree = new TTree("mcTracks", "contains geant tracks");
     mcOutputTree->Branch("mcTrack", &mcTrack );
     //reco hit tree
     TTree *recoHitOutputTree = new TTree("recoHits", "contains reco hits");
     recoHitOutputTree->Branch("recoHit", &recoHit );
     //reco track tree
     TTree *recoOutputTree = new TTree("recoTracks", "contains reco tracks");
     recoOutputTree->Branch("recoTrack", &recoTrack );

     //here I define the parser object
     LArParser *mcParser = new LArParser();
     LArParser *recoParser = new LArParser();

     //define the Run object using the mockRun flag (always the same in this case)
     Run *run = new Run(mockRun, "metadata/test.db");

  mcParser->setTTree(mcTree);
  mcParser->setRun(run);

  recoParser->setTTree(recoTree);
  recoParser->setRun(run);

  //loop over the events in mcTree. Should be the same for also recoTree
  for(int evt=0; evt<mcTree->GetEntries(); evt++){
     //check if the tree exists and has been correctly set
     if( !mcParser->isTreeGood() || !recoParser->isTreeGood() )
     {
       cout << "Invalid ttree" << endl;
       return 1;
     }

    //check if the two trees have the same number of entries. I ideally want to loop over one of them
    if( mcTree->GetEntries() != recoTree->GetEntries() )
    {
      cout << "Not the same number of events" << endl;
      return 1;
    }

    //loop over the events in mcTree. Should be the same for also recoTree
    for(int evt=0; evt<mcTree->GetEntries(); evt++)
    {

      //data structures array
      vector<MCTrack> mcTracks;
      vector<Track> recoTracks;
      vector<Hit> recoHits; //hits not associated to a track

      mcParser->getMCTracksEvent( mcTracks, evt);
      recoParser->getRecoTracksEvent(recoTracks, evt);
      recoParser->getRecoHitsEvent( recoHits, evt );

      for( auto track : mcTracks ){
        mcTrack = track;
        mcOutputTree->Fill();
      }

      for( auto track : recoTracks ){
        recoTrack = track;
        recoOutputTree->Fill();
      }

      for(auto hit : recoHits){
        recoHit = hit;
        recoHitOutputTree->Fill();
      }
    }//end event loop

  ofile->cd();

  mcOutputTree->Write();
  recoOutputTree->Write();
  recoHitOutputTree->Write();

  ofile->Close();

 return 0;
}
