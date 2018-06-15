////////////////////////////////////////////////////////////////////////////////
// Select and study different track cuts
//
//mailto:andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////

//NB: the code is a bit slow

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

vector<int> runList = {840}; //runs to process
string path="/Users/scarpell/cernbox/311data/"; //target path to input data ( e.g. ntuples on /eos/ )
string suffix="Ana";

////////////////////////////////////////////////////////////////////////////////

inline bool ExistTest (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

////////////////////////////////////////////////////////////////////////////////

void trackCuts(){
  //main macro filling up the data products and appying some cuts

  //output file
  Track testTrack;
  TFile *ofile = new TFile("cutOutput.root", "RECREATE");
  TTree *tpcTree = new TTree("tpcCut", "contains tracks that cross the tpc");
  tpcTree->Branch("selectedTracks", &testTrack );

  //here I define the parser object
  LArParser *parser = new LArParser();

  //here I define the cuts objects to test
  Cut *testCut = new Cut( {2,2,2,2,2,2}, 100 );

  //start loop over the runs in the runList
  for(int runNumber : runList){ //loop over runsclear

    //define the run object and fill up the data products
    Run *run = new Run(runNumber, "metadata/test.db");

    //start loop over the subRuns
    int subRunCount = run->getNumberOfSubruns(); //keep track of the effective subrun number
    for(int subrun=0; subrun< 4; subrun++){

      string folder=to_string(runNumber)+"/";
      string file=to_string(runNumber)+"-"+to_string(subrun)+"-"+suffix+".root";
      string filename=path+folder+file;

      if( ExistTest(filename) ){

        cout << "Processing file: " << filename << endl;
        TFile *file = new TFile(filename.c_str());

        if( file->IsOpen() ){
          TTree *ttree = (TTree*)file->Get("analysistree/anatree");

          //set the parser with the correct run and tree
          parser->setTTree( ttree );
          parser->setRun( run );

          //check if the tree exists and has been correctly set
          if( parser->isTreeGood() ){

            //now I loop over the events of the tree
            for(int evt=0; evt< ttree->GetEntries(); evt++ ){

              //data structure arrays
              vector<Hit> hits;
              vector<Track> tracks;

              //use the parser to fill up the data structures
              parser->getRecoHitsEvent( hits, evt );
              parser->getRecoTracksEvent( tracks, evt );

              //here I loop over the track elements and save only the tracks meeting the cuts requirements
              for(auto track : tracks){
                //check the cuts and save to file
                testTrack = track;
                if( testCut->isCrossingTPC( testTrack ) ){
                  tpcTree->Fill();
                }
              }

              //open a new tfile and save the objects there (one tree for every cut)
              //define cuts here
            }
          }else{
            cout << "Invalid TTree in file: " << filename << endl;
          }
        }
        file->Close();
      }else{
        cout << "File " << filename << " not found." << endl;
      }

    }//end subrun
  }//end for run

  ofile->cd();
  tpcTree->Write();
  ofile->Close();

}//end macro
