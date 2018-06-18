////////////////////////////////////////////////////////////////////////////////
// Study reconstruction efficiency in the MCSample
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

vector<int> runList = {0}; //runs to process
string path="/Users/scarpell/cernbox/311/simulation/g4detsim/"; //target path to input data ( e.g. ntuples on /eos/ )
string prefix="";
string suffix="-G4Detsim-Parser.root";

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

void fillRunList( ){
  //fill up the runList with all the file if empty
  if (runList.size() == 0){
    #if verbose
    cout << "Processing all runs in " << path << "*..." << endl;
    #endif
    string wildcard_path = path + "*";
    for( auto irun : glob(wildcard_path) ){
      runList.push_back(atoi(irun.data()));
    }
  }
  return;
}

////////////////////////////////////////////////////////////////////////////////
// Main macro

void recoEff(){

  fillRunList();

  //output file
  Track testTrack;
  //TFile *ofile = new TFile("cutOutput.root", "RECREATE");
  //TTree *tpcTree = new TTree("tpcCut", "contains tracks that cross the tpc");
  //tpcTree->Branch("selectedTracks", &testTrack );

  //here I define the parser object
  LArParser *parser = new LArParser();

  //start loop over the runs in the runList
  for(int runNumber : runList){ //loop over runList

    //define the run object and fill up the data products
    Run *run = new Run(mockRun, "metadata/test.db");

    //no need to run over subruns (are mc). Just jump to filename definition
    string file=prefix+to_string(runNumber)+suffix;
    string filename=path+file;

    if( ExistTest(filename) ){

        cout << "Processing file: " << filename << endl;
        TFile *file = new TFile(filename.c_str());

        if( file->IsOpen() ){
          TTree *ttree = (TTree*)file->Get("analysistree/anatree");

          //set the parser with the correct run and tree
          parser->setTTree(ttree);
          parser->setRun(run);

          //check if the tree exists and has been correctly set
          if( parser->isTreeGood() ){

            //now I loop over the events of the tree
            for(int evt=0; evt< ttree->GetEntries(); evt++ ){

              //data structure arrays
              vector<MCTrack> mctracks;

              //use the parser to fill up the data structures
              parser->getMCTracksEvent(mctracks, evt);

              for( auto mctrack : mctracks ){
                cout << mctrack.startE << endl;
              }

            }
          }else{
            cout << "Invalid TTree in file: " << filename << endl;
          }
        }
        file->Close();
    }else{
        cout << "File " << filename << " not found." << endl;
    }

  }//end for run

  //ofile->cd();
  //tpcTree->Write();
  //ofile->Close();

}//end macro
