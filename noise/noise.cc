////////////////////////////////////////////////////////////////////////////////
// Macro to calculate the rms for the given run and subrun
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

bool hasChannel( vector<Hit> hits ){
  //check if the vector exists

  if( !hits.size() ){ return false; }
  else { return true; }

}

bool inROI( double t, vector<Hit> hits ){
  //check if every time deposit is in a ROI

  int padSize = 10; //ticks to ignore before and after each hit

  //now check if the time falls in the ROI identified by one hits
  vector<Hit>::iterator pos = find_if( hits.begin(), hits.end(),
  [t, padSize]( Hit & h )->bool
    {
      return ( (h.startTime-padSize <= t) && (t <= h.endTime+padSize) );
    }
  );

  if ( pos == hits.end() ) { return false; }
  else { return true; }

}

double getMeanVector( vector<double> v ){
  //quick evaluate the mean of one vector
  double mean = 0;

  if( v.size() == 0 ){
    return mean;
  }
  else{
    double sum = accumulate(std::begin(v), std::end(v), 0.0);
    mean = sum/(double)v.size();
    return mean;
  }
}

////////////////////////////////////////////////////////////////////////////////
// Main macro

int main( int argc, char* argv[] ){

  //parse the input argument ===================================================

  int runNum; // run to process
  int subrun; // subrun to process
  string filter; //incoherent or coherent

  for ( int i=1; i<argc; i=i+2 ) {
   if      ( string(argv[i]) == "-run" ) runNum= atoi(argv[i+1]);
   else if ( string(argv[i]) == "-subrun" ) subrun = atoi(argv[i+1]);
   else if ( string(argv[i]) == "-filter" ) filter = argv[i+1];
   else {
     //PrintUsage();
     return 1;
   }
  }

  bool filterOn = false;
  if( filter == "yes" ){
    cout << " calculate noise after coherent noise removal " << endl;
    filterOn = true;
  }

  //define and variables =======================================================

  int fRun;
  int fSubrun;
  int fEvent;
  int fEventSeconds;
  int fEventNanoSeconds;
  double fMeanView1, fMeanView0;
  double fRMSView0, fRMSView1;

  string filterName = "";
  if( filterOn )
    filterName = "coherent";

  string outputFilename = "/eos/user/a/ascarpel/Noise/rms/"+filterName+"NoiseRms-"+to_string(runNum)+"-"+to_string(subrun)+".root";
  TFile *ofile = new TFile(outputFilename.c_str(), "RECREATE");
  if(!ofile->IsOpen()){
    cout << "File: " << outputFilename << " cannot be open!" << endl;
    return 1;
  }

  TTree *noiseTree = new TTree("noise", "noise mean and rms");
  noiseTree->Branch("Run", &fRun, "Run/I");
  noiseTree->Branch("Subrun", &fSubrun, "Subrun/I");
  noiseTree->Branch("Event", &fEvent, "Event/I");
  noiseTree->Branch("EventSeconds", &fEventSeconds, "EventSeconds/I");
  noiseTree->Branch("EventNanoSeconds", &fEventNanoSeconds, "EventNanoSeconds/I");
  noiseTree->Branch("MeanView0", &fMeanView0, "MeanView0/D");
  noiseTree->Branch("RMSView0", &fRMSView0, "RMSView0/D");
  noiseTree->Branch("MeanView1", &fMeanView1, "MeanView1/D");
  noiseTree->Branch("RMSView1", &fRMSView1, "RMSView1/D");

  map<int, vector<double> > view2mean;
  map<int, vector<double> > view2rms;

  //Prepare inputs =============================================================

  string commonPath = "/eos/experiment/wa105/offline/LArSoft/Data/";
  string recoFile =  commonPath+"Reco/2018_June_24/ROOT/recofull/"+to_string(runNum)+"/"+to_string(runNum)+"-"+to_string(subrun)+"-RecoFull-Parser.root";
  string rawFile =  commonPath+"Raw/ROOT/"+to_string(runNum)+"/"+to_string(runNum)+"-"+to_string(subrun)+"-RawWaveform-Parser.root" ;

  LArParser *rawParser = new LArParser();
  LArParser *recoParser = new LArParser();

  TTree *rawTree = getTTree( rawFile );
  TTree *recoTree = getTTree( recoFile );

  //check if the tree has been correctly set ========================

  if( recoTree->GetEntries() != rawTree->GetEntries() ){
    cout << "Error! Trees havent the same number of entries! " << endl;
    return 1;
  }

//Event looper ===============================================================

cout << "Start event looper" << endl;

int mod = 6;
if( rawTree->GetEntries() < 300 ){
  mod = rawTree->GetEntries()/50; //it will process event by event only max 100 events
}

  for(int evt=0; evt<rawTree->GetEntries(); evt++){

    if( evt % mod !=0 ){ continue; } //process one event every 6 ( about 50 events per subrun w 335 events )
      cout << " Processing event " << evt;

      //reco objects ----------------------------------------------------------
      vector<Hit> recoHits;
      recoParser->getRecoHitsEvent( recoTree,  recoHits, evt );

      //make an hit channels map
      map<int, vector<Hit> > ch2hits;
      for(auto recoHit : recoHits)
        ch2hits[ recoHit.channel ].push_back( recoHit );

      //raw objects
      vector<Channel> rawChannels;
      rawParser->getRawChannelsEvent( rawTree, rawChannels, evt );

      if(!rawChannels.size()){
        cout << "ERROR:No rawChannels object created. Break event loop" << endl;
        break;
      }

      cout << "..";

      //loop over channels -----------------------------------------------------
      for(auto rawChannel : rawChannels){

        rawChannel.subtractPedestal(true);

        if( !hasChannel( ch2hits[ rawChannel.channel ] ) ){ continue; }
        if ( rawChannel.isDead() ){ continue; }

        vector<double> adc = rawChannel.signal; //real ADCs
        vector<double> noiseAdc; // only noisy ADCs

        for( int t=2; t<(int)adc.size(); t++  ){

          if ( !inROI( t, ch2hits[ rawChannel.channel ] ) ){
            noiseAdc.push_back( adc.at(t) );
          }
        } //end adc loop

        //calculate Mean and RMS here --------------------------------------------

        double mean, rms;
        GetMeanAndRMS( noiseAdc, mean, rms );

        if(rms > 0){
          view2mean[ rawChannel.view ].push_back( mean );
          view2rms[ rawChannel.view ].push_back( rms );
        }

        noiseAdc.clear();

      } //end ch loop

      cout << "..";

      fRun = rawChannels.at(1).run;
      fSubrun = rawChannels.at(1).subRun;
      fEvent = rawChannels.at(1).event;
      fEventSeconds = rawChannels.at(1).timeSeconds;
      fEventNanoSeconds = rawChannels.at(1).timeNanoSeconds;

      //calulate the mean mean and the meam rms
      fMeanView0 = getMeanVector( view2mean[0] );
      fRMSView0 = getMeanVector( view2rms[0] );
      fMeanView1 = getMeanVector( view2mean[1] );
      fRMSView1 = getMeanVector( view2rms[1] );

      cout << "..";

      //fill the tree
      noiseTree->Fill();

      view2mean[0].clear();
      view2rms[0].clear();
      view2mean[1].clear();
      view2rms[1].clear();

      cout << ".." << endl;
  }//end event loop

  //fill and close file ========================================================
  ofile->cd();

  noiseTree->Write();

  ofile->Close();

  cout << "All done! " << endl;

  return 0;

}//end main
