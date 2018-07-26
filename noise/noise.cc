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
#include "TTree.h"
#include "TH1.h"

//projects include
#include "Run.hh"
#include "DataStructure.hh"
#include "Cuts.hh"
#include "Utils.hh"

using namespace std;


void GetMeanAndRMS(vector<double> array, double & mean, double & rms){

  //use welford method..consisted with qScan

  mean = 0;
  rms  = 0;

  if( array.size() == 0 )
    return;

  double A = 0;
  double Q = 0;

  size_t istart = 3; //exclude first 3 point, can give crazy results

  for(size_t i=istart;i<array.size();i++){

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
// Main macro

int main( int argc, char* argv[] ){

  //parse the input argument ===================================================
  int runNum, nSubruns;

  for ( int i=1; i<argc; i=i+2 ) {
   if      ( string(argv[i]) == "-run" ) runNum= atoi(argv[i+1]);
   else if ( string(argv[i]) == "-subrun" ) nSubruns = atoi(argv[i+1]);
   else {
     //PrintUsage();
     return 1;
   }
  }

  //define output ==============================================================

  int fRun;
  int fSubrun;
  int fEvent;
  int fEventSeconds;
  int fEventNanoSeconds;
  double fMeanView1, fMeanView0;
  double fRMSView0, fRMSView1;
  int fChannelsCoverage;

  string outputFilename = "/eos/user/a/ascarpel/Noise/rms/noiseRms-"+to_string(runNum)+".root";
  TFile *ofile = new TFile(outputFilename.c_str(), "RECREATE");
  if(!ofile->IsOpen()){
    cout << "File: " << outputFilename << " cannot be open!" << endl;
    return 1;
  }

  TTree *noiseTree = new TTree("noise", "noise mean and rms");
  noiseTree->Branch("Run", &fRun, "Run/I");
  noiseTree->Branch("Subrun", &fSubrun, "Subrun/I");
  noiseTree->Branch("Event", &fEvent, "Event/I");
  noiseTree->Branch("EventSeconds", fEventSeconds, "EventSeconds/I");
  noiseTree->Branch("EventNanoSeconds", fEventNanoSeconds, "EventNanoSeconds/I");
  noiseTree->Branch("MeanView0", &fMeanView0, "MeanView0/D");
  noiseTree->Branch("RMSView0", &fRMSView0, "RMSView0/D");
  noiseTree->Branch("MeanView1", &fMeanView1, "MeanView1/D");
  noiseTree->Branch("RMSView1", &fRMSView1, "RMSView1/D");
  noiseTree->Branch("ChannelsCoverage", &fChannelsCoverage, "ChannelsCoverage/I" );

  TH1D *hChMean = new TH1D("hChMean", ";Channel number; Mean (adc counts)", Ch_0+Ch_1, 0, Ch_0+Ch_1);
  TH1D *hChRMS = new TH1D("hChRMS", ";Channel number; Rms (adc counts)", Ch_0+Ch_1, 0, Ch_0+Ch_1);

  //here I initialize the parser object ========================================

  //subrun looper ==============================================================

  map<int, vector<double> > view2adc;
  map<int, int > channleMultiplicity;

  //for( size_t subrun = 0; subrun < nSubruns; subrun++ ){

    size_t subrun = nSubruns;
    string commonPath = "/eos/experiment/wa105/offline/LArSoft/Data/";
    string rawFile =  commonPath+"/Raw/ROOT/"+to_string(runNum)+"/"+to_string(runNum)+"-"+to_string(subrun)+"-RawWaveform-Parser.root" ;
    string recoFile =  commonPath+"Reco/2018_June_24/ROOT/recofast/"+to_string(runNum)+"/"+to_string(runNum)+"-"+to_string(subrun)+"-RecoFast-Parser.root";

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
      vector<Hit> recoHits;
      recoParser->getRecoHitsEvent( recoHits, evt );

      //raw objects
      vector<Channel> rawChannels;
      rawParser->getChannelsEvent( rawChannels, evt );

      //loop over channels
      for(auto rawChannel : rawChannels){

        fRun = rawChannel.run.getRunNumber();
        fSubrun = rawChannel.subRun;
        fEvent = rawChannel.event;
        fEventSeconds = rawChannel.timeSeconds;
        fEventNanoSeconds = rawChannel.timeNanoSeconds;

        int ch = rawChannel.channel;

        //find if the channel has a reco hit
        vector<Hit>:: iterator pos = find_if(recoHits.begin(), recoHits.end(),
        [ch]( Hit & recoHit ) -> bool { return  recoHit.channel == ch; });

        if (  pos == recoHits.end()  ){

          if ( !rawChannel.isDead() && channleMultiplicity[ ch ] == 0){

            double mean, rms;
            GetMeanAndRMS( rawChannel.signal, mean, rms);

            hChMean->Fill( ch, mean );
            hChRMS->Fill( ch, rms );

            //add the adc vector to the vector calculating the rms
            for ( auto adc : rawChannel.signal ){
              if ( adc < 150 ){
                view2adc[ rawChannel.view ].push_back( adc );
              }
            }

            //add the channel multiplicity if the condition is statisified
            channleMultiplicity[ ch ] = 1;

          } //end if bad and dead channels
        }//end if criteria
      } //end ch

    }//end event loop
  //}//end subrun loops

  //Get mean and rms =========================================================

  //calculate mean and rms for each view
  GetMeanAndRMS( view2adc[0], fMeanView0, fRMSView0 );
  GetMeanAndRMS( view2adc[1], fMeanView1, fRMSView1 );

  //get channel coverage ( number of unique channels without hits )
  for( auto it : channleMultiplicity ){
    if ( it.second > 0 )
      fChannelsCoverage++;
  }

  //fill the tree
  noiseTree->Fill();

  //fill and close file
  ofile->cd();

  ofile->Write();
  hChMean->Write();
  hChRMS->Write();

  ofile->Close();

  cout << "All done! " << endl;

  return 0;
}//end main
