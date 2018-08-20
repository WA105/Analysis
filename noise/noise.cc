////////////////////////////////////////////////////////////////////////////////
// Routine to calculate the noise mean and rms for a given subrun
//
// Input is a recofull file with noise hits
// Output is a TTree with the average noise rms and mean for that event in ADCs
// Output are a TH1D with the rms and the mean as function of the channel number
// Usage: ./build/noise -i path/to/input.root -o path/to/output.root -h no -p no
// -i input
// -o output
// -h exclude ROI checking at hits ( "yes" or "no" )
// -p remove pedestal ( "yes" or "no" )
//
// mailto:andrea.scarpelli@cern.ch
//
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

void GetMeanAndRMS(vector<float> array, float & mean, float & rms){
  //use welford method..consisted with qScan

  mean = 0;
  rms  = 0;

  if( array.size() == 0 )
    return;

  float A = 0;
  float Q = 0;

  for(size_t i=0;i<array.size();i++){

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

bool hasChannel( vector<Hit> hits ){
  //check if the vector exists

  if( !hits.size() ){ return false; }
  else { return true; }

}

bool inROI( float t, vector<Hit> hits ){
  //check if every time deposit is in a ROI

  //if no hit, also no ROI
  if( hits.size() == 0 ){
    return false;
  }

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

string getFileNumber( string filename ){
  //assuming the filename encorded in a format /path/to/file/run-subrun-Parser.root

  //isolate filenumber from path
  string s = filename;
  string delimiter = "/";

  size_t pos = 0;
  string token;

  while ((pos = s.find(delimiter)) != string::npos) {
    token = s.substr(0, pos);
    s.erase(0, pos + delimiter.length());
  }

  pos = s.find("-");
  string number = s.substr(0, pos);

  return number;
}

////////////////////////////////////////////////////////////////////////////////
// Main macro

int main( int argc, char* argv[] ){

  //parse the input argument ===================================================
  string recoFile;
  string outputFilename;
  string filter; //incoherent or coherent
  string useHits;
  string removePedestal;

  for ( int i=1; i<argc; i=i+2 )
  {
   if ( string(argv[i]) == "-i" ) recoFile = argv[i+1];
   else if ( string(argv[i]) == "-o" ) outputFilename = argv[i+1];
   else if ( string(argv[i]) == "-h" ) useHits = argv[i+1];
   else if ( string(argv[i]) == "-p" ) removePedestal = argv[i+1];
   else {
     cout << "Invalid option!" << endl;
     return 1;
   }
  }

  //define and variables =======================================================

  int fRun;
  int fSubrun;
  int fEvent;
  int fEventSeconds;
  int fEventNanoSeconds;
  float fMeanView1, fMeanView0;
  float fRMSView0, fRMSView1;

  //string outputFilename = filterName+"NoiseRms.root";
  TFile *ofile = new TFile(outputFilename.c_str(), "RECREATE");
  if(!ofile->IsOpen())
  {
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

  map<int, vector<float> > view2mean;
  map<int, vector<float> > view2rms;

  string runNum = getFileNumber( recoFile );

  string name = "_"+ runNum;
  string title = "Run: " + runNum;

  TH1D *hChMean = new TH1D( ("hChMean"+name).c_str(), (""+title).c_str(),
                                                    Ch_0+Ch_1, -Ch_0, Ch_1 );
  TH1D *hChRMS = new TH1D( ("hChRMS"+name).c_str(), (""+title).c_str(),
                                                    Ch_0+Ch_1, -Ch_0, Ch_1 );

  //Prepare inputs =============================================================

  LArParser *recoParser = new LArParser();
  TTree *recoTree = getTTree( recoFile );

  //check if the tree has been correctly set ===================================

  if( !recoTree){
    cout << "Error! Trees doesn't exist! " << endl;
    return 1;
  }

//Event looper =================================================================

map<int, vector<float>> fCh2Mean;
map<int, vector<float>> fCh2RMS;

int mod = 6;
if( recoTree->GetEntries() < 300 )
{
  mod = recoTree->GetEntries()/50; //it will process event by event only max 100 events

  //protect the case when mod=0
  if(mod == 0) { mod = 1; }

}

for(int evt=0; evt<recoTree->GetEntries(); evt++)
{
    if( evt % mod !=0 ){ continue; } //process one event every 6 ( about 50 events per subrun w 335 events )
      cout << " Processing event " << evt << endl;

      //reco objects -----------------------------------------------------------
      bool hasHits;

      if(useHits == "yes")
        hasHits=true;
      else
        hasHits=false;

      vector<Hit> recoHits;
      map<int, vector<Hit> > ch2hits;

      if(hasHits)
        recoParser->getRecoHitsEvent( recoTree,  recoHits, evt );

      if(!recoHits.size())
      {
        cout << "WARNING: No recoHits object created. ROIs may be accounted in the noise " << endl;
        hasHits=false;
      }
      else
      {
        //make an hit channels map
        for(auto recoHit : recoHits){
          ch2hits[ recoHit.channel ].push_back( recoHit );
        }
      }

      //raw objects
      vector<Channel> rawChannels;
      recoParser->getRecoChannelsEvent( recoTree, rawChannels, evt );

      if(!rawChannels.size()){
        //a channel info must be found in every event
        cout << "ERROR: No rawChannels object created." << endl;
        return 1;
      }

      //loop over channels -----------------------------------------------------
      for(auto rawChannel : rawChannels){

      if(removePedestal=="yes")
        rawChannel.subtractPedestal(true);
      else
        rawChannel.subtractPedestal(false);

         //skip dead or bad channels
        if ( rawChannel.isDead() || rawChannel.isBad() ){  continue; }

        vector<float> adc = rawChannel.signal; //real ADCs
        vector<float> noiseAdc; // only noise ADCs after ROI excluded

        //loop over the timesignal, skip the first two digits
        for( int t=2; t<(int)adc.size(); t++  )
        {
          if(hasHits)
          {
            if (  !inROI( t, ch2hits[ rawChannel.channel ] ) ) //exclude ROIs
            {
              noiseAdc.push_back( adc.at(t) );
            }
          }
          else //if no hits are found, there are no ROIs
          {
            noiseAdc.push_back( adc.at(t) );
          }
        } //end adc loop

        //calculate Mean and RMS here --------------------------------------------

        float mean, rms;
        GetMeanAndRMS( noiseAdc, mean, rms );

        if(rms > 0){
          view2mean[ rawChannel.view ].push_back( mean );
          view2rms[ rawChannel.view ].push_back( rms );

          fCh2Mean[rawChannel.channel].push_back(mean);
          fCh2RMS[rawChannel.channel].push_back(rms);

        }

        map<int, vector<float>>::iterator mapIt;
        for( mapIt = fCh2Mean.begin(); mapIt != fCh2Mean.end(); ++mapIt ){

          float mean = getMeanVector( fCh2Mean[ mapIt->first ] );
          float rms = getMeanVector( fCh2RMS[ mapIt->first ] );

          int channel = ViewToDAQChan( mapIt->first );

          if ( channel < 320 ){ channel = -channel; }
          hChMean->SetBinContent( channel, mean );
          hChRMS->SetBinContent( channel, rms );

        }

        noiseAdc.clear();

      } //end ch loop

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

      //fill the tree
      noiseTree->Fill();

      view2mean[0].clear();
      view2rms[0].clear();
      view2mean[1].clear();
      view2rms[1].clear();

  }//end event loop

  //fill and close file ========================================================
  ofile->cd();

  noiseTree->Write();
  hChMean->Write();
  hChRMS->Write();
  ofile->Close();

  cout << "All done! " << endl;

  return 0;

}//end main
