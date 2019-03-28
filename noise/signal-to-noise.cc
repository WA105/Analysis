////////////////////////////////////////////////////////////////////////////////
// Routine to calculate the noise mean and rms for a given subrun
//
// Input is a recofull file with noise hits
// Output is a TTree with the average noise rms and mean for that event in ADCs
// Output are a TH1D with the rms and the mean as function of the channel number
// Usage: ./build/noise -i path/to/input.root -o path/to/output.root -h no -p no
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
#include "TH2.h"

//projects include
#include "Run.hh"
#include "DataStructure.hh"
#include "Cuts.hh"
#include "Utils.hh"
#include "Geometry.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Main macro

int main( int argc, char* argv[] ){

  //parse the input argument ===================================================
  string recoFile, rawFile;
  string outputFilename;
  string removePedestal;

  for ( int i=1; i<argc; i=i+2 )
  {
   if ( string(argv[i]) == "-raw" ) rawFile = argv[i+1];
   else if ( string(argv[i]) == "-reco" ) recoFile = argv[i+1];
   else if ( string(argv[i]) == "-o" ) outputFilename = argv[i+1];
   else if ( string(argv[i]) == "-p" ) removePedestal = argv[i+1];
   else {
     cout << "Invalid option!" << endl;
     return 1;
   }
  }

  //define and variables =======================================================

  //string outputFilename = filterName+"NoiseRms.root";
  TFile *ofile = new TFile(outputFilename.c_str(), "RECREATE");
  if(!ofile->IsOpen())
  {
    cout << "File: " << outputFilename << " cannot be open!" << endl;
    return 1;
  }

  TH2D *hDirSummedADC = new TH2D("hDirSummedADC", "dir summed adc", 100, 0, 1, 100, 0, 1 );
  TH2D *hDirIntegral = new TH2D("hDirIntegral", " dir integral ", 100, 0, 1, 100, 0, 1 );

  TH2D* hChanSummedADC = new TH2D("hChanSummedADC", "chan summedAdc", Ch_0+Ch_1, 0, Ch_0+Ch_1, 100, 0, 1);
  TH2D* hChanIntegral = new TH2D("hChanIntegral", " chan integral", Ch_0+Ch_1, 0, Ch_0+Ch_1, 100, 0, 1 );

  //Prepare inputs =============================================================

  //raw
  LArParser *rawParser = new LArParser();
  TTree *rawTree = getTTree( rawFile );

  //reco
  LArParser *recoParser = new LArParser();
  TTree *recoTree = getTTree( recoFile );

  //check if the tree has been correctly set ===================================

  if( !recoTree && !rawTree){
    cout << "Error! Trees don't exist! " << endl;
    return 1;
  }

  if( recoTree->GetEntries() != rawTree->GetEntries() ){
    cout << "Error! Trees don't have the same number of entries! " << endl;
    return 1;
  }

//Event looper =================================================================

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

      //raw objects
      vector<Channel> rawChannels;
      rawParser->getRawChannelsEvent( rawTree, rawChannels, evt );

      if(!rawChannels.size()){
        //a channel info must be found in every event
        cout << "ERROR: No rawChannels object created." << endl;
        return 1;
      }

      //reco objects -----------------------------------------------------------
      vector<Hit> recoHits;
      vector<Track> recoTracks;
      map<int, vector<Hit> > ch2hits;

      recoParser->getRecoTracksEvent( recoTree,  recoTracks, evt );

      //loop over tracks -------------------------------------------------------
      for( auto recoTrack : recoTracks ){


        //dedicated track selection (truth matching and high purity-completeness)
        if( recoTrack.length < 100 ){ continue; }

        //fetch the hits in the track
        recoHits = recoTrack.hitsTrk;

        //get the track direction
        double dir = acos( abs(recoTrack.startPointX-recoTrack.endPointX)/recoTrack.length  );

        //loop over hits -------------------------------------------------------
        for( auto hit : recoHits ){

          Channel rawChannel;
          int ch = hit.channel;

          //find the signal array associated to the channel where the hit is found
          //now check if the time falls in the ROI identified by one hits
          vector<Channel>::iterator pos = find_if( rawChannels.begin(), rawChannels.end(),
          [ ch ]( Channel & rawChannel )->bool
            {
              return rawChannel.channel == ch;
            }
          );

          if( pos != rawChannels.end() ){
            rawChannel = (*pos);
          }else{
            cout << "ERROR: channel not found! check next hit" << endl;
            break;
          }

          //match hit with starting and ending time in ROI
          //(Bad, dead channels already selected in LArSoft sim)
          int startTime = hit.startTime;
          int endTime = hit.endTime;
          double integral = hit.chargeIntegral;
          double summedAdc = hit.chargeSummedADC;
          double rawSummedAdc = rawChannel.sumAdcInROI( startTime,  endTime );

          hDirSummedADC->Fill( dir, ( rawSummedAdc-summedAdc )/ rawSummedAdc );
          hDirIntegral->Fill( dir, ( rawSummedAdc-integral ) / rawSummedAdc );

          hChanSummedADC->Fill( ch, ( rawSummedAdc-summedAdc ) / rawSummedAdc );
          hChanIntegral->Fill( ch, ( rawSummedAdc-integral ) / rawSummedAdc );

          //fill plot with dqdx

        }//end loop hits
      }//end loop tracks

  }//end event loop

  //fill and close file ========================================================
  ofile->cd();

  hDirSummedADC->Write();
  hDirIntegral->Write();
  hChanSummedADC->Write();
  hChanIntegral->Write();

  ofile->Close();

  cout << "All done! " << endl;

  return 0;

}//end main
