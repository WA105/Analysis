
////////////////////////////////////////////////////////////////////////////////
// Routine to calculate the fft of a given run
//
// Input is a recofull file
// Output is a TFile with the average fft map
// Usage: ./build/fft -i path/to/input.root -o path/to/output.root -m no -p no
// -i input
// -o output
// -m mirror waveform: if yes artificially extends to match the first power of two after 1667. ( "yes" or "no" )
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
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2.h"
#include "TComplex.h"

//projects include
#include "Run.hh"
#include "DataStructure.hh"
#include "Geometry.hh"
#include "Utils.hh"

using namespace std;

int main( int argc, char* argv[] ){

  //parse the input argument ===================================================

  string recoFile, outFile;
  string mirrorWaveform;
  string removePedestal;

  int pattern = 32; //channel periodicity pattern

  for ( int i=1; i<argc; i=i+2 )
  {
   if      ( string(argv[i]) == "-i" ) recoFile= argv[i+1];
   else if ( string(argv[i]) == "-o" ) outFile = argv[i+1];
   else if ( string(argv[i]) == "-m" ) mirrorWaveform = argv[i+1];
   else if ( string(argv[i]) == "-p" ) removePedestal = argv[i+1];
   else {
     //PrintUsage();
     return 1;
   }
  }

  //Define FFT object holders ==================================================

  //fft
  bool mirror=false;
  if (mirrorWaveform == "yes"){ mirror = true; }

  fftUtils *fftUtil = new fftUtils(tdc, sampling_freq, mirror);

  int fft_points = fftUtil->GetTimeSamples();
  int freq_size = fftUtil->GetFrequencySize();
  float df = sampling_freq/fft_points;

  cout << "Max num of channesl: " << maxNumChannels << endl;
  cout << "Frequency size: " << freq_size << endl;
  cout << "Sampling frequency (MHz): " << sampling_freq << endl;
  cout << "Number of time samples: " << tdc  << endl;
  cout << "Number of fft points: " << fft_points << endl;
  cout << "resolution [MHz]: " << df << endl;

  // variables =================================================================

  TProfile *hPowerSpectrum[NUM_OF_VIEWS];
  TProfile2D *hPatternProfile2D[NUM_OF_VIEWS];
  TProfile2D *hFrequencyProfile2D[NUM_OF_VIEWS];
  TProfile2D *hFrequencyChannelMap;
  TH2D *hFreqencyRMSChannelMap;

  hFrequencyProfile2D[0] = new TProfile2D( "hFrequencyProfile2D_0" , "View 0; channels;f [MHz]" , Ch_0, 0, Ch_0, freq_size, 0., freq_size*df );
  hFrequencyProfile2D[1] = new TProfile2D( "hFrequencyProfile2D_1" , "View 1; channels;f [MHz]" , Ch_1, 0, Ch_1, freq_size, 0., freq_size*df );

  hPowerSpectrum[0]  = new TProfile( "PowerSpectrum_0" , "View 0; f [MHz]", freq_size, 0, freq_size*df );
  hPowerSpectrum[1]  = new TProfile( "PowerSpectrum_1" , "View 1; f [MHz]", freq_size, 0, freq_size*df );

  hPatternProfile2D[0] = new TProfile2D( "hPatternProfile2D_0" , "View 0 (32 ch); channels;f [MHz]" , (int)Ch_0/pattern, 0, (int)Ch_0/pattern, freq_size, 0., freq_size*df );
  hPatternProfile2D[1] = new TProfile2D( "hPatternProfile2D_1" , "View 1 (32 ch); channels;f [MHz]" , (int)Ch_1/pattern, 0, (int)Ch_1/pattern, freq_size, 0., freq_size*df );

  hFrequencyChannelMap = new TProfile2D( "hFrequencyChannelMap" , "Frequency vs. channel map; channels;f [MHz]" , Ch_0+Ch_1, 0, Ch_0+Ch_1, freq_size, 0., freq_size*df );
  hFreqencyRMSChannelMap = new TProfile2D( "hFrequencyRMSChannelMap" , "Frequency rms vs. channel map; channels;f [MHz]" , Ch_0+Ch_1, 0, Ch_0+Ch_1, freq_size, 0., freq_size*df );
  //Prepare inputs =============================================================

  LArParser *recoParser = new LArParser();
  TTree *recoTree = getTTree( recoFile );

  //check if the tree has been correctly set ===================================

  if( !recoTree )
  {
    cout << "Error! Trees doesn't exist! " << endl;
    return 1;
  }

//Event looper =================================================================

vector<float> timeVector;
vector<TComplex> fftdata;
vector<float> freq;
vector<float> spectrum;

int mod = 6;
if( recoTree->GetEntries() < 300 )
{
  mod = recoTree->GetEntries()/50; //it will process event by event only max 100 events

  //protect the case when mod=0
  if(mod == 0) { mod = 1; }

}

cout << "Start event looper" << endl;

for(int evt=0; evt<recoTree->GetEntries(); evt++)
{

    if( evt % mod != 0 ){ continue; } //process one event every 6 ( about 50 events per subrun w 335 events )
      cout << " Processing event " << evt << endl;

      //raw objects
      vector<Channel> channels;
      recoParser->getRecoChannelsEvent( recoTree, channels, evt );

      if(!channels.size()){
        cout << "ERROR:No rawChannels object created. Break event loop" << endl;
        break;
      }

      //loop over channels -----------------------------------------------------

      for(auto channel : channels)
      {

        timeVector.clear();
        fftdata.clear();
        freq.clear();
        spectrum.clear();

        int ch = channel.channel;

        if(removePedestal=="yes")
          channel.subtractPedestal(true);
        else
          channel.subtractPedestal(false);

        if ( channel.isDead() || channel.isBad() ){ continue; }

        timeVector.resize(channel.signal.size()-2);

        for(size_t t=2; t<channel.signal.size(); t++)
          timeVector.at(t-2) = channel.signal.at(t);

        fftUtil->DoTimeArray(timeVector);

        fftUtil->FFT(timeVector, fftdata);

        fftUtil->PowerSpectrum( fftdata, freq, spectrum, true);

        if (channel.view == 1){ ch = ch-320; }

        //loop over frequecies and fill histograms
        for(int f=0; f<(int)freq.size(); f++){

          //skip here unwanted frequecies or spectrum values
          if( !isnormal(spectrum.at(f)) ){ continue; }

          hPatternProfile2D[ channel.view ]->Fill( (int)ch/pattern , freq.at(f), spectrum.at(f) );
          hFrequencyProfile2D[channel.view]->Fill( ch , freq.at(f), spectrum.at(f) );
          hPowerSpectrum[channel.view]->Fill( freq.at(f), spectrum.at(f) );

          if (channel.view == 1)
            hFrequencyChannelMap->Fill( ch+320 , freq.at(f), spectrum.at(f) );
          else
            hFrequencyChannelMap->Fill( ch , freq.at(f), spectrum.at(f) );

        }//end f

      } //end ch loop

 }//end event loop

 //get the errors of the TProfile2D
 for(int xx = 0; xx<hFrequencyChannelMap->GetNbinsX(); xx++)
 {
   for(int yy = 0; yy<hFrequencyChannelMap->GetNbinsY(); yy++)
   {
     double ee = hFrequencyChannelMap->GetBinError( xx, yy );
     hFreqencyRMSChannelMap->SetBinContent(xx, yy, ee);
   }
 }

 //fill and close file =========================================================

 TFile *ofile = new TFile(outFile.c_str(), "RECREATE");
 if(!ofile->IsOpen())
 {
    cout << "File: " << outFile << " cannot be open!" << endl;
    return 1;
  }

 ofile->cd();

 hFrequencyProfile2D[0]->Write();
 hFrequencyProfile2D[1]->Write();

 hPatternProfile2D[0]->Write();
 hPatternProfile2D[1]->Write();

 hPowerSpectrum[0]->Write();
 hPowerSpectrum[1]->Write();

 hFrequencyChannelMap->Write();
 hFreqencyRMSChannelMap->Write();

 ofile->Close();

 return 0;

}//end main
