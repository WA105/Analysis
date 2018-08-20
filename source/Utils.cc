////////////////////////////////////////////////////////////////////////////////
// Utilities functions and classes method implementation
//
//mailto:andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////

#include <glob.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h> //to use struct stat

#include "TFile.h"
#include "Utils.hh"
#include "TMath.h"
#include "TComplex.h"
#include "TFFTRealComplex.h"
#include "TFFTComplexReal.h"
#include "TVirtualFFT.h"

using namespace std;

// Spare functions /////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////////

TTree *getTTree( string filename ){
  //find the ttree associated to that specific filename: assuming fixed tree name

  TTree *ttree;

  if( ExistTest(filename) ){

      cout << "Processing file: " << filename << endl;
      TFile *file = new TFile( filename.c_str(), "READ" );

      if( file->IsOpen() )
        ttree = (TTree*)file->Get("analysistree/anatree");
      else
        cout << "getTTree::Error: Invalid TFile" << endl;

    }else{
      cout << "getTTree::Error: file doesn't exist " << filename << endl;
    }

    return ttree;
}

////////////////////////////////////////////////////////////////////////////////

fftUtils::fftUtils(int time_samples, float sampling_freq, bool roundup){

  fDetectorTimeSize = time_samples;
  if(roundup)
    fSize = (int)( pow(2, ceil(log(fDetectorTimeSize)/log(2))) );
  else
    fSize = (int)( pow(2, floor(log(fDetectorTimeSize)/log(2))) );
  fSampligFreq = sampling_freq;
  fFreqSize = fSize/2+1;
  nhalf = fSize/2;
  df = fSampligFreq/fSize;

  fFFT = new TFFTRealComplex(fSize, false);
  int dummy[1] = {0};

  fFFT->Init("  ", -1, dummy);
}

fftUtils::~fftUtils(){ delete fFFT; }

int fftUtils::GetTimeSamples(){
  return fSize;
}

float fftUtils::GetSamplingFreq(){
  return fSampligFreq;
}

int fftUtils::GetFrequencySize(){
  return fFreqSize;
}

void fftUtils::DoTimeArray(vector<float> &time_vector){
  /* Fill up the time array up to match the number of fft points (for 311:
  1024 if roundup is false, 2048 if roundup is true  */

  int time_points = (int) time_vector.size();

  if(time_points > fSize){ return; } //nothing else to do here.

  // x2 mirror waveform ////////////////////////////////////////////////////////
  //fetch the first n_window entry of the array and calculate mean in the window

  int n_window = 200; double sum=0;
  for(int i = time_points-n_window; i<time_points; i++){ sum+=time_vector[i]; }
  double shift = sum/n_window; //shift for the mean within the window
  shift *=2;

  time_vector.resize(fSize);

  //do x2 mirror (default)
  for(int t = time_points; t<2*time_points; t++)
    time_vector[t] = -time_vector[2*time_points - t - 1 ] + shift;

  return;
}

void fftUtils::FFT(vector<float> input, vector<TComplex> & output){

    //Does the fft for the input at channel ch
    double real      = 0.;
    double imaginary = 0.;

    for(size_t p = 0; p < (size_t)fSize ; ++p)
      {
        if( p < (size_t)fSize )
          fFFT->SetPoint(p, input[p]);
        else // zero pad
          fFFT->SetPoint(p, 0);
      }

    fFFT->Transform();

    if( (int)output.size() < fFreqSize) output.resize( fFreqSize );

    for(int i = 0; i < fFreqSize; ++i){
      fFFT->GetPointComplex(i, real, imaginary);
      output[i] = TComplex(real, imaginary);
    }

    return;
}

void fftUtils::PowerSpectrum(vector<TComplex> &fftdata, vector<float> &freq,
                                          vector<float> &spectrum, bool skipdc){

  spectrum.clear();
  freq.clear();

  double p0 =0;
  size_t istart = 0;
  if(skipdc) istart++; //skip the first tdc value (sometimes does crazy things)

  for( size_t i=istart;i<fftdata.size();i++)
    {
      double p = fftdata[i].Rho();
      double f = df*i;
    //  if( i!=0 && i!=nhalf ) p *= 2;
      if( i==0 ) p0 = p;

      freq.push_back(f);
      spectrum.push_back(p);
    }
    return;
}
