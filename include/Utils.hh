////////////////////////////////////////////////////////////////////////////////
// Utilities functions and classes
//
//mailto:andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////

#ifndef  __UTILS_H
#define __UTILS_H

#include <glob.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h> //to use struct stat

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TComplex.h"
#include "TFFTRealComplex.h"
#include "TFFTComplexReal.h"
#include "TVirtualFFT.h"


using namespace std;

// Spare functions /////////////////////////////////////////////////////////////

inline bool ExistTest (const std::string& name);

TTree *getTTree( string filename );

string getFileNumber( string filename );

double getModule( double x, double y , double z );

// FFT Class ///////////////////////////////////////////////////////////////////
class fftUtils{
  public:
    fftUtils(int time_samples, float sampling_freq, bool roundup = false);
    ~fftUtils();

    //getters
    int   GetTimeSamples();
    float GetSamplingFreq();
    int   GetFrequencySize();

    void DoTimeArray(vector<float> &time_vector);
    void FFT(vector<float> input, vector<TComplex> & output);
    void PowerSpectrum(vector<TComplex> &fftdata, vector<float> &freq,
                                          vector<float> &spectrum, bool skipdc);
  private:
    int fDetectorTimeSize;
    int fSize;
    float fSampligFreq; //in MHz
    int fFreqSize;
    size_t nhalf;
    float df;

    TFFTRealComplex *fFFT;
};

#endif // __UTILS_H
