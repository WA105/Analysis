////////////////////////////////////////////////////////////////////////////////
//
// 311 common analysis function and variables
//
////////////////////////////////////////////////////////////////////////////////
#ifndef  __311LIB_H
#define __311LIB_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile2D.h>
#include <TF1.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>

using namespace std;

//Geometry**********************************************************************

const int NUM_OF_VIEWS = 2;
const int NUM_OF_LEMS = 12;
int tpc_boundaries[6] = {-50, 50, -50, 50, 0, 300}; //minx,maxx,miny,maxy,minz,max
double Ch_0 = 320;
double Ch_1 = 960;
double pitch = 0.3;
int tdc = 1667;
double ADC2CHARGE = 45.31875; //ADC*ticks (from qScan)

//Data products ****************************************************************

class free_hit{
  public:
    int run;
    int subrun;
    int event;
    int view;
    int track_id;
    int channel;
    double adc_sum;
    double adc_integral;
    double peak_amp;
    double peak_time;
    double chi2;
};

class hit{
  public:
    int run;
    int subrun;
    int event;
    int view;
    int track_id;
    double adc_sum;
    double adc_integral;
    double peak_amp;
    double peak_time;
    double chi2;
    double sp_x;
    double sp_y;
    double sp_z;
    double dq;
    double dqdx;
    int lem;
};

class track{
  public:
    int run;
    int subrun;
    int event;
    int id;
    double start_x;
    double start_y;
    double start_z;
    double end_x;
    double end_y;
    double end_z;
    double length;
    double phi;
    double theta;
    int nhits;
    vector<free_hit> free_hits_trk;
    vector<hit> hits_trk;
};

//Functions declaration*********************************************************

int find_lem(double y, double z);

bool isGood_lem( int lem );

bool isGood_lem(vector<int> lems, int lem);

double find_projection(hit h);

void read_tree(TChain *rTree, vector<track> & tracks);

void drays_mitigation(track & t);

void select_mip(vector<track> tracks, vector<track> & mips,
                         vector<int> vol_cut, int length_cut, double angle_cut);

void select_tracks(vector<track> tracks, vector<track> & mips,
                         vector<int> vol_cut, int length_cut, double angle_cut);

double get_theta(track t);

double get_phi(track t);

double truncated_mean(std::vector<double> &vec, double rmlow, double rmhigh );

void FWHM(double &st, double &end, TH1D *h);

double find_max_bin(TH1D *hist);

double langaufun(double *x, double *par);

TF1 *langaufit(TH1D *his, double *fitrange, double *startvalues,
 double *parlimitslo, double *parlimitshi, double *fitparams, double *fiterrors,
                                      double *ChiSqr, int *NDF, bool find_best);

#endif // __311LIB_H
