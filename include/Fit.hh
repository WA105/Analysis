////////////////////////////////////////////////////////////////////////////////
// Class handling the fit and the graph making for the field scans
//
// mailto: andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////

#ifndef  __FIT_H
#define __FIT_H

#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"

double langaufun(double *x, double *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

    // Numeric constants
    double invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
    double mpshift  = -0.22278298;       // Landau maximum location

    // Control constants
    double np = 100.0;      // number of convolution steps
    double sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

    // Variables
    double xx;
    double mpc;
    double fland;
    double sum = 0.0;
    double xlow,xupp;
    double step;
    double i;

    // MP shift correction
    mpc = par[1] - mpshift * par[0];

    // Range of convolution integral
    xlow = x[0] - sc * par[3];
    xupp = x[0] + sc * par[3];

    step = (xupp-xlow) / np;

    // Convolution integral of Landau and Gaussian by sum
    for(i=1.0; i<=np/2; i++) {
        xx = xlow + (i-.5) * step;
        fland = TMath::Landau(xx,mpc,par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0],xx,par[3]);

        xx = xupp - (i-.5) * step;
        fland = TMath::Landau(xx,mpc,par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0],xx,par[3]);
    }

    return (par[2] * step * sum * invsq2pi / par[3]);
}

TF1 *langaufit(TH1D *his, double *fitrange, double *startvalues, double *parlimitslo, double *parlimitshi, double *fitparams, double *fiterrors, double *ChiSqr, Int_t *NDF)
{
   // Once again, here are the Landau * Gaussian parameters:
   //   par[0]=Width (scale) parameter of Landau density
   //   par[1]=Most Probable (MP, location) parameter of Landau density
   //   par[2]=Total area (integral -inf to inf, normalization constant)
   //   par[3]=Width (sigma) of convoluted Gaussian function
   //
   // Variables for langaufit call:
   //   his             histogram to fit
   //   fitrange[2]     lo and hi boundaries of fit range
   //   startvalues[4]  reasonable start values for the fit
   //   parlimitslo[4]  lower parameter limits
   //   parlimitshi[4]  upper parameter limits
   //   fitparams[4]    returns the final fit parameters
   //   fiterrors[4]    returns the final fit errors
   //   ChiSqr          returns the chi square
   //   NDF             returns ndf

   Int_t i;
   Char_t FunName[100];

   sprintf(FunName,"Fitfcn_%s",his->GetName());

   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
   if (ffitold) delete ffitold;

   TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
   ffit->SetParameters(startvalues);
   ffit->SetParNames("Width","MP","Area","GSigma");

   for (i=0; i<4; i++) {
      ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
   }

   his->Fit(FunName,"RB0");   // fit within specified range, use ParLimits, do not plot

   ffit->GetParameters(fitparams);    // obtain fit parameters
   for (i=0; i<4; i++) {
      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
   }
   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
   NDF[0] = ffit->GetNDF();           // obtain ndf

   return (ffit);              // return fit function
}

////////////////////////////////////////////////////////////////////////////////

TF1* dofit( TH1D *hist, double lowlim, double uplim, bool norm)
{
  //Actual function implementing the fit described before

  Double_t fr[2];
  Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
  fr[0]=lowlim;
  fr[1]=uplim;

  cout << fr[0] << endl;
  cout << fr[1] << endl;

  pllo[0]=0.1*hist->GetStdDev();
  pllo[1]=0.1*hist->GetMean();
  if(norm){ pllo[2]=0.1; }
  else{ pllo[2]=  pllo[2]=0.1*hist->Integral(); }
  pllo[3]=0.1*hist->GetStdDev();

  plhi[0]=10*hist->GetStdDev();
  plhi[1]=10*hist->GetMean();
  if(norm){ plhi[2]=10; }
  else{ plhi[2]=10*hist->Integral(); }
  plhi[3]=10*hist->GetStdDev();

  sv[0]=hist->GetStdDev();
  sv[1]=hist->GetMean();
  if(norm){ sv[2]=1; }
  else{ sv[2]=hist->Integral(); }
  sv[3]=hist->GetStdDev();

  Double_t chisqr;
  Int_t    ndf;
  TF1 *fitsnr = langaufit(hist,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);

  return fitsnr;
}

//==============================================================================

class fitLandau{

  public:
    fitLandau();
    ~fitLandau();

    void setTitle(string title);
    void setName(string name);

    string getTitle();
    string getName();
    TH1D* getHist(int view);
    TF1* getFunction(int view);

    void initHist(string name, string title, int bin, double low, double high);
    void fillTH1(int view, double fillin);
    void doFit(int view, double start, double end, bool norm);
    double getFWHM(int view);

  protected:
    TH1D *dqds[2];
    TF1 *fitf[2];

    string fName;
    string fTitle;
};

//==============================================================================

class makeGraphs{

  public:
    makeGraphs(int n);
    ~makeGraphs();

    void setPointFit(int n, double binCenter, double binError, fitLandau* myHist);

    TGraphErrors *getGraphMPV(int view){ return gMPV[view]; };
    TGraphErrors *getGraphResolution(int view){ return gResolution[view]; };
    TGraphErrors *getGraphGain(){ return gGain; };
    TGraphErrors *getGraphAsymmetry(){ return gAsymmetry; };

    void fillGraphs(int n, double binCenter, double binError);

  protected:
    const double mipMPV=7.75;
    const double mipMEAN=10;

    TGraphErrors *gMPV[2];
    TGraphErrors *gResolution[2];
    TGraphErrors *gGain;
    TGraphErrors *gAsymmetry;

    TGraphErrors *gSlow;
    TGraphErrors *gFast;
    //TGraphErrors *gSum;

    TF1 *fSlow;
    TF1 *fFast;
    //TF1 *fSum;

    TGraphErrors *slowComponent;
    TGraphErrors *fastComponent;
    TGraphErrors *sumComponent;

    TF1 *fitf[2];

    double mpv[2];
    double sigma[2];
    double width[2];
    double gain;
    double asymmetry;
    double resolution[2];
    double empv[2];
    double esigma[2];
    double ewidth[2];
    double egain;
    double eresolution[2];
    double easymmetry;
};

//==============================================================================

#endif // __FIT_H
