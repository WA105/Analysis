////////////////////////////////////////////////////////////////////////////////
// The macro processes the input ROOT file specified in the variable filename
// and fits the dQ/ds distribution using the langauss function
//
// mailto:andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////

#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"

#include "TF1.h"
#include "Fit.hh"
#include "Geometry.hh"

TH1D *dqds[2];
TF1 *fit[2];

double phi[2] = {5,85}; //phi bin

//calorimeteric constants (uncomment the one you need)
double calo[2] = {59.859, 66.6802}; //ETHZ
//double calo[2] = {45.31875, 45.31875}; //IPNL
//double calo[2] = {42.2759, 59.3642}; //Monte Carlo

//Select here your filename
string filename = "/Users/scarpell/cernbox/signal_to_noise/840-dqds.root";

//==============================================================================
//Estetic functions

string round_and_convert( double num, int precision )
{

    string str = to_string(num);
    size_t pos = str.find(".");
    string str1 = str.substr( 0, pos+1 );
    string str2 = str.substr( pos+1, pos+precision );

    return str1+str2;
}

void makePretty(TCanvas *c, int canvasNum,  TH1D *hist, TF1 *ff, string view, Color_t color, double *trunc)
{
  c->cd(canvasNum);
  hist->SetStats(00000);
  hist->SetMarkerSize(0.6);
  hist->SetLineWidth(2.0);

  hist->GetXaxis()->SetTitle("dQ/ds [fC/cm]");
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->SetTitle("Hits");
  hist->GetYaxis()->CenterTitle();

  ff->SetLineWidth(3.0);
  ff->SetLineColor(color);

  hist->Draw("E");
  ff->Draw("SAME");

  double mpv = ff->GetParameter(1);
  double empv = ff->GetParError(1);
  double width = ff->GetParameter(0);
  double ewidth = ff->GetParError(0);
  double sigma = ff->GetParameter(3);
  double esigma = ff->GetParError(3);

  string smpv = "#scale[2.2]{MPV = "+round_and_convert(mpv, 1)+" #pm "+round_and_convert(empv, 1)+" fC/cm}";
  string swidth ="#scale[2.2]{Width = "+round_and_convert(width, 1)+" #pm "+round_and_convert(ewidth, 1)+" fC/cm}";
  string ssigma ="#scale[2.2]{Sigma = "+round_and_convert(sigma, 1)+" #pm "+round_and_convert(esigma, 1)+" fC/cm}";
  string strunc = "#scale[2.2]{Tr. Mean = "+round_and_convert(trunc[0], 1)+" #pm "+round_and_convert(trunc[1], 1)+" fC/cm}";

  TLegend *l1 = new TLegend(0.457206, 0.463265, 0.90602, 0.906803);
  l1->SetHeader(("#scale[3.0]{#bf{"+view+"}}").c_str(), "C");
  l1->AddEntry((TObject*)0,smpv.c_str(), "");
  l1->AddEntry((TObject*)0,swidth.c_str(), "");
  l1->AddEntry((TObject*)0,ssigma.c_str(), "");
  l1->AddEntry((TObject*)0,strunc.c_str(), "");
  l1->Draw("");
}

//==============================================================================
// Loop over the ttree entries

void loopTree(TTree *tree, fitLandau *hist, double philim[2])
{
  int view;
  double dqds;
  double integral;
  double ds;
  double phi;

  tree->SetBranchAddress("dQds", &dqds);
  tree->SetBranchAddress("Integral", &integral);
  tree->SetBranchAddress("fLocaldS", &ds);
  tree->SetBranchAddress("View", &view);
  tree->SetBranchAddress("StartPhi", &phi);
  for(int i=0; i<tree->GetEntries(); i++)
  {
    tree->GetEntry(i);
    if(rescalePhi(phi)>philim[0] && rescalePhi(phi)<philim[1])
    {
      hist->fillTH1(view, dqds/calo[view]);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
//Here's the actual macro

void dQds()
{
  bool norm=false;
  gStyle->SetPalette(kLightTemperature);

  TFile *_file0 = TFile::Open(filename.c_str());
  TTree *ttree = (TTree*) _file0->Get("hits");

  //Initialize the histogram using the fit class fitLandau in Fit.hh
  fitLandau *myHist = new fitLandau();
  string name = "dQds";
  string title = ";dQ/ds [fC/cm];Arb. Units";
  myHist->initHist(name.c_str(), title.c_str(), 50, 0, 30);

  //fill the histogram
  loopTree(ttree, myHist, phi);

  myHist->doFit(0, 0.5, 1.5, norm);
  myHist->doFit(1, 0.5, 1.5, norm);

  //Draw =======================================================================

  TCanvas *c = new TCanvas("c", "", 1000, 400);
  c->Divide(2,1);

  //not correct values
  double trunc0[2] = {492.2158241909493/calo[0], 0.844720554183213/calo[0]};
  double trunc1[2] = {483.2593816255622/calo[1], 0.48157343787501106/calo[1]};

  makePretty(c, 1, myHist->getHist(0), myHist->getFunction(0), "View 0", kRed+1, trunc0);
  makePretty(c, 2, myHist->getHist(1), myHist->getFunction(1), "View 1", kBlue+1, trunc1);
}
