////////////////////////////////////////////////////////////////////////////////
//
// 311 common analysis function and variables
//
////////////////////////////////////////////////////////////////////////////////
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

//project libraries
#include <311Lib.h>

using namespace std;

//Utilities ********************************************************************

int find_lem(double y, double z){
  //for the 311 geometry retun the lem nubmer associated to a certain geometry
  int lem =-99999;
  int lem_in_module = 4; //num of minimal modules of 4 lems

 //fist of all check if the set of coodinate is valid (tolerance on boundaries)
 if( fabs(z) > 300 || fabs(y) > 50){
   //cout << "z " << z << " and y " << y <<" are unknown coordinates" << endl;
   return lem;
 }

  int module = (int)z/100;
  int pos = (int)z % 100; //I base the classification on the first module and then rescale

  if(pos < 50 && y >= 0 )       { lem = (module*lem_in_module) + 1; }
  else if( pos < 50 && y < 0 )  { lem = (module*lem_in_module) + 3; }
  else if( pos >= 50 && y >= 0 ){ lem = (module*lem_in_module) + 2; }
  else if( pos >= 50 && y < 0 ) { lem = (module*lem_in_module) + 4; }
  else {
    cout << "z " << pos << " y " << y <<" are unknown coordinate: LEM not found " << endl;
    return lem;
  }

  return lem;
}

bool isGood_lem( int lem ){
  //select only good lems among a list

  vector<int> lems = { 2, 4, 5, 6, 7, 8, 9, 11}; //active lems

   if (find(lems.begin(),lems.end(), lem) != lems.end())
     return true;
   else
      return false;
}

bool isGood_lem(vector<int> lems, int lem){
  //select only good lems among a list
  if (find(lems.begin(),lems.end(), lem) != lems.end())
     return true;
  else
      return false;
}

double find_projection(hit h){
  //return the correct variable projection on the view
  double coor = -99999;
  switch (h.view) {
    case 0:
      coor = h.sp_y;
      break;
    case 1:
      coor = h.sp_z;
      break;
    default:
      cout << " Unknown projection " << endl;
  }
  return coor;
}

double get_theta(track t){
  return acos( abs( t.end_x - t.start_x )/abs( t.length ) );
}

double get_phi(track t){
  if(abs(t.end_z - t.start_z)!=0)
    return atan( abs(t.end_y - t.start_y)/abs(t.end_z - t.start_z) );
  else
    return TMath::Pi()/2;
}

void read_tree(TChain *rTree, vector<track> & tracks){
  //read the tree and store all the tree information in the respective variables

  const int NMaxHitsPerEvent=100000;
  const int NMaxClustersPerEvent=10000;
  const int NMaxTracksPerEvent=1000;
  const int NMaxTracksPerEventTimesNViews=NMaxTracksPerEvent*NUM_OF_VIEWS;
  int NEventsPerRun=335;

  //Load ROOT file and tree
  int NEntries = (int)rTree->GetEntries();

  //Define variables to store the data of the ROOT file
  //Metadata
  int tRun;
  int tSubrun;
  int tEvent;
  int  tEventTime;
  char tIsData;

  //Hit variables
  int tNumberOfHits;
  short tHit_TrackID[NMaxHitsPerEvent];
  short tHit_View[NMaxHitsPerEvent];
  short tHit_Channel[NMaxHitsPerEvent];
  float tHit_ChargeIntegral[NMaxHitsPerEvent];
  float tHit_ChargeSummedADC[NMaxHitsPerEvent];
  float tHit_PeakHeight[NMaxHitsPerEvent];
  float tHit_PeakTime[NMaxHitsPerEvent];
  float tHit_GoodnessOfFit[NMaxHitsPerEvent];


  //Cluster variables
  short tNumberOfClusters;
  short tCluster_NumberOfHits[NMaxClustersPerEvent];

  //Track variables
  short tNumberOfTracks;
  float tTrack_Theta_pmtrack[NMaxTracksPerEvent];
  short tTrack_NumberOfHitsPerView_pmtrack[NMaxTracksPerEventTimesNViews];
  float tTrack_dQdx_pmtrack[NMaxHitsPerEvent];

  float tTrack_StartX_pmtrack[NMaxTracksPerEvent];
  float tTrack_StartY_pmtrack[NMaxTracksPerEvent];
  float tTrack_StartZ_pmtrack[NMaxTracksPerEvent];

  float tTrack_EndX_pmtrack[NMaxTracksPerEvent];
  float tTrack_EndY_pmtrack[NMaxTracksPerEvent];
  float tTrack_EndZ_pmtrack[NMaxTracksPerEvent];

  float tTrack_HitX_pmtrack[NMaxHitsPerEvent];
  float tTrack_HitY_pmtrack[NMaxHitsPerEvent];
  float tTrack_HitZ_pmtrack[NMaxHitsPerEvent];

  //Link branches in the ROOT file to variables
  //Metadata
  rTree->SetBranchAddress("Run",&tRun);
  rTree->SetBranchAddress("Subrun",&tSubrun);
  rTree->SetBranchAddress("Event",&tEvent);
  //rTree->SetBranchAddress("EventTimeSeconds",&tEventTime);
  rTree->SetBranchAddress("IsData",&tIsData);

  //Hit variables
  rTree->SetBranchAddress("NumberOfHits",&tNumberOfHits);
  rTree->SetBranchAddress("Hit_TrackID",&tHit_TrackID);
  rTree->SetBranchAddress("Hit_View",&tHit_View);
  rTree->SetBranchAddress("Hit_Channel",&tHit_Channel);
  rTree->SetBranchAddress("Hit_ChargeIntegral",&tHit_ChargeIntegral);
  rTree->SetBranchAddress("Hit_ChargeSummedADC",&tHit_ChargeSummedADC);
  rTree->SetBranchAddress("Hit_PeakHeight",&tHit_PeakHeight);
  rTree->SetBranchAddress("Hit_PeakTime",&tHit_PeakTime);
  rTree->SetBranchAddress("Hit_GoodnessOfFit",&tHit_GoodnessOfFit);

  //Cluster variables
  rTree->SetBranchAddress("NumberOfClusters",&tNumberOfClusters);
  rTree->SetBranchAddress("Cluster_NumberOfHits",&tCluster_NumberOfHits);

  //Track variables
  rTree->SetBranchAddress("NumberOfTracks_pmtrack",&tNumberOfTracks);
  rTree->SetBranchAddress("Track_Theta_pmtrack",&tTrack_Theta_pmtrack);
  rTree->SetBranchAddress("Track_NumberOfHitsPerView_pmtrack",&tTrack_NumberOfHitsPerView_pmtrack);
  rTree->SetBranchAddress("Track_dQdx_pmtrack",&tTrack_dQdx_pmtrack);
  rTree->SetBranchAddress("Track_StartX_pmtrack", &tTrack_StartX_pmtrack);
  rTree->SetBranchAddress("Track_StartY_pmtrack", &tTrack_StartY_pmtrack);
  rTree->SetBranchAddress("Track_StartZ_pmtrack", &tTrack_StartZ_pmtrack);
  rTree->SetBranchAddress("Track_EndX_pmtrack", &tTrack_EndX_pmtrack);
  rTree->SetBranchAddress("Track_EndY_pmtrack", &tTrack_EndY_pmtrack);
  rTree->SetBranchAddress("Track_EndZ_pmtrack", &tTrack_EndZ_pmtrack);
  rTree->SetBranchAddress("Track_HitX_pmtrack", &tTrack_HitX_pmtrack);
  rTree->SetBranchAddress("Track_HitY_pmtrack", &tTrack_HitY_pmtrack);
  rTree->SetBranchAddress("Track_HitZ_pmtrack", &tTrack_HitZ_pmtrack);

  for(int i=0; i<NEntries; i++) //Event loop
  {

    rTree->GetEntry(i);

    //initialize classes (not pointers for the moment)

    int a=0; //Need this counter to remember where we are in this array.
    for(int j=0; j<tNumberOfTracks; j++) //Track loop
    {
      track dummy_track;

      dummy_track.run = tRun;
      dummy_track.subrun = tSubrun;
      dummy_track.event = tEvent;
      dummy_track.id = j;
      dummy_track.start_x =tTrack_StartX_pmtrack[j];
      dummy_track.start_y =tTrack_StartY_pmtrack[j];
      dummy_track.start_z =tTrack_StartZ_pmtrack[j];

      dummy_track.end_x =tTrack_EndX_pmtrack[j];
      dummy_track.end_y =tTrack_EndY_pmtrack[j];
      dummy_track.end_z =tTrack_EndZ_pmtrack[j];

      double mag_x = pow( (dummy_track.end_x - dummy_track.start_x), 2);
      double mag_y = pow( (dummy_track.end_y - dummy_track.start_y) ,2);
      double mag_z = pow( (dummy_track.end_z - dummy_track.start_z) ,2);
      dummy_track.length = sqrt( mag_x + mag_y + mag_z );

      dummy_track.theta = get_theta(dummy_track);
      dummy_track.phi = get_phi(dummy_track);

      for(int k=0; k<tNumberOfHits; k++){

        if(dummy_track.id == tHit_TrackID[k] ){

          free_hit dummy_free_hit;

          dummy_free_hit.run = tRun;
          dummy_free_hit.subrun = tSubrun;
          dummy_free_hit.event = tEvent;
          dummy_free_hit.track_id = tHit_TrackID[k];
          dummy_free_hit.view = tHit_View[k];
          dummy_free_hit.channel = tHit_Channel[k];
          dummy_free_hit.adc_sum = tHit_ChargeSummedADC[k];
          dummy_free_hit.adc_integral = tHit_ChargeIntegral[k];
          dummy_free_hit.peak_amp = tHit_PeakHeight[k];
          dummy_free_hit.peak_time = tHit_PeakTime[k];
          dummy_free_hit.chi2 = tHit_GoodnessOfFit[k];

          dummy_track.free_hits_trk.push_back(dummy_free_hit);
        }
      }

      for(int k=0; k<NUM_OF_VIEWS; k++) //View loop (for this track)
      {
        for(int l=a; l<a+tTrack_NumberOfHitsPerView_pmtrack[2*j+k]; l++) //Hit loop
        {
          hit dummy_hits;

          dummy_hits.run = tRun;
          dummy_hits.subrun = tSubrun;
          dummy_hits.event = tEvent;
          dummy_hits.view = k;
          dummy_hits.track_id = dummy_track.id;
          dummy_hits.dqdx = tTrack_dQdx_pmtrack[l]/ADC2CHARGE;
          dummy_hits.sp_x = tTrack_HitX_pmtrack[l];
          dummy_hits.sp_y = tTrack_HitY_pmtrack[l];
          dummy_hits.sp_z = tTrack_HitZ_pmtrack[l];
          dummy_hits.lem  = find_lem(dummy_hits.sp_y, dummy_hits.sp_z);

          dummy_track.hits_trk.push_back(dummy_hits);
        } //end hits
        a+=tTrack_NumberOfHitsPerView_pmtrack[2*j+k];
      } //end view

    dummy_track.nhits = dummy_track.hits_trk.size();
    tracks.push_back(dummy_track);
    }//end tracks
  }//Event loop

}

//Cuts *************************************************************************

void drays_mitigation(track & t){

  //mitigate the effects of delta rays recontructed with the track removing
  //consecutve high dqds hits
  //take track as input, return smae track, but with delta rays hits padded to 0

  map<size_t, vector<double>> hits2view;

  unsigned int n_consecutive =15;
  double dqds_cut = 20;
  int initial_cut = 5;

  unsigned int c=0;
  double sum=0;
  int view = 0;
  int ii=0;

  vector<double> hit_list;

  for(unsigned int hh=0; hh<t.hits_trk.size(); hh++){

    auto h = t.hits_trk.at(hh);
    double dqds = h.dqdx;

    //assuming ordering first all hits view 0 then all hist view 1

    if(ii<initial_cut && h.view == view){ dqds=0; ++ii; }
    else if( ii==initial_cut ){ view++; ii=0; }

    hit_list.push_back( dqds );

  }

  for(unsigned int hh=0; hh<hit_list.size(); hh++){

    //0 pad consecutive hits above dqds_cut
    if(c<n_consecutive){
      sum+=hit_list.at(hh);
      ++c;
    }
    else if(c==n_consecutive){
      if(sum>dqds_cut*n_consecutive){
        fill(hit_list.begin()+hh-n_consecutive, hit_list.begin()+hh, 0);
      }
      c=0; sum=0;
    }
  }

  for( unsigned int hh=0; hh<t.hits_trk.size(); hh++ ){
    t.hits_trk.at(hh).dqdx = hit_list.at( hh );
  }

  return;
}//end function

void select_mip(vector<track> tracks, vector<track> & mips,
                vector<int> vol_cut = {2, 2, 2, 2, 2, 2}, int length_cut = 100,
                                                     double angle_cut = 3.15/2){
  //select particles crossing anode-cathode

  int count_mip=0;

  for(auto t : tracks){

    int minx = tpc_boundaries[0] + vol_cut[0];
    int maxx = tpc_boundaries[1] - vol_cut[1];
    int miny = tpc_boundaries[2] + vol_cut[2];
    int maxy = tpc_boundaries[3] - vol_cut[3];
    int minz = tpc_boundaries[4] + vol_cut[4];
    int maxz = tpc_boundaries[5] - vol_cut[5];

    double mag = sqrt( pow( (t.end_x - t.start_x) ,2)
            + pow( (t.end_y - t.start_y) ,2) + pow( (t.end_z - t.start_z) ,2) );

    if( mag > length_cut ) {
      //cut on the track angle
      if( t.theta > angle_cut ){ continue; }
      //check if the track is flipped on X
      if( max(t.end_x, t.start_x) > maxx && min(t.end_x, t.start_x) < minx ){

        //drays_mitigation(t);

        mips.push_back(t);
        count_mip++;
      } //end if x
      /*
      else if( max(t.end_y, t.start_y) > maxy && min(t.end_y, t.start_y) < miny ){
        mips.push_back(t);
        count_mip++;
      } //end if z
      else if( max(t.end_z, t.start_z) > maxz && min(t.end_z, t.start_z) < minz ){
        mips.push_back(t);
        count_mip++;
      } //end if y
      */
  } //end mag
 }//end for

  cout << " Selected " << count_mip << " mips over " << tracks.size()
                                                         << " tracks " <<  endl;
 return;

}

void select_tracks(vector<track> tracks, vector<track> & mips,
                vector<int> vol_cut = {2, 2, 2, 2, 2, 2}, int length_cut = 100,
                                                     double angle_cut = 3.15/2){
  //select particles crossing the detector in any direction

  int count_mip=0;

  for(auto t : tracks){

    int minx = tpc_boundaries[0] + vol_cut[0];
    int maxx = tpc_boundaries[1] - vol_cut[1];
    int miny = tpc_boundaries[2] + vol_cut[2];
    int maxy = tpc_boundaries[3] - vol_cut[3];
    int minz = tpc_boundaries[4] + vol_cut[4];
    int maxz = tpc_boundaries[5] - vol_cut[5];

    double mag = sqrt( pow( (t.end_x - t.start_x) ,2)
            + pow( (t.end_y - t.start_y) ,2) + pow( (t.end_z - t.start_z) ,2) );

    if( mag > length_cut ) {
      //cut on the track angle
      if( t.theta > angle_cut ){ continue; }
      //check if the track is flipped on X
      if( max(t.end_x, t.start_x) > maxx && min(t.end_x, t.start_x) < minx ){
        mips.push_back(t);
        count_mip++;
      } //end if x
      else if( max(t.end_y, t.start_y) > maxy && min(t.end_y, t.start_y) < miny ){
        mips.push_back(t);
        count_mip++;
      } //end if z
      else if( max(t.end_z, t.start_z) > maxz && min(t.end_z, t.start_z) < minz ){
        mips.push_back(t);
        count_mip++;
      } //end if y
  } //end mag
 }//end for

  cout << " Selected " << count_mip << " mips over " << tracks.size()
                                                         << " tracks " <<  endl;
 return;

}

//Analysis *********************************************************************

void FWHM(double &st, double &end, TH1D *h){
  //find fit paramter using FWHM. (same approach as uBoone)
  double x_max = find_max_bin( h );
  double bin_max = h->FindBin( x_max );
  double n_max = h->GetBinContent( bin_max );

  double th=0.5*n_max;

  for(int bin=0; bin<bin_max; bin++ ){
    if( th < h->GetBinContent(bin) ){
        st=h->GetBinCenter( bin-1 );
      break;
    }
  }

  for(int bin=bin_max; bin<h->GetNbinsX(); bin++ ){
    if( th > h->GetBinContent(bin) ){
      end=h->GetBinCenter( bin+1 );
      break;
    }
  }
  return;
}

double truncated_mean(std::vector<double> &vec, double rmlow, double rmhigh ){
  //**************************************************************************//
  // truncated mean
  //   rmlow  - fraction of events to remove at the low part of the distribution
  //   rmhigh - fraction of events to remove at the high part of the distribution
  //

 if(vec.empty()) return -999;

 std::vector<double> vtmp(vec);

 // sort it
 std::sort(vtmp.begin(), vtmp.end());

 size_t n = vtmp.size();
 size_t nlow, nhigh;

  //evaluate the number of entries to trim from the start and from the end
 // TODO: proceeding in this way, rounding the entry, we might be biased
 nlow  = (size_t)(rmlow * n);
 nhigh = n - (size_t)(rmhigh * n);

 // some basic checks
 if(rmlow <= 0 || rmlow >= 1) nlow = 0;
 if(rmhigh <= 0 || rmhigh >= 1) nhigh = n;

 if( (rmlow + rmhigh) >= 0.99)
   {
     cout<<"ERROR: truncated_mean() the fractions to remove are to high"<<endl;
     nlow  = 0;
     nhigh = n;
   }

 //cout<<“Subsample “<<nlow<<“, “<<nhigh<<” out of “<<n<<endl;

 // should not happen though
 if( nlow >= n ) nlow = 0;
 if( nhigh >= n ) nhigh = n;


 double trmean = 0.0;

 for(size_t i=nlow;i<nhigh;i++)
   trmean += vtmp[i];

 trmean /= ((float)nhigh - (float)nlow);
 return trmean;
}

double find_max_bin(TH1D *hist){
  //retunt the position of the bin with the largest number of entries
  int bin_max =0; int n_bin_max =0; double x_max =0.;

  for(int bin = 1; bin < (hist->GetSize()-2); bin++ ){
    if(hist->GetBinContent(bin) > n_bin_max){
      n_bin_max = hist->GetBinContent(bin);
      bin_max = bin;
    } //end if
  } //end for bin

  x_max = hist->GetXaxis()->GetBinCenter(bin_max);
  return x_max;
}

double langaufun(double *x, double *par){

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

TF1 *langaufit(TH1D *his, double *fitrange, double *startvalues,
 double *parlimitslo, double *parlimitshi, double *fitparams, double *fiterrors,
                              double *ChiSqr, int *NDF, bool find_best = false){
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
   //   If find_best == true, is returned the fit with the lower chisquare within multiple ranges.

   int i;
   char FunName[100];

   sprintf(FunName,"Fitfcn_%s",his->GetName());

   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
   if (ffitold) delete ffitold;

   TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
   ffit->SetParameters(startvalues);
   ffit->SetParNames("Width","MP","Area","GSigma");

   for (i=0; i<4; i++) {
      ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
   }

   //lower bound varies from 0.0*mean to 0.5*mean in septs of 0.05
   //upper bound varies from 0.9*mean to 3.1*mean in
   double min_chi =9999.; double min_j =0.; double min_i =0.;
   double x_max = find_max_bin(his);

   if(find_best){

     for(double i=0.1; i < fitrange[1]; i+=0.05){
       for(double j=fitrange[0]; j<3.0; j+=0.15){
         ffit->SetRange( i*his->GetMean(),j*his->GetMean() );
         his->Fit(FunName,"RBSQ");   // fit within specified range, use ParLimits, do not plot
         if( i*his->GetMean() < x_max && j*his->GetMean() > x_max ){
           if( ( ffit->GetChisquare()/ffit->GetNDF() ) < min_chi ){
             min_chi = ffit->GetChisquare()/ffit->GetNDF();
             min_j = j; min_i =i;
           }
         }
       }//end j
     }//end i
     ffit->SetRange( min_i*his->GetMean(),min_j*his->GetMean() );
     his->Fit(FunName,"LRBSQ");
   }//end find_best
   else{
     his->Fit(FunName,"LRBSQ");
   }

   ffit->GetParameters(fitparams);    // obtain fit parameters
   for (i=0; i<4; i++) {
      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
   }
   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
   NDF[0] = ffit->GetNDF();           // obtain ndf

   return (ffit);              // return fit function
}
