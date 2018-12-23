////////////////////////////////////////////////////////////////////////////////
// Data structures methods implementation
//
//mailto:andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include "TVector3.h"

#include "Run.hh"
#include "DataStructure.hh"
#include "Geometry.hh"
#include "Utils.hh"
//#include "Cuts.hh"

#include "TH1D.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
Channel::Channel(){}
Channel::~Channel(){}

bool Channel::isDead(){
  //tag for dead channels (channel containing only 0)

  vector<float>:: iterator pos = find_if( signal.begin(), signal.end(),
  []( float & adc ) -> bool { return  adc != 0; });

  return ( pos == signal.end() );

}

bool Channel::isBad(){

  //NB: ok for now...but this value shoudl depend on the pedestal

  //tag for bad channels: channels with weird adc values (+/- 100 adc everywhere)
  bool forAll = all_of( signal.begin(), signal.end(),
  []( float & adc ) -> bool { return  abs(adc) > 100; });

  return forAll;
}

void Channel::subtractPedestal( bool subtractPedestal ){
  //subtract pedestal

  int ipedfinal = 0;  // an integer for raw digits
	if (subtractPedestal){

		   const int PedestalIterations = 3;
		   double ped = 0.;
		   double pedthreshold[PedestalIterations] = {10};

		   for(int i=0; i < PedestalIterations; i++){
		      double csum = 0.;
		      double tickcounter = 0;
		      for (size_t itick=0;itick<tdc;++itick){

		          if( i == 0 || this->signal.at(itick) < ped + pedthreshold[i]){
			         csum += this->signal.at(itick);
			         tickcounter+=1;
			        }
		       }
		       if(tickcounter < 100) break;
		       csum /= tickcounter;
		       ped = csum;
		   }
		   ipedfinal = ped;

       for (size_t itick=0;itick<tdc;++itick){
         this->signal.at(itick) -= ipedfinal;
       }

	  }
    return;
  }

double Channel::sumAdcInROI( int startTime, int endTime ){
  //return the summed adc of the signal between startTime and endTime

  double sum=0;
  for( int t=startTime; t< endTime; t++ ){
    sum+= signal.at(t);
  }

  return sum;

}

////////////////////////////////////////////////////////////////////////////////
Hit::Hit(){}
Hit::~Hit(){}

void Hit::calibrateCharge( string type){
  //correct the charge with the appropriate calibration constants from ADC to fC

  //apply the nominal calibration if calo* is 0
  double calo0, calo1;
  if( calo0 == 0 ){ calo0 = getCalAmpConstant(0, type); }
  if( calo1 == 0 ){ calo1 = getCalAmpConstant(1, type); }

  if( view == 0 ){ chargeIntegral *= calo0; chargeSummedADC *= calo0; }
  else if( view == 1 ){ chargeIntegral *= calo1; chargeSummedADC *= calo1; }
  else{
    //don't correct and prompt an error message
    cout << "Hit::calibrateCharge ERROR! Invalid view " << endl;
  }

  return;
}

double Hit::getdQdx(){
  //get the dqdx

  double dqdx=0;
  if( dxLocalTrackDirection ==0 ){ return dqdx; }
  else{
    dqdx = chargeIntegral/dxLocalTrackDirection;
  }

  return dqdx;
}

void Hit::findLem(){
  //for the 311 geometry retun the lem nubmer associated to a certain geometry
  lem =-99999;


  //TODO Hardcoded z coodinates
 //fist of all check if the set of coodinate is valid (tolerance on boundaries)
 if( fabs(Z) > 300 || fabs(Y) > 50){
   cout << "Unknown coordinates: LEM not found!" << endl;
   return;
 }

  int module = (int)Z/100;
  int pos = (int)Z % 100; //I base the classification on the first module and then rescale

  if(pos < lem_size && Y >= 0 )       { lem = (module*lem_in_module) + 1; }
  else if( pos < lem_size && Y < 0 )  { lem = (module*lem_in_module) + 3; }
  else if( pos >= lem_size && Y >= 0 ){ lem = (module*lem_in_module) + 2; }
  else if( pos >= lem_size && Y < 0 ) { lem = (module*lem_in_module) + 4; }
  else {
    //cout << "Unknown coordinates: LEM not found!" << endl;
    return;
  }

  return;
}

bool Hit::isGoodLem( vector<int> lems ){
  //select only good lems among a list
  if (find(lems.begin(),lems.end(), lem) != lems.end())
     return true;
  else
      return false;
}

double Hit::electronsFromSumADC(string type)
{
  //covert ADC from sumADC into electrons
  double adc2fc = getCalAmpConstant(view, type);
  double charge = 0;
  charge = chargeSummedADC*adc2fc;

  return charge/electronCharge; //return the number of electrons
}

double Hit::electronsFromIntegral(string type)
{
  //covert ADC from Integral into electrons
  double adc2fc = getCalAmpConstant(view, type);
  double charge = 0;
  charge = chargeIntegral*adc2fc;

  return charge/electronCharge; //return the number of electrons
}

////////////////////////////////////////////////////////////////////////////////
Track::Track(){}
Track::~Track(){}

////////////////////////////////////////////////////////////////////////////////
MCTrack::MCTrack(){}
MCTrack::~MCTrack(){}

////////////////////////////////////////////////////////////////////////////////
LArParser::LArParser(){}

LArParser::~LArParser(){}

void LArParser::setChannelNoise( string model )
{
  fUseChannelNoise=true;

  //open noise model tfile
  TFile *file = new TFile(model.c_str(), "READ");

  if( !file->IsOpen() )
  {
    cout << "setChannelNoise(): Cannot open noise model" << endl;
    fUseChannelNoise = false;
    file->Close();
  }
  else
  {
    TH1D *noiseHist = (TH1D*)file->Get("hChRMS_729");
    if( !noiseHist )
    {
      cout << "setChannelNoise(): Cannot open noise model" << endl;
      fUseChannelNoise = false;
      file->Close();
    }
    else
    {
      for( int bin = 0; bin<noiseHist->GetNbinsX(); bin++ )
      {
        fChannelNoise[ bin ] = noiseHist->GetBinContent(bin);
      }
      file->Close();
    }
  }

  file->Close();
}

void LArParser::setRawBranches(TTree *tree){

  //Link brances in root file containg raw data

  tree->SetBranchAddress("Run",&tRun);
  tree->SetBranchAddress("Subrun",&tSubrun);
  tree->SetBranchAddress("EventNumberInRun",&tEventNumberInRun);
  tree->SetBranchAddress("EventTimeSeconds",&tEventTimeSeconds);
  tree->SetBranchAddress("EventTimeNanoseconds", &tEventTimeNanoseconds);
  tree->SetBranchAddress("RawWaveform_NumberOfChannels",&tRawWaveform_NumberOfChannels);
  tree->SetBranchAddress("RawWaveform_NumberOfTicks",&tRawWaveform_NumberOfTicks);
  tree->SetBranchAddress("RawWaveform_NumberOfTicksInAllChannels",&tRawWaveform_NumberOfTicksInAllChannels);
  tree->SetBranchAddress("RawWaveform_Channel",&tRawWaveform_Channel);
  tree->SetBranchAddress("RawWaveform_ADC",&tRawWaveform_ADC);

}

void LArParser::setMCBranches(TTree *tree){

  //Link branches in the ROOT file to variables
  //Metadata

  //g4 particles
  tree->SetBranchAddress("Run",&tRun);
  tree->SetBranchAddress("Subrun",&tSubrun);
  tree->SetBranchAddress("EventNumberInRun",&tEventNumberInRun);
  tree->SetBranchAddress("EventTimeSeconds",&tEventTimeSeconds);
  tree->SetBranchAddress("EventTimeNanoseconds", &tEventTimeNanoseconds);
  tree->SetBranchAddress("MCTruth_GEANT4_NumberOfParticles", &tNGeantTrackPerEvent);
  tree->SetBranchAddress("MCTruth_GEANT4_PDGCode",&tPdg);
  tree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_ParticleID",&tParticleId);
  tree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_StartEnergy",&tStartE);
  tree->SetBranchAddress("MCTruth_GEANT4_IsInTPCAV",&tIsInTPCAV);
  tree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_StartTime",&tStartTime);
  tree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_StartMomentum",&tMom);
  tree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_StartMomentum_X",&tMomX);
  tree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_StartMomentum_Y",&tMomY);
  tree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_StartMomentum_Z",&tMomZ);
  tree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_StartPoint_X",&tStartX);
  tree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_StartPoint_Y",&tStartY);
  tree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_StartPoint_Z",&tStartZ);
  tree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_EndPoint_X",&tEndX);
  tree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_EndPoint_Y",&tEndY);
  tree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_EndPoint_Z",&tEndZ);
  tree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_StartDirection_Phi",&tStartPhi);
  tree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_StartDirection_Theta",&tStartTheta);
  tree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_Pathlength",&tLengthAV);
  tree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_EndEnergy",&tEndE);

}

void LArParser::setRecoChannelBranches(TTree *tree){

  //event metadata
  tree->SetBranchAddress("Run",&tRun);
  tree->SetBranchAddress("Subrun",&tSubrun);
  tree->SetBranchAddress("EventNumberInRun",&tEventNumberInRun);
  tree->SetBranchAddress("EventTimeSeconds",&tEventTimeSeconds);
  tree->SetBranchAddress("EventTimeNanoseconds", &tEventTimeNanoseconds);

  tree->SetBranchAddress("RecoWaveforms_NumberOfChannels",&tRecoWaveform_NumberOfChannels);
  tree->SetBranchAddress("RecoWaveform_NTicks",&tRecoWaveform_NumberOfTicks);
  tree->SetBranchAddress("RecoWaveform_NumberOfTicksInAllChannels",&tRecoWaveform_NumberOfTicksInAllChannels);
  tree->SetBranchAddress("RecoWaveform_Channel",&tRecoWaveform_Channel);
  tree->SetBranchAddress("RecoWaveform_ADC",&tRecoWaveform_ADC);

}

void LArParser::setRecoHitBranches(TTree *tree){

    //event metadata
    tree->SetBranchAddress("Run",&tRun);
    tree->SetBranchAddress("Subrun",&tSubrun);
    tree->SetBranchAddress("EventNumberInRun",&tEventNumberInRun);
    tree->SetBranchAddress("EventTimeSeconds",&tEventTimeSeconds);
    tree->SetBranchAddress("EventTimeNanoseconds", &tEventTimeNanoseconds);

    tree->SetBranchAddress("NumberOfHits",&tNumberOfHits);
    tree->SetBranchAddress("Hit_TPC",&tHit_TPC);
    tree->SetBranchAddress("Hit_View",&tHit_View);
    tree->SetBranchAddress("Hit_Channel",&tHit_Channel);
    tree->SetBranchAddress("Hit_Amplitude",&tHit_Amplitude);
    tree->SetBranchAddress("Hit_ChargeSummedADC",&tHit_ChargeSummedADC);
    tree->SetBranchAddress("Hit_ChargeIntegral",&tHit_ChargeIntegral);
    tree->SetBranchAddress("Hit_PeakTime",&tHit_PeakTime);
    tree->SetBranchAddress("Hit_StartTime",&tHit_StartTime);
    tree->SetBranchAddress("Hit_EndTime",&tHit_EndTime);
    tree->SetBranchAddress("Hit_Width",&tHit_Width);
    tree->SetBranchAddress("Hit_GoodnessOfFit",&tHit_GoodnessOfFit);
    tree->SetBranchAddress("Hit_Multiplicity",&tHit_Multiplicity);
    tree->SetBranchAddress("Hit_TrackID",&tHit_TrackID);
    tree->SetBranchAddress("Hit_trueID",&tHit_particleID);
    tree->SetBranchAddress("Hit_trueEnergyMax",&tHit_TrueEnergy);
    tree->SetBranchAddress("Hit_trueEnergyFraction",&tHit_TrueEnergyFraction);
    tree->SetBranchAddress("Hit_ClusterID",&tHit_ClusterID);
    tree->SetBranchAddress("hit_chanElectrons", &tHit_chanElectrons);
    tree->SetBranchAddress("hit_chanMaxElectrons", &tHit_chanMaxElectrons);
    tree->SetBranchAddress("hit_bestChanIDE", &tHit_bestChanIDE);

    return;
}

void LArParser::setRecoClusterBranches( TTree *tree ){

  //event metadata
  tree->SetBranchAddress("Run",&tRun);
  tree->SetBranchAddress("Subrun",&tSubrun);
  tree->SetBranchAddress("EventNumberInRun",&tEventNumberInRun);
  tree->SetBranchAddress("EventTimeSeconds",&tEventTimeSeconds);
  tree->SetBranchAddress("EventTimeNanoseconds", &tEventTimeNanoseconds);

  /*
  treeTree->SetBranchAddress("NumberOfClusters",&tNumberOfClusters);
  treeTree->SetBranchAddress("ClusterID",&tClusterID);
  treeTree->SetBranchAddress("Cluster_View",&tCluster_View);
  treeTree->SetBranchAddress("Cluster_NumberOfHits",&tCluster_NumberOfHits);
  treeTree->SetBranchAddress("Cluster_ChargeIntegral",&tCluster_ChargeIntegral);
  treeTree->SetBranchAddress("Cluster_StartChannel",&tCluster_StartChannel);
  treeTree->SetBranchAddress("Cluster_StartTick",&tCluster_StartTick);
  treeTree->SetBranchAddress("Cluster_EndChannel",&tCluster_EndChannel);
  treeTree->SetBranchAddress("Cluster_EndTick",&tCluster_EndTick);
  treeTree->SetBranchAddress("Cluster_StartCharge",&tCluster_StartCharge);
  treeTree->SetBranchAddress("Cluster_StartAngle",&tCluster_StartAngle);
  treeTree->SetBranchAddress("Cluster_EndCharge",&tCluster_EndCharge);
  treeTree->SetBranchAddress("Cluster_EndAngle",&tCluster_EndAngle);
  */

  return;
}

void LArParser::setRecoTrackBranches(TTree *tree){

  //event metadata
  tree->SetBranchAddress("Run",&tRun);
  tree->SetBranchAddress("Subrun",&tSubrun);
  tree->SetBranchAddress("EventNumberInRun",&tEventNumberInRun);
  tree->SetBranchAddress("EventTimeSeconds",&tEventTimeSeconds);
  tree->SetBranchAddress("EventTimeNanoseconds", &tEventTimeNanoseconds);

  tree->SetBranchAddress("NumberOfTracks",&tNumberOfTracks);
  tree->SetBranchAddress("TrackID",&tTrackID);
  tree->SetBranchAddress("Track_NumberOfHits",&tTrack_NumberOfHits);
  tree->SetBranchAddress("Track_BestParticleID",&tTrack_BestParticleID);
  tree->SetBranchAddress("Track_Purity",&tTrack_Purity);
  tree->SetBranchAddress("Track_Completeness",&tTrack_Completeness);
  tree->SetBranchAddress("Track_Length_Trajectory",&tTrack_Length);
  tree->SetBranchAddress("Track_StartPoint_X", &tTrack_StartPoint_X);
  tree->SetBranchAddress("Track_StartPoint_Y", &tTrack_StartPoint_Y);
  tree->SetBranchAddress("Track_StartPoint_Z", &tTrack_StartPoint_Z);
  tree->SetBranchAddress("Track_StartPoint_DistanceToBoundary", &tTrack_StartPoint_DistanceToBoundary);
  tree->SetBranchAddress("Track_EndPoint_X", &tTrack_EndPoint_X);
  tree->SetBranchAddress("Track_EndPoint_Y", &tTrack_EndPoint_Y);
  tree->SetBranchAddress("Track_EndPoint_Z", &tTrack_EndPoint_Z);
  tree->SetBranchAddress("Track_EndPoint_DistanceToBoundary",&tTrack_EndPoint_DistanceToBoundary);
  tree->SetBranchAddress("Track_StartDirection_Theta",&tTrack_StartDirection_Theta);
  tree->SetBranchAddress("Track_StartDirection_Phi",&tTrack_StartDirection_Phi);
  tree->SetBranchAddress("Track_StartDirection_X", &tTrack_StartDirection_X);
  tree->SetBranchAddress("Track_StartDirection_Y", &tTrack_StartDirection_Y);
  tree->SetBranchAddress("Track_StartDirection_Z", &tTrack_StartDirection_Z);
  tree->SetBranchAddress("Track_EndDirection_Theta",&tTrack_EndDirection_Theta);
  tree->SetBranchAddress("Track_EndDirection_Phi",&tTrack_EndDirection_Phi);
  tree->SetBranchAddress("Track_EndDirection_X", &tTrack_EndDirection_X);
  tree->SetBranchAddress("Track_EndDirection_Y", &tTrack_EndDirection_Y);
  tree->SetBranchAddress("Track_EndDirection_Z", &tTrack_EndDirection_Z);
  tree->SetBranchAddress("Track_PitchInViews", &tTrack_PitchInViews);
  tree->SetBranchAddress("Track_NumberOfHitsPerView",&tTrack_NumberOfHitsPerView);
  tree->SetBranchAddress("Track_Hit_X", &tTrack_Hit_X);
  tree->SetBranchAddress("Track_Hit_Y", &tTrack_Hit_Y);
  tree->SetBranchAddress("Track_Hit_Z", &tTrack_Hit_Z);
  tree->SetBranchAddress("Track_Hit_ds_LocalTrackDirection", &tTrack_dx_LocalTrackDirection);
  tree->SetBranchAddress("Track_Hit_ds_3DPosition", &tTrack_dx_3DPosition);
  tree->SetBranchAddress("Track_Hit_TPC", &tTrack_Hit_TPC);
  tree->SetBranchAddress("Track_Hit_View", &tTrack_Hit_View);
  tree->SetBranchAddress("Track_Hit_Channel", &tTrack_Hit_Channel);
  tree->SetBranchAddress("Track_Hit_PeakTime", &tTrack_Hit_PeakTime);
  tree->SetBranchAddress("Track_Hit_Amplitude", &tTrack_Hit_Amplitude);
  tree->SetBranchAddress("Track_Hit_ChargeSummedADC", &tTrack_Hit_ChargeSummedADC);
  tree->SetBranchAddress("Track_Hit_ChargeIntegral", &tTrack_Hit_ChargeIntegral);
  tree->SetBranchAddress("Track_Hit_StartTime", &tTrack_Hit_StartTime);
  tree->SetBranchAddress("Track_Hit_EndTime", &tTrack_Hit_EndTime);
  tree->SetBranchAddress("Track_Hit_Width", &tTrack_Hit_Width);
  tree->SetBranchAddress("Track_Hit_GoodnessOfFit", &tTrack_Hit_GoodnessOfFit);
  tree->SetBranchAddress("Track_Hit_Multiplicity", &tTrack_Hit_Multiplicity);
  tree->SetBranchAddress("Track_Hit_trueID", &tTrack_Hit_particleID);
  tree->SetBranchAddress("Track_Hit_trueEnergyMax", &tTrack_Hit_TrueEnergy);
  tree->SetBranchAddress("Track_Hit_trueEnergyFraction", &tTrack_Hit_TrueEnergyFraction);
  tree->SetBranchAddress("Track_Hit_ChannelElectrons", &tTrack_Hit_ChannelElectrons);
  tree->SetBranchAddress("Track_Hit_MaxChannelElectrons", &tTrack_Hit_MaxChannelElectrons);
  tree->SetBranchAddress("Track_Hit_BestChanIDE", &tTrack_Hit_BestChanIDE);

  return;
}

bool LArParser::isTreeGood(TTree *tree){
  //check if the fTree object exists

  if(!tree)
    return false;
  else
    return true;

}

void LArParser::fillRawChannels( vector<Channel> & channels ){

  for(int j=0; j < tRawWaveform_NumberOfChannels; j++){

    Channel dummyChannel;

    for(int k=0; k < tRawWaveform_NumberOfTicks; k++){

      dummyChannel.run = tRun;
      dummyChannel.subRun = tSubrun;
      dummyChannel.event = tEventNumberInRun;
      dummyChannel.timeSeconds = tEventTimeSeconds;
      dummyChannel.timeNanoSeconds = tEventTimeNanoseconds;
      dummyChannel.channel = tRawWaveform_Channel[j];
      dummyChannel.nTicks = tRawWaveform_NumberOfTicks;
      dummyChannel.signal.push_back( tRawWaveform_ADC[j*tRawWaveform_NumberOfTicks+k] );

      if(dummyChannel.channel < Ch_0){
        dummyChannel.view = 0;
      }
      else{
        dummyChannel.view = 1;
      }
    } //end loop on tdc ticks

    channels.push_back( dummyChannel );
  }

  return;
}

void LArParser::fillRecoChannels( vector<Channel> & channels ){

  for(int j=0; j < tRecoWaveform_NumberOfChannels; j++){

    Channel dummyChannel;

    for(int k=0; k < tRecoWaveform_NumberOfTicks[j]; k++){

      dummyChannel.run = tRun;
      dummyChannel.subRun = tSubrun;
      dummyChannel.event = tEventNumberInRun;
      dummyChannel.timeSeconds = tEventTimeSeconds;
      dummyChannel.timeNanoSeconds = tEventTimeNanoseconds;
      dummyChannel.channel = tRecoWaveform_Channel[j];
      dummyChannel.nTicks = tRecoWaveform_NumberOfTicks[j];
      dummyChannel.signal.push_back( tRecoWaveform_ADC[j*tRecoWaveform_NumberOfTicks[j]+k] );

      if(dummyChannel.channel < Ch_0){
        dummyChannel.view = 0;
      }
      else{
        dummyChannel.view = 1;
      }
    } //end loop on tdc ticks

    channels.push_back( dummyChannel );
  }

  return;
}

void LArParser::fillMCTrack( vector<MCTrack> & tracks ){
  //fill the reco MCTrack data product

  for( int l=0; l<tNGeantTrackPerEvent; l++  ){

    MCTrack dummyTrack;

    dummyTrack.run = tRun;
    dummyTrack.subrun = tSubrun;
    dummyTrack.eventNumber = tEventNumberInRun;
    dummyTrack.pdgCode=tPdg[l];
    dummyTrack.particleID=tParticleId[l];
    dummyTrack.isInTPCAV=tIsInTPCAV[l];
    dummyTrack.startE=tStartE[l];
    dummyTrack.mom=tMom[l];
    dummyTrack.momX=tMomX[l];
    dummyTrack.momY=tMomY[l];
    dummyTrack.momZ=tMomZ[l];
    dummyTrack.startTime=tStartTime[l];


    if ( dummyTrack.isInTPCAV == 0 )
    {
      //those quantities are not inside the active volume initialize with -999
      dummyTrack.startX=-999;
      dummyTrack.startY=-999;
      dummyTrack.startZ=-999;
      dummyTrack.startPhi=-999;
      dummyTrack.startTheta=-999;
      dummyTrack.lengthAV=-999;
      dummyTrack.endE=-999;
      dummyTrack.endX=-999;
      dummyTrack.endY=-999;
      dummyTrack.endZ=-999;
      dummyTrack.length =-999;
      dummyTrack.startDirection.SetXYZ(-999,-999,-999);
    }
    else
    {
      dummyTrack.startX=tStartX[l];
      dummyTrack.startY=tStartY[l];
      dummyTrack.startZ=tStartZ[l];
      dummyTrack.startPhi=tStartPhi[l];
      dummyTrack.startTheta=tStartTheta[l];
      dummyTrack.lengthAV=tLengthAV[l];
      dummyTrack.endE=tEndE[l];
      dummyTrack.endX=tEndX[l];
      dummyTrack.endY=tEndY[l];
      dummyTrack.endZ=tEndZ[l];
      dummyTrack.length = getModule( abs( dummyTrack.endX - dummyTrack.startX ),
                                     abs( dummyTrack.endY - dummyTrack.startY ),
                                     abs( dummyTrack.endZ - dummyTrack.startZ )
                                   );
      dummyTrack.startDirection.SetXYZ( abs( dummyTrack.endX - dummyTrack.startX ) / dummyTrack.length,
                                        abs( dummyTrack.endY - dummyTrack.startY ) / dummyTrack.length,
                                        abs( dummyTrack.endZ - dummyTrack.startZ ) / dummyTrack.length
                                       );
    }


    tracks.push_back(dummyTrack);
  }

  return;
}

void LArParser::fillRecoHits( vector<Hit> & hits ){
  //fill the reco hits object data product

  for(int l=0; l<tNumberOfHits; l++){

    Hit dummyHits;

    dummyHits.run=tRun;
    dummyHits.subRun=tSubrun;
    dummyHits.event=tEventNumberInRun;
    dummyHits.TPC=tHit_TPC[l];
    dummyHits.view=tHit_View[l];
    dummyHits.channel=tHit_Channel[l];
    dummyHits.peakTime=tHit_PeakTime[l];
    dummyHits.amplitude=tHit_Amplitude[l];
    dummyHits.chargeSummedADC=tHit_ChargeSummedADC[l];
    dummyHits.chargeIntegral=tHit_ChargeIntegral[l];
    dummyHits.startTime=tHit_StartTime[l];
    dummyHits.endTime=tHit_EndTime[l];
    dummyHits.width=tHit_Width[l];
    dummyHits.goodnessOfFit=tHit_GoodnessOfFit[l];
    dummyHits.multiplicity=tHit_Multiplicity[l];
    dummyHits.trackID=tHit_TrackID[l];
    dummyHits.particleID= tHit_particleID[l];
    dummyHits.trueEnergy = tHit_TrueEnergy[l];
    dummyHits.trueEnergyFraction = tHit_TrueEnergyFraction[l];
    dummyHits.chanElectrons = tHit_chanElectrons[l];
    dummyHits.chanMaxElectrons = tHit_chanMaxElectrons[l];
    dummyHits.chanBestID = tHit_bestChanIDE[l];
    //dummyHit.lem=dummyHits.findLEM(dummyHits.Y, dummyHits.Z);

    if( fUseChannelNoise )
    {
      dummyHits.channelNoise = fChannelNoise[ dummyHits.channel ];
    }
    else{ dummyHits.channelNoise=0; }

    hits.push_back(dummyHits);
  }

  return;
}

void LArParser::fillRecoTrack( vector<Track> & tracks ){
  //Loop over track element and fill the data products

  int a=0; //Need this counter to remember where we are in this array.
  for(int j=0; j<tNumberOfTracks; j++) //Track loop
  {
    Track dummyTrack;

    dummyTrack.run = tRun;
    dummyTrack.subRun=tSubrun;
    dummyTrack.event=tEventNumberInRun;
    dummyTrack.trackID=tTrackID[j];
    dummyTrack.numberOfHits=tTrack_NumberOfHits[j];

    dummyTrack.bestParticleID = tTrack_BestParticleID[j];
    dummyTrack.completeness = tTrack_Completeness[j];
    dummyTrack.purity = tTrack_Purity[j];

    dummyTrack.length=tTrack_Length[j];
    dummyTrack.startPointX=tTrack_StartPoint_X[j];
    dummyTrack.startPointY=tTrack_StartPoint_Y[j];
    dummyTrack.startPointZ=tTrack_StartPoint_Z[j];
    dummyTrack.startPointDistanceToBoundary=tTrack_StartPoint_DistanceToBoundary[j];
    dummyTrack.endPointX=tTrack_EndPoint_X[j];
    dummyTrack.endPointY=tTrack_EndPoint_Y[j];
    dummyTrack.endPointZ=tTrack_EndPoint_Z[j];
    dummyTrack.endPointDistanceToBoundary=tTrack_EndPoint_DistanceToBoundary[j];
    dummyTrack.startDirectionTheta=tTrack_StartDirection_Theta[j];
    dummyTrack.startDirectionPhi=tTrack_StartDirection_Phi[j];
    dummyTrack.startDirectionX=tTrack_StartDirection_X[j];
    dummyTrack.startDirectionY=tTrack_StartDirection_Y[j];
    dummyTrack.startDirectionZ=tTrack_StartDirection_Z[j];
    dummyTrack.endDirectionTheta=tTrack_EndDirection_Theta[j];
    dummyTrack.endDirectionPhi=tTrack_EndDirection_Phi[j];
    dummyTrack.endDirectionX=tTrack_EndDirection_X[j];
    dummyTrack.endDirectionY=tTrack_EndDirection_Y[j];
    dummyTrack.endDirectionZ=tTrack_EndDirection_Z[j];
    dummyTrack.startDirection.SetXYZ( dummyTrack.startDirectionX,
                                      dummyTrack.startDirectionY,
                                      dummyTrack.startDirectionZ
                                    );
    dummyTrack.endDirection.SetXYZ( dummyTrack.endDirectionX,
                                      dummyTrack.endDirectionY,
                                      dummyTrack.endDirectionZ
                                    );

    for(int k=0; k<NUM_OF_VIEWS; k++) //View loop (for this track)
    {
      for(int l=a; l<a+tTrack_NumberOfHitsPerView[j][k]; l++) //Hit loop
      {
        Hit dummyHits;

        dummyHits.run = tRun;
        dummyHits.subRun=tSubrun;
        dummyHits.event=tEventNumberInRun;
        dummyHits.X=tTrack_Hit_X[l];
        dummyHits.Y=tTrack_Hit_Y[l];
        dummyHits.Z=tTrack_Hit_Z[l];
        dummyHits.dxLocalTrackDirection=tTrack_dx_LocalTrackDirection[l];
        dummyHits.dx3DPosition=tTrack_dx_3DPosition[l];
        dummyHits.TPC=tTrack_Hit_TPC[l];
        dummyHits.view=tTrack_Hit_View[l];
        dummyHits.channel=tTrack_Hit_Channel[l];
        dummyHits.peakTime=tTrack_Hit_PeakTime[l];
        dummyHits.amplitude=tTrack_Hit_Amplitude[l];
        dummyHits.chargeSummedADC=tTrack_Hit_ChargeSummedADC[l];
        dummyHits.chargeIntegral=tTrack_Hit_ChargeIntegral[l];
        dummyHits.startTime=tTrack_Hit_StartTime[l];
        dummyHits.endTime=tTrack_Hit_EndTime[l];
        dummyHits.width=tTrack_Hit_Width[l];
        dummyHits.goodnessOfFit=tTrack_Hit_GoodnessOfFit[l];
        dummyHits.multiplicity=tTrack_Hit_Multiplicity[l];
        dummyHits.particleID= tTrack_Hit_particleID[l];
        dummyHits.trueEnergy = tTrack_Hit_TrueEnergy[l];
        dummyHits.trueEnergyFraction = tTrack_Hit_TrueEnergyFraction[l];
        dummyHits.chanElectrons = tTrack_Hit_ChannelElectrons[l];
        dummyHits.chanMaxElectrons = tTrack_Hit_MaxChannelElectrons[l];
        dummyHits.chanBestID = tTrack_Hit_BestChanIDE[l];

        dummyHits.findLem(); //use 3D info to fetch

        if( fUseChannelNoise )
        {
          dummyHits.channelNoise = fChannelNoise[ dummyHits.channel ];
        }
        else{ dummyHits.channelNoise=0; }

        dummyTrack.hitsTrk.push_back(dummyHits);
      } //end hits
      a+=tTrack_NumberOfHitsPerView[j][k];;
    } //end view

  tracks.push_back(dummyTrack);
  }//end tracks

  return;
}

void LArParser::getRawChannelsEvent( TTree *tree,  vector<Channel> & channels, int event  ){
  //just fill the hit array for a specific event

  if ( this->isTreeGood(tree) ){

    this->setRawBranches(tree);

    tree->GetEntry(event);
    this->fillRawChannels( channels );

    //this->clean();

    return;

  }else{
    cout << "LArParser::getRawChannelsEvent ERROR:Tree doesn't exist!" << endl;
    return;
  }

}

void LArParser::getRecoChannelsEvent( TTree *tree,  vector<Channel> & channels, int event  ){
  //just fill the hit array for a specific event

  //fTree = tree;

  if ( this->isTreeGood( tree ) ){

      this->setRecoChannelBranches(tree);
      tree->GetEntry(event);
      this->fillRecoChannels( channels );
      //this->clean();

      return;

  }else{
    cout << "LArParser::getRecoChannelsEvent ERROR:Tree doesn't exist!" << endl;
    return;
  }

}

void LArParser::getMCTracksEvent(TTree *tree, vector<MCTrack> & tracks, int event  ){
  //just fill the hit array for a specific event

  if ( this->isTreeGood(tree) ){

    this->setMCBranches(tree);

    tree->GetEntry(event);
    this->fillMCTrack( tracks );

    return;

  }else{
    cout << "LArParser::getMCTracksEvent ERROR:Tree doesn't exist!" << endl;
    return;
  }

}

void LArParser::getRecoHitsEvent( TTree *tree,  vector<Hit> & hits, int event  ){
  //just fill the hit array for a specific event

  if ( this->isTreeGood(tree) ){

      this->setRecoHitBranches(tree);

      tree->GetEntry(event);
      this->fillRecoHits( hits );

      return;
  }else{
    cout << "LArParser::getRecoHitsEvent ERROR:Tree doesn't exist!" << endl;
    return;
  }

}

void LArParser::getRecoTracksEvent(TTree *tree, vector<Track> & tracks, int event  ){
  //just fill the hit array for a specific event

  this->setRecoTrackBranches(tree);

  tree->GetEntry(event);
  this->fillRecoTrack( tracks );

  return;

}

void LArParser::getMCTracks(TTree *tree, vector<MCTrack> & tracks ){
  //return an array holding the reconstruced track

  this->setMCBranches(tree);

  for(int i=0; i<tree->GetEntries(); i++) //Event loop
  {
    tree->GetEntry(i);
    this->fillMCTrack( tracks );
  }//Event loop

  return;
}

void LArParser::getRecoHits(TTree *tree, vector<Hit> & hits ){
  //return a vector holding the free reconstructed hits

  this->setRecoHitBranches(tree);

  for(int i=0; i<tree->GetEntries(); i++){ //Event loop

    tree->GetEntry(i);
    this->fillRecoHits( hits );
  }//Event loop

  return;
}

void LArParser::getRecoTracks(TTree *tree, vector<Track> & tracks ){
  //return an array holding the reconstruced track

  this->setRecoTrackBranches(tree);

  for(int i=0; i<tree->GetEntries(); i++) //Event loop
  {
    tree->GetEntry(i);
    this->fillRecoTrack( tracks );
  }//Event loop

  return;
}

//void LArParser::clean(){
  //fTree->Delete();
//}
