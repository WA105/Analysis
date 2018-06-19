////////////////////////////////////////////////////////////////////////////////
// Data structures methods implementation
//
//mailto:andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////

//TODO add missing branches to MCTrack
//TODO add backtracker quantities to Track

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include "Run.hh"
#include "DataStructure.hh"
//#include "Cuts.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
Hit::Hit(){}
Hit::~Hit(){}

////////////////////////////////////////////////////////////////////////////////
Track::Track(){}
Track::~Track(){}

////////////////////////////////////////////////////////////////////////////////
MCTrack::MCTrack(){}
MCTrack::~MCTrack(){}

////////////////////////////////////////////////////////////////////////////////
LArParser::LArParser(){}

LArParser::LArParser( TTree *tree, Run *run ){
  //alternative constructor needs a tree and a run object
  fTree = tree;
  fRun = run;
}

LArParser::~LArParser(){}

void LArParser::setMCBranches(){

  //Link branches in the ROOT file to variables
  //Metadata
  fTree->SetBranchAddress("Subrun", &tSubrun);
  fTree->SetBranchAddress("EventNumberInRun", &tEventNumberInRun);

  //g4 particles
  fTree->SetBranchAddress("MCTruth_GEANT4_NumberOfParticles", &tNGeantTrackPerEvent);
  fTree->SetBranchAddress("MCTruth_GEANT4_PDGCode",&tPdg);
  fTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_ParticleID",&tParticleId);
  fTree->SetBranchAddress("MCTruth_GEANT4_StartEnergy",&tStartE);
  fTree->SetBranchAddress("MCTruth_GEANT4_IsInTPCAV",&tIsInTPCAV);
  fTree->SetBranchAddress("MCTruth_GEANT4_StartTime",&tStartTime);
  fTree->SetBranchAddress("MCTruth_GEANT4_StartMomentum",&tMom);
  fTree->SetBranchAddress("MCTruth_GEANT4_StartMomentum_X",&tMomX);
  fTree->SetBranchAddress("MCTruth_GEANT4_StartMomentum_Y",&tMomY);
  fTree->SetBranchAddress("MCTruth_GEANT4_StartMomentum_Z",&tMomZ);
  fTree->SetBranchAddress("MCTruth_GEANT4_StartPoint_X",&tStartX);
  fTree->SetBranchAddress("MCTruth_GEANT4_StartPoint_Y",&tStartY);
  fTree->SetBranchAddress("MCTruth_GEANT4_StartPoint_Z",&tStartZ);
  fTree->SetBranchAddress("MCTruth_GEANT4_InTPCAV_Pathlength",&tLengthAV);

}

void LArParser::setRecoBranches(){

    //Link branches in the ROOT file to variables
    //Metadata
    fTree->SetBranchAddress("Run",&tRun);
    fTree->SetBranchAddress("Subrun",&tSubrun);
    fTree->SetBranchAddress("EventNumberInRun",&tEventNumberInRun);
    //fTree->SetBranchAddress("EventTimeSeconds",&tEventTimeSeconds);
    //fTree->SetBranchAddress("EventTimeNanoseconds",&tEventTimeNanoseconds);
    //fTree->SetBranchAddress("IsData",&tIsData);

    //Hit variables
    fTree->SetBranchAddress("NumberOfHits",&tNumberOfHits);
    fTree->SetBranchAddress("Hit_TPC",&tHit_TPC);
    fTree->SetBranchAddress("Hit_View",&tHit_View);
    fTree->SetBranchAddress("Hit_Channel",&tHit_Channel);
    fTree->SetBranchAddress("Hit_ChargeSummedADC",&tHit_ChargeSummedADC);
    fTree->SetBranchAddress("Hit_ChargeIntegral",&tHit_ChargeIntegral);
    fTree->SetBranchAddress("Hit_PeakTime",&tHit_PeakTime);
    fTree->SetBranchAddress("Hit_StartTime",&tHit_StartTime);
    fTree->SetBranchAddress("Hit_EndTime",&tHit_EndTime);
    fTree->SetBranchAddress("Hit_Width",&tHit_Width);
    fTree->SetBranchAddress("Hit_GoodnessOfFit",&tHit_GoodnessOfFit);
    fTree->SetBranchAddress("Hit_Multiplicity",&tHit_Multiplicity);
    fTree->SetBranchAddress("Hit_TrackID",&tHit_TrackID);
    fTree->SetBranchAddress("Hit_trueID",&tHit_particleID);
    fTree->SetBranchAddress("hit_trueEnergyMax",&tHit_TrueEnergy);
    fTree->SetBranchAddress("hit_trueEnergyFraction",&tHit_TrueEnergy);
    //fTree->SetBranchAddress("Hit_ClusterID",&tHit_ClusterID);

    //  //Cluster variables
    //  fTree->SetBranchAddress("NumberOfClusters",&tNumberOfClusters);
    //  fTree->SetBranchAddress("ClusterID",&tClusterID);
    //  fTree->SetBranchAddress("Cluster_View",&tCluster_View);
    //  fTree->SetBranchAddress("Cluster_NumberOfHits",&tCluster_NumberOfHits);
    //  fTree->SetBranchAddress("Cluster_ChargeIntegral",&tCluster_ChargeIntegral);
    //  fTree->SetBranchAddress("Cluster_StartChannel",&tCluster_StartChannel);
    //  fTree->SetBranchAddress("Cluster_StartTick",&tCluster_StartTick);
    //  fTree->SetBranchAddress("Cluster_EndChannel",&tCluster_EndChannel);
    //  fTree->SetBranchAddress("Cluster_EndTick",&tCluster_EndTick);
    //  fTree->SetBranchAddress("Cluster_StartCharge",&tCluster_StartCharge);
    //  fTree->SetBranchAddress("Cluster_StartAngle",&tCluster_StartAngle);
    //  fTree->SetBranchAddress("Cluster_EndCharge",&tCluster_EndCharge);
    //  fTree->SetBranchAddress("Cluster_EndAngle",&tCluster_EndAngle);

    //  //Track variables
    fTree->SetBranchAddress("NumberOfTracks",&tNumberOfTracks);
    fTree->SetBranchAddress("TrackID",&tTrackID);
    fTree->SetBranchAddress("Track_NumberOfHits",&tTrack_NumberOfHits);
    fTree->SetBranchAddress("Track_Length_Trajectory",&tTrack_Length);

    fTree->SetBranchAddress("Track_StartPoint_X", &tTrack_StartPoint_X);
    fTree->SetBranchAddress("Track_StartPoint_Y", &tTrack_StartPoint_Y);
    fTree->SetBranchAddress("Track_StartPoint_Z", &tTrack_StartPoint_Z);
    fTree->SetBranchAddress("Track_StartPoint_DistanceToBoundary", &tTrack_StartPoint_DistanceToBoundary);
    fTree->SetBranchAddress("Track_EndPoint_X", &tTrack_EndPoint_X);
    fTree->SetBranchAddress("Track_EndPoint_Y", &tTrack_EndPoint_Y);
    fTree->SetBranchAddress("Track_EndPoint_Z", &tTrack_EndPoint_Z);
    fTree->SetBranchAddress("Track_EndPoint_DistanceToBoundary",&tTrack_EndPoint_DistanceToBoundary);

    fTree->SetBranchAddress("Track_StartDirection_Theta",&tTrack_StartDirection_Theta);
    fTree->SetBranchAddress("Track_StartDirection_Phi",&tTrack_StartDirection_Phi);
    fTree->SetBranchAddress("Track_StartDirection_X", &tTrack_StartDirection_X);
    fTree->SetBranchAddress("Track_StartDirection_Y", &tTrack_StartDirection_Y);
    fTree->SetBranchAddress("Track_StartDirection_Z", &tTrack_StartDirection_Z);

    fTree->SetBranchAddress("Track_EndDirection_Theta",&tTrack_EndDirection_Theta);
    fTree->SetBranchAddress("Track_EndDirection_Phi",&tTrack_EndDirection_Phi);
    fTree->SetBranchAddress("Track_EndDirection_X", &tTrack_EndDirection_X);
    fTree->SetBranchAddress("Track_EndDirection_Y", &tTrack_EndDirection_Y);
    fTree->SetBranchAddress("Track_EndDirection_Z", &tTrack_EndDirection_Z);

    fTree->SetBranchAddress("Track_PitchInViews", &tTrack_PitchInViews);
    fTree->SetBranchAddress("Track_NumberOfHitsPerView",&tTrack_NumberOfHitsPerView);

  //  //Track hit variables
    fTree->SetBranchAddress("Track_Hit_X", &tTrack_Hit_X);
    fTree->SetBranchAddress("Track_Hit_Y", &tTrack_Hit_Y);
    fTree->SetBranchAddress("Track_Hit_Z", &tTrack_Hit_Z);
    fTree->SetBranchAddress("Track_Hit_ds_LocalTrackDirection", &tTrack_dx_LocalTrackDirection);
    fTree->SetBranchAddress("Track_Hit_ds_3DPosition", &tTrack_dx_3DPosition);
    fTree->SetBranchAddress("Track_Hit_TPC", &tTrack_Hit_TPC);
    fTree->SetBranchAddress("Track_Hit_View", &tTrack_Hit_View);
    fTree->SetBranchAddress("Track_Hit_Channel", &tTrack_Hit_Channel);
    fTree->SetBranchAddress("Track_Hit_PeakTime", &tTrack_Hit_PeakTime);
    fTree->SetBranchAddress("Track_Hit_ChargeSummedADC", &tTrack_Hit_ChargeSummedADC);
    fTree->SetBranchAddress("Track_Hit_ChargeIntegral", &tTrack_Hit_ChargeIntegral);
    fTree->SetBranchAddress("Track_Hit_StartTime", &tTrack_Hit_StartTime);
    fTree->SetBranchAddress("Track_Hit_EndTime", &tTrack_Hit_EndTime);
    fTree->SetBranchAddress("Track_Hit_Width", &tTrack_Hit_Width);
    fTree->SetBranchAddress("Track_Hit_GoodnessOfFit", &tTrack_Hit_GoodnessOfFit);
    fTree->SetBranchAddress("Track_Hit_Multiplicity", &tTrack_Hit_Multiplicity);
    fTree->SetBranchAddress("Track_Hit_trueID", &tTrack_Hit_particleID);
    fTree->SetBranchAddress("Track_Hit_trueEnergyMax", &tTrack_Hit_TrueEnergy);
    fTree->SetBranchAddress("Track_Hit_trueEnergyFraction", &tTrack_Hit_TrueEnergyFraction);

    return;
}

bool LArParser::isTreeGood(){
  //check if the fTree object exists

  if(!fTree)
    return false;
  else
    return true;

}

void LArParser::fillMCTrack( vector<MCTrack> & tracks ){
  //fill the reco MCTrack data product

  for( int l=0; l<tNGeantTrackPerEvent; l++  ){

    MCTrack dummyTrack;

    dummyTrack.subrun = tSubrun;
    dummyTrack.event = tEventNumberInRun;

    dummyTrack.pdgCode=tPdg[l];
    dummyTrack.particleId=tParticleId[l];
    dummyTrack.isInTPCAV=tIsInTPCAV[l];
    dummyTrack.startE=tStartE[l];
    dummyTrack.mom=tMom[l];
    dummyTrack.momX=tMomX[l];
    dummyTrack.momY=tMomY[l];
    dummyTrack.momZ=tMomZ[l];
    dummyTrack.startTime=tStartTime[l];
    dummyTrack.startX=tStartX[l];
    dummyTrack.startY=tStartY[l];
    dummyTrack.startZ=tStartZ[l];
    dummyTrack.lengthAV=tLengthAV[l];
    dummyTrack.endX=tEndX[l];
    dummyTrack.endY=tEndY[l];
    dummyTrack.endZ=tEndZ[l];

    tracks.push_back(dummyTrack);
  }

  return;
}

void LArParser::fillRecoHits( vector<Hit> & hits ){
  //fill the reco hits object data product

  for(int l=0; l<tNumberOfHits; l++){

    Hit dummyHits;

    dummyHits.run=*fRun;
    dummyHits.subRun=tSubrun;
    dummyHits.event=tEventNumberInRun;
    dummyHits.TPC=tHit_TPC[l];
    dummyHits.view=tHit_View[l];
    dummyHits.channel=tHit_Channel[l];
    dummyHits.peakTime=tHit_PeakTime[l];
    dummyHits.chargeSummedADC=tHit_ChargeSummedADC[l];
    dummyHits.chargeIntegral=tHit_ChargeIntegral[l];
    dummyHits.startTime=tHit_StartTime[l];
    dummyHits.endTime=tHit_EndTime[l];
    dummyHits.width=tHit_Width[l];
    dummyHits.goodnessOfFit=tHit_GoodnessOfFit[l];
    dummyHits.multiplicity=tHit_Multiplicity[l];
    dummyHits.trackID=tHit_TrackID[l];
    dummyHits.particleID= tHit_TrackID[l];
    dummyHits.trueEnergy = tHit_TrueEnergy[l];
    dummyHits.trueEnergyFraction = tHit_TrueEnergyFraction[l];
    //dummyHit.lem=dummyHits.findLEM(dummyHits.Y, dummyHits.Z);

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

    dummyTrack.run = *fRun;
    dummyTrack.subRun=tSubrun;
    dummyTrack.event=tEventNumberInRun;
    dummyTrack.trackID=tTrackID[j];
    dummyTrack.numberOfHits=tTrack_NumberOfHits[j];
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

    for(int k=0; k<NUM_OF_VIEWS; k++) //View loop (for this track)
    {
      for(int l=a; l<a+tTrack_NumberOfHitsPerView[j][k]; l++) //Hit loop
      {
        Hit dummyHits;

        dummyHits.run = *fRun;
        dummyHits.subRun=tSubrun;
        dummyHits.event=tEventNumberInRun;
        dummyHits.X=tTrack_Hit_X[l];
        dummyHits.Y=tTrack_Hit_Y[l];
        dummyHits.Z=tTrack_Hit_X[l];
        dummyHits.dxLocalTrackDirection=tTrack_dx_LocalTrackDirection[l];
        dummyHits.dx3DPosition=tTrack_dx_3DPosition[l];
        dummyHits.TPC=tTrack_Hit_TPC[l];
        dummyHits.view=tTrack_Hit_View[l];
        dummyHits.channel=tTrack_Hit_Channel[l];
        dummyHits.peakTime=tTrack_Hit_PeakTime[l];
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
        //dummyHits.lem=dummyHits.findLEM(dummyHits.Y, dummyHits.Z);

        dummyTrack.hitsTrk.push_back(dummyHits);
      } //end hits
      a+=tTrack_NumberOfHitsPerView[j][k];;
    } //end view

  tracks.push_back(dummyTrack);
  }//end tracks

  return;
}

void LArParser::getMCTracksEvent( vector<MCTrack> & tracks, int event  ){
  //just fill the hit array for a specific event

  this->setMCBranches();

  fTree->GetEntry(event);
  this->fillMCTrack( tracks );

  return;

}

void LArParser::getRecoHitsEvent( vector<Hit> & hits, int event  ){
  //just fill the hit array for a specific event

  this->setRecoBranches();

  fTree->GetEntry(event);
  this->fillRecoHits( hits );

  return;

}

void LArParser::getRecoTracksEvent( vector<Track> & tracks, int event  ){
  //just fill the hit array for a specific event

  this->setRecoBranches();

  fTree->GetEntry(event);
  this->fillRecoTrack( tracks );

  return;

}

void LArParser::getMCTracks( vector<MCTrack> & tracks ){
  //return an array holding the reconstruced track

  this->setMCBranches();

  for(int i=0; i<fTree->GetEntries(); i++) //Event loop
  {
    fTree->GetEntry(i);
    this->fillMCTrack( tracks );
  }//Event loop

  return;
}

void LArParser::getRecoHits( vector<Hit> & hits ){
  //return a vector holding the free reconstructed hits

  this->setRecoBranches();

  for(int i=0; i<fTree->GetEntries(); i++){ //Event loop

    fTree->GetEntry(i);
    this->fillRecoHits( hits );
  }//Event loop

  return;
}

void LArParser::getRecoTracks( vector<Track> & tracks ){
  //return an array holding the reconstruced track

  this->setRecoBranches();

  for(int i=0; i<fTree->GetEntries(); i++) //Event loop
  {
    fTree->GetEntry(i);
    this->fillRecoTrack( tracks );
  }//Event loop

  return;
}
