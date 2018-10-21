////////////////////////////////////////////////////////////////////////////////
// Reconstruction efficiency methods inplementation
//
//mailto:andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////

#include <glob.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h> //to use struct stat

#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"

#include "Run.hh"
#include "DataStructure.hh"
#include "TEfficiency.h"

#include "Efficiency.hh"

Efficiency::Efficiency(){
  this->initClass();
}

Efficiency::Efficiency( int fileNumber ){
  fFileNumber = fileNumber;
  this->initClass();
}

Efficiency::~Efficiency(){}

void Efficiency::initClass(){
  //constructor of the class: initialize here the TTree
  fMcTTree = new TTree("mctree", "Truth info");
  fRecoTTree = new TTree("recotree", "Reco information");

  //set branches
  fMcTTree->Branch("FileNumber", &fFileNumber, "FileNumber/I");
  fMcTTree->Branch("Event", &fEvent, "Event/I");
  fMcTTree->Branch("UniqueEvent", &fUniqueEventLabel, "UniqueEvent/I");
  fMcTTree->Branch("ParticleId", &fParticleId, "ParticleId/I");
  fMcTTree->Branch("Pdg", &fPdg, "Pdg/D");
  fMcTTree->Branch("Theta", &fTrueTheta, "Theta/D");
  fMcTTree->Branch("Phi", &fTruePhi, "Phi/D");
  fMcTTree->Branch("Energy", &fTrueE, "Energy/D");
  fMcTTree->Branch("TrueStartX", &fTrueStartX, "TrueStartX/D");
  fMcTTree->Branch("TrueStartY", &fTrueStartY, "TrueStartY/D");
  fMcTTree->Branch("TrueStartZ", &fTrueStartZ, "TrueStartZ/D");
  fMcTTree->Branch("TrueEndX", &fTrueEndX, "TrueEndX/D");
  fMcTTree->Branch("TrueEndY", &fTrueEndY, "TrueEndY/D");
  fMcTTree->Branch("TrueEndZ", &fTrueEndZ, "TrueEndZ/D");
  fMcTTree->Branch("TrueStartDirectionX", &fTrueStartDirectionX, "TrueStartDirectionX/D");
  fMcTTree->Branch("TrueStartDirectionY", &fTrueStartDirectionY, "TrueStartDirectionY/D");
  fMcTTree->Branch("TrueStartDirectionZ", &fTrueStartDirectionZ, "TrueStartDirectionZ/D");
  fMcTTree->Branch("TrueEndDirectionX", &fTrueEndDirectionX, "TrueEndDirectionX/D");
  fMcTTree->Branch("TrueEndDirectionY", &fTrueEndDirectionY, "TrueEndDirectionY/D");
  fMcTTree->Branch("TrueEndDirectionZ", &fTrueEndDirectionZ, "TrueEndDirectionZ/D");

  fRecoTTree->Branch("FileNumber", &fFileNumber, "FileNumber/I");
  fRecoTTree->Branch("Event", &fEvent, "Event/I");
  fRecoTTree->Branch("UniqueEvent", &fUniqueEventLabel, "UniqueEvent/I");
  fRecoTTree->Branch("TrackId", &fTrackId, "TrackId/I");
  fRecoTTree->Branch("NtracksEvent", &fNtracksEvent, "NtracksEvent/I");
  fRecoTTree->Branch("ParticleId", &fParticleId, "ParticleId/I");
  fRecoTTree->Branch("Pdg", &fPdg, "Pdg/D");
  fRecoTTree->Branch("Theta", &fRecoTheta, "Theta/D");
  fRecoTTree->Branch("Phi", &fRecoPhi, "Phi/D");
  fRecoTTree->Branch("Energy", &fRecoE, "Energy/D");
  fRecoTTree->Branch("Direction", &fDirection, "Direction/D");
  fRecoTTree->Branch("DiffStartX", &fDiffStartX, "DiffStartX/D");
  fRecoTTree->Branch("DiffStartY", &fDiffStartY, "DiffStartY/D");
  fRecoTTree->Branch("DiffStartZ", &fDiffStartZ, "DiffStartZ/D");
  fRecoTTree->Branch("StartX", &fStartX, "StartX/D");
  fRecoTTree->Branch("StartY", &fStartY, "StartY/D");
  fRecoTTree->Branch("StartZ", &fStartZ, "StartZ/D");
  fRecoTTree->Branch("EndX", &fEndX, "EndX/D");
  fRecoTTree->Branch("EndY", &fEndY, "EndY/D");
  fRecoTTree->Branch("EndZ", &fEndZ, "EndZ/D");
  fRecoTTree->Branch("StartDirectionX", &fStartDirectionX, "StartDirectionX/D");
  fRecoTTree->Branch("StartDirectionY", &fStartDirectionY, "StartDirectionY/D");
  fRecoTTree->Branch("StartDirectionZ", &fStartDirectionZ, "StartDirectionZ/D");
  fRecoTTree->Branch("EndDirectionX", &fEndDirectionX, "EndDirectionX/D");
  fRecoTTree->Branch("EndDirectionY", &fEndDirectionY, "EndDirectionY/D");
  fRecoTTree->Branch("EndDirectionZ", &fEndDirectionZ, "EndDirectionZ/D");
  fRecoTTree->Branch("RecoLength", &fRecoLength, "RecoLength/D");
  fRecoTTree->Branch("TrueLength", &fTrueLength, "TrueLength/D");
  fRecoTTree->Branch("Purity", &fPurirty, "Purity/D");
  fRecoTTree->Branch("Completeness", &fCompleteness, "Completeness/D");


  //histograms
  size_t arraySize = max(pdgCode.size(), pdgNames.size() );

  for( size_t i =0; i< arraySize; i++ ){

    TH1D* hThetaG4 = new TH1D(("hThetaG4"+pdgNames.at(i)).c_str(), ("#theta (geant) "+pdgNames.at(i)+";#theta (deg)").c_str() , nBinsTheta, thetaStart, thetaEnd);
    TH1D* hPhiG4 = new TH1D(("hPhiG4"+pdgNames.at(i)).c_str(), ("#phi (geant) "+pdgNames.at(i)+";#phi (deg)").c_str() , nBinsPhi, phiStart, phiEnd);
    //more geant info ?

    TH1D* hThetaTrue = new TH1D(("hThetaTrue"+pdgNames.at(i)).c_str(), ("#theta (true) "+pdgNames.at(i)+";#theta (deg)").c_str() , nBinsTheta, thetaStart, thetaEnd);
    TH1D* hPhiTrue = new TH1D(("hPhiTrue"+pdgNames.at(i)).c_str(), ("#phi (true) "+pdgNames.at(i)+";#phi (deg)").c_str() , nBinsPhi, phiStart, phiEnd);
    TH2D* hPhiThetaTrue = new TH2D(("hPhiThetaTrue"+pdgNames.at(i)).c_str(), ("#phi vs. #theta (true) "+pdgNames.at(i)+"; #theta (deg);#phi (deg)").c_str(), nBinsTheta, thetaStart, thetaEnd, nBinsPhi, phiStart, phiEnd); ;

    TH1D* hThetaReco = new TH1D(("hThetaReco"+pdgNames.at(i)).c_str(), ("#theta (reco) "+pdgNames.at(i)+";#theta (deg)").c_str() , nBinsTheta, thetaStart, thetaEnd);;
    TH1D* hPhiReco = new TH1D(("hPhiReco"+pdgNames.at(i)).c_str(), ("#phi (reco) "+pdgNames.at(i)+";#phi (deg)").c_str() , nBinsPhi, phiStart, phiEnd);
    TH2D* hPhiThetaReco = new TH2D(("hPhiThetaReco"+pdgNames.at(i)).c_str(), ("#phi vs. #theta (true) "+pdgNames.at(i)+"; #theta (deg);#phi (deg)").c_str(), nBinsTheta, thetaStart, thetaEnd, nBinsPhi, phiStart, phiEnd);

    TH1D *hDir = new TH1D(("hDir"+pdgNames.at(i)).c_str(), ("Direction "+pdgNames.at(i)+";true dir.dot( reco dir )").c_str() , nBinsDir, -1, 1);
    TH2D *hDirTheta = new TH2D(("hDirTheta"+pdgNames.at(i)).c_str(), ("Direction  vs #theta"+pdgNames.at(i)+";true dir.dot( reco dir ); #theta (deg)").c_str() , nBinsDir, -1, 1, nBinsTheta, thetaStart, thetaEnd);
    TH2D *hDirPhi = new TH2D(("hDirPhi"+pdgNames.at(i)).c_str(), ("Direction  vs #phi"+pdgNames.at(i)+";true dir.dot( reco dir ); #phi (deg)").c_str() , nBinsDir, -1, 1, nBinsPhi, phiStart, phiEnd);

    TH1D *hPosX = new TH1D(("hPosX"+pdgNames.at(i)).c_str(), ("Start position X"+pdgNames.at(i)+";true start - reco start (cm) )").c_str() , nBinsPos, -5, 5);
    TH2D *hPosXTheta = new TH2D(("hPosXTheta"+pdgNames.at(i)).c_str(), ("Start position X vs #theta"+pdgNames.at(i)+";true start - reco start (cm) ); #theta (deg)").c_str() , nBinsPos, -5, 5, nBinsTheta, thetaStart, thetaEnd);
    TH2D *hPosXPhi = new TH2D(("hPosXPhi"+pdgNames.at(i)).c_str(), ("Start position X vs #phi"+pdgNames.at(i)+";true start - reco start (cm) ); #phi (deg)").c_str() , nBinsPos, -5, 5, nBinsPhi, phiStart, phiEnd);

    TH1D *hPosY = new TH1D(("hPosY"+pdgNames.at(i)).c_str(), ("Start position Y"+pdgNames.at(i)+";true start - reco start (cm) )").c_str() , nBinsPos, -5, 5);
    TH2D *hPosYTheta = new TH2D(("hPosYTheta"+pdgNames.at(i)).c_str(), ("Start position Y vs #theta"+pdgNames.at(i)+";true start - reco start (cm) ); #theta (deg)").c_str() , nBinsPos, -5, 5, nBinsTheta, thetaStart, thetaEnd);
    TH2D *hPosYPhi = new TH2D(("hPosYPhi"+pdgNames.at(i)).c_str(), ("Start position Y vs #phi"+pdgNames.at(i)+";true start - reco start (cm) ); #phi (deg)").c_str() , nBinsPos, -5, 5, nBinsPhi, phiStart, phiEnd);

    TH1D *hPosZ = new TH1D(("hPosZ"+pdgNames.at(i)).c_str(), ("Start position Z"+pdgNames.at(i)+";true start - reco start (cm) )").c_str() , nBinsPos, -5, 5);
    TH2D *hPosZTheta = new TH2D(("hPosZTheta"+pdgNames.at(i)).c_str(), ("Start position Z vs #theta"+pdgNames.at(i)+";true start - reco start (cm) ); #theta (deg)").c_str() , nBinsPos, -5, 5, nBinsTheta, thetaStart, thetaEnd);
    TH2D *hPosZPhi = new TH2D(("hPosZPhi"+pdgNames.at(i)).c_str(), ("Start position Z vs #phi"+pdgNames.at(i)+";true start - reco start (cm) ); #phi (deg)").c_str() , nBinsPos, -5, 5, nBinsPhi, phiStart, phiEnd);

    fThetaG4Map[pdgCode.at(i)] = hThetaG4;
    fPhiG4Map[pdgCode.at(i)] = hPhiG4;
    fDirMap[pdgCode.at(i)] = hDir;
    fDirThetaMap[pdgCode.at(i)] = hDirTheta;
    fDirPhiMap[pdgCode.at(i)] = hDirPhi;
    fPosXMap[pdgCode.at(i)] = hPosX;
    fPosXThetaMap[pdgCode.at(i)] = hPosXTheta;
    fPosXPhiMap[pdgCode.at(i)] = hPosXPhi;
    fPosYMap[pdgCode.at(i)] = hPosY;
    fPosYThetaMap[pdgCode.at(i)] = hPosYTheta;
    fPosYPhiMap[pdgCode.at(i)] = hPosYPhi;
    fPosZMap[pdgCode.at(i)] = hPosZ;
    fPosZThetaMap[pdgCode.at(i)] = hPosZTheta;
    fPosZPhiMap[pdgCode.at(i)] = hPosZPhi;
    fThetaTrueMap[pdgCode.at(i)] = hThetaTrue;
    fThetaRecoMap[pdgCode.at(i)] = hThetaReco;
    fPhiTrueMap[pdgCode.at(i)] = hPhiTrue;
    fPhiRecoMap[pdgCode.at(i)] = hPhiReco;
    fPhiThetaTrueMap[pdgCode.at(i)] = hPhiThetaTrue;
    fPhiThetaRecoMap[pdgCode.at(i)] = hPhiThetaReco;

  }
}

void Efficiency::matchTruth(){
  //match reco and truth, calculate purity and completeness of a reco track

  //first of all reset all the variables
  fEnergyMap.clear();
  fPurirty=0;
  fCompleteness=0;

  //loop over the hits in tracks and associate every energy deposit to the correct particleID
  double energyTrk = 0.;
  for( auto hit : fTrack.hitsTrk ){
    fEnergyMap[ hit.particleID ] += hit.trueEnergy;
    energyTrk += (hit.trueEnergy/hit.trueEnergyFraction);
  }

  //find the best particle ID (the one that contribute the most in the track energy account )
  fBestTrackID =0.;
  double maxe = 0.;
  for( auto const & val : fEnergyMap  ){
    if( maxe < val.second ){
      maxe = val.second;
      fBestTrackID= val.first;
    }
  }

  //calculate the total energy of the best track in hits
  double totalEnergy = 0;
  for( auto hit : fHits ){
    if( fBestTrackID == hit.particleID )
      totalEnergy+=hit.trueEnergy;
  }

  fPurirty = (fEnergyMap[ fBestTrackID ])/totalEnergy;
  fCompleteness = fEnergyMap[ fBestTrackID ]/energyTrk;

}

void Efficiency::fillMap1D(int pdg, map<int, TH1D*> map, double fillIn ){
  //fill the map if the pdg code of the best tParticleId
    if( map.find(pdg) != map.end() )
      map[pdg]->Fill( fillIn );
    else
      map[0]->Fill( fillIn );
}

void Efficiency::fillMap2D(int pdg, map<int, TH2D*> map, double fillZ, double fillY ){
  //fill the map if the pdg code of the best tParticleId
  if( map.find(pdg) != map.end() )
    map[pdg]->Fill( fillZ, fillY );
  else
    map[0]->Fill( fillZ, fillY );
}

void Efficiency::setMapEntry(int id, MCTrack mctrack ){

  fParticleMap[id] = mctrack;

  //fill g4 histograms right afterwards
  fillMap1D( abs( fParticleMap[id].pdgCode ), fThetaG4Map, fParticleMap[id].startTheta );
  fillMap1D( abs( fParticleMap[id].pdgCode ), fPhiG4Map, fParticleMap[id].startPhi );

  //fFileNumber = fParticleMap[ id ].run;

  fEvent = fParticleMap[ id ].eventNumber; //get to the event a progressive number across files
  fUniqueEventLabel = fFileNumber*100 + fEvent;
  fParticleId = id;
  fPdg = fParticleMap[ id ].pdgCode;
  fTrueTheta = fParticleMap[ id ].startTheta;
  fTruePhi = fParticleMap[ id ].startPhi;
  fTrueE = fParticleMap[ id ].startE;
  fTrueStartX = fParticleMap[ id ].startX;
  fTrueStartY = fParticleMap[ id ].startY;
  fTrueStartZ = fParticleMap[ id ].startZ;
  fTrueEndX = fParticleMap[ id ].endX;
  fTrueEndY = fParticleMap[ id ].endY;
  fTrueEndZ = fParticleMap[ id ].endZ;

  fTrueStartDirectionX = fParticleMap[ id ].startDirection.X();
  fTrueStartDirectionY = fParticleMap[ id ].startDirection.Y();
  fTrueStartDirectionZ = fParticleMap[ id ].startDirection.Z();
  fTrueEndDirectionX = fParticleMap[ id ].endDirection.X();
  fTrueEndDirectionY = fParticleMap[ id ].endDirection.Y();
  fTrueEndDirectionZ = fParticleMap[ id ].endDirection.Z();

  fMcTTree->Fill();

}

void Efficiency::fill(){
  //fill the histograms

  //match with truth
  this->matchTruth();

  fEvent = fParticleMap[ fBestTrackID ].eventNumber; //get to the event a progressive number across files
  fUniqueEventLabel = fFileNumber*100 + fEvent; //assumung 100 events per file

  fParticleId = fBestTrackID;
  fTrackId = fTrack.trackID;
  fPdg = fParticleMap[ fBestTrackID ].pdgCode;
  fRecoTheta = fParticleMap[fBestTrackID].startTheta;
  fRecoPhi = fParticleMap[fBestTrackID].startPhi;
  fRecoE = fParticleMap[fBestTrackID].startE;

  fDirection = fTrack.startDirection.Dot( fParticleMap[fBestTrackID].startDirection );
  fStartDirectionX = fTrack.startDirection.X();
  fStartDirectionY = fTrack.startDirection.Y();
  fStartDirectionZ = fTrack.startDirection.Z();
  fEndDirectionX = fTrack.endDirection.X();
  fEndDirectionY = fTrack.endDirection.Y();
  fEndDirectionZ = fTrack.endDirection.Z();

  fDiffStartX = min(fTrack.startPointX - fParticleMap[fBestTrackID].startX, fTrack.endPointX - fParticleMap[fBestTrackID].startX);
  fDiffStartY = min(fTrack.startPointY - fParticleMap[fBestTrackID].startY, fTrack.endPointY - fParticleMap[fBestTrackID].startY);
  fDiffStartZ = min(fTrack.startPointZ - fParticleMap[fBestTrackID].startZ, fTrack.endPointZ - fParticleMap[fBestTrackID].startZ);

  fStartX = fTrack.startPointX;
  fStartY = fTrack.startPointY;
  fStartZ = fTrack.startPointZ;

  fEndX = fTrack.endPointX;
  fEndY = fTrack.endPointY;
  fEndZ = fTrack.endPointZ;

  fRecoLength = fTrack.length;
  fTrueLength = fParticleMap[fBestTrackID].length;



  fRecoTTree->Fill();

  //fill first the mc quanties
  fillMap1D( abs(fPdg), fDirMap, fDirection);
  fillMap2D( abs(fPdg), fDirThetaMap, fDirection, fParticleMap[fBestTrackID].startTheta);
  fillMap2D( abs(fPdg), fDirPhiMap, fDirection, fParticleMap[fBestTrackID].startPhi);

  fillMap1D( abs(fPdg), fPosXMap, fDiffStartX);
  fillMap2D( abs(fPdg), fPosXThetaMap, fDiffStartX, fParticleMap[fBestTrackID].startTheta);
  fillMap2D( abs(fPdg), fPosXPhiMap, fDiffStartX, fParticleMap[fBestTrackID].startPhi);

  fillMap1D( abs(fPdg), fPosYMap, fDiffStartY);
  fillMap2D( abs(fPdg), fPosYThetaMap, fDiffStartY, fParticleMap[fBestTrackID].startTheta);
  fillMap2D( abs(fPdg), fPosYPhiMap, fDiffStartY, fParticleMap[fBestTrackID].startPhi);

  fillMap1D( abs(fPdg), fPosZMap, fDiffStartZ);
  fillMap2D( abs(fPdg), fPosZThetaMap, fDiffStartZ, fParticleMap[fBestTrackID].startTheta);
  fillMap2D( abs(fPdg), fPosZPhiMap, fDiffStartZ, fParticleMap[fBestTrackID].startPhi);

  fillMap1D( abs(fPdg), fThetaTrueMap, fParticleMap[fBestTrackID].startTheta);
  fillMap1D( abs(fPdg), fPhiTrueMap, fParticleMap[fBestTrackID].startPhi );
  fillMap2D( abs(fPdg), fPhiThetaTrueMap, fParticleMap[fBestTrackID].startTheta, fParticleMap[fBestTrackID].startPhi);

  if(fCompleteness>0.5 && fPurirty>0.5){

    //fill the reco quantities
    this->fillMap1D( abs(fPdg), fThetaRecoMap, fParticleMap[fBestTrackID].startTheta);
    this->fillMap1D( abs(fPdg), fPhiRecoMap, fParticleMap[fBestTrackID].startPhi );
    this->fillMap2D( abs(fPdg), fPhiThetaRecoMap, fParticleMap[fBestTrackID].startTheta, fParticleMap[fBestTrackID].startPhi);
  }

}


void Efficiency::makeEfficiencyPlot(){
  //calculate the efficieny plots

  size_t arraySize = max(pdgCode.size(), pdgNames.size() );

  for( size_t i =0; i< arraySize; i++ ){

    //Only for matched tracks
    if( TEfficiency::CheckConsistency(*fPhiRecoMap[pdgCode.at(i)], *fPhiTrueMap[pdgCode.at(i)]) )
      fPhiEfficiency[pdgCode.at(i)] = new TEfficiency(*fPhiRecoMap[pdgCode.at(i)], *fPhiTrueMap[pdgCode.at(i)]);

    if( TEfficiency::CheckConsistency(*fThetaRecoMap[pdgCode.at(i)], *fThetaTrueMap[pdgCode.at(i)]) )
      fThetaEfficiency[pdgCode.at(i)] = new TEfficiency(*fThetaRecoMap[pdgCode.at(i)], *fThetaTrueMap[pdgCode.at(i)]);

    if( TEfficiency::CheckConsistency(*fPhiThetaRecoMap[pdgCode.at(i)], *fPhiThetaTrueMap[pdgCode.at(i)]) )
      fPhiThetaEfficiency[pdgCode.at(i)] = new TEfficiency(*fPhiThetaRecoMap[pdgCode.at(i)], *fPhiThetaTrueMap[pdgCode.at(i)]);

    fPhiEfficiency[pdgCode.at(i)]->SetName( (pdgNames.at(i)+"Phi_Efficiency").c_str() );
    fThetaEfficiency[pdgCode.at(i)]->SetName( (pdgNames.at(i)+"Theta_Efficiency").c_str() );
    fPhiThetaEfficiency[pdgCode.at(i)]->SetName( (pdgNames.at(i)+"PhiTheta_Efficiency").c_str() );
    fPhiEfficiency[pdgCode.at(i)]->SetTitle( (pdgNames.at(i)+"Phi_Efficiency").c_str() );
    fThetaEfficiency[pdgCode.at(i)]->SetTitle( (pdgNames.at(i)+"Theta_Efficiency").c_str() );
    fPhiThetaEfficiency[pdgCode.at(i)]->SetTitle( (pdgNames.at(i)+"PhiTheta_Efficiency").c_str() );

  }
}

void Efficiency::write(){

  size_t arraySize = max(pdgCode.size(), pdgNames.size() );

  for( size_t i =0; i< arraySize; i++ ){

    fThetaG4Map[pdgCode.at(i)]->Write();
    fPhiG4Map[pdgCode.at(i)]->Write();
    fDirMap[pdgCode.at(i)]->Write();
    fDirThetaMap[pdgCode.at(i)]->Write();
    fDirPhiMap[pdgCode.at(i)]->Write();
    fPosXMap[pdgCode.at(i)]->Write();
    fPosXThetaMap[pdgCode.at(i)]->Write();
    fPosXPhiMap[pdgCode.at(i)]->Write();
    fPosYMap[pdgCode.at(i)]->Write();
    fPosYThetaMap[pdgCode.at(i)]->Write();
    fPosYPhiMap[pdgCode.at(i)]->Write();
    fPosZMap[pdgCode.at(i)]->Write();
    fPosZThetaMap[pdgCode.at(i)]->Write();
    fPosZPhiMap[pdgCode.at(i)]->Write();
    fThetaTrueMap[pdgCode.at(i)]->Write();
    fThetaRecoMap[pdgCode.at(i)]->Write();
    fPhiTrueMap[pdgCode.at(i)]->Write();
    fPhiRecoMap[pdgCode.at(i)]->Write();
    fPhiThetaTrueMap[pdgCode.at(i)]->Write();
    fPhiThetaRecoMap[pdgCode.at(i)]->Write();
    fThetaEfficiency[pdgCode.at(i)]->Write();
    fPhiEfficiency[pdgCode.at(i)]->Write();
    fPhiThetaEfficiency[pdgCode.at(i)]->Write();

  }

  //write tree
  fMcTTree->Write();
  fRecoTTree->Write();
}

void Efficiency::clean(){

  //quantities to be reset after the event
  fParticleMap.clear();
  fHits.clear();
}
