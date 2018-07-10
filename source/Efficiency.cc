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

  //constructor of the class: initialize here the TTree
  fMcTTree = new TTree("mctree", "Truth info");
  fRecoTTree = new TTree("recotree", "Reco information");

  //set branches
  fMcTTree->Branch("FileNumber", &fFileNumber, "FileNumber/I");
  fMcTTree->Branch("Event", &fEvent, "Event/I");
  fMcTTree->Branch("ParticleId", &fParticleId, "ParticleId/I");
  fMcTTree->Branch("Pdg", &fPdg, "Pdg/D");
  fMcTTree->Branch("Theta", &fTrueTheta, "Theta/D");
  fMcTTree->Branch("Phi", &fTruePhi, "Phi/D");
  fMcTTree->Branch("Energy", &fTrueE, "Energy/D");

  fRecoTTree->Branch("FileNumber", &fFileNumber, "FileNumber/I");
  fRecoTTree->Branch("Event", &fEvent, "Event/I");
  fRecoTTree->Branch("TrackId", &fTrackId, "TrackId/I");
  fRecoTTree->Branch("NtracksEvent", &fNtracksEvent, "NtracksEvent/I");
  fRecoTTree->Branch("ParticleId", &fParticleId, "ParticleId/I");
  fRecoTTree->Branch("Pdg", &fPdg, "Pdg/D");
  fRecoTTree->Branch("Theta", &fRecoTheta, "Theta/D");
  fRecoTTree->Branch("Phi", &fRecoPhi, "Phi/D");
  fRecoTTree->Branch("Energy", &fRecoE, "Energy/D");
  fRecoTTree->Branch("Completeness", &fCompleteness, "Completeness/D");
  fRecoTTree->Branch("Purity", &fPurirty, "Purity/D");


  //histograms
  size_t arraySize = max(pdgCode.size(), pdgNames.size() );

  for( size_t i =0; i< arraySize; i++ ){

    TH1D* hThetaTrue = new TH1D(("hThetaTrue"+pdgNames.at(i)).c_str(), ("#theta (true) "+pdgNames.at(i)+";#theta (deg)").c_str() , nBinsTheta, thetaStart, thetaEnd);
    TH1D* hPhiTrue = new TH1D(("hPhiTrue"+pdgNames.at(i)).c_str(), ("#phi (true) "+pdgNames.at(i)+";#phi (deg)").c_str() , nBinsPhi, phiStart, phiEnd);
    TH2D* hPhiThetaTrue = new TH2D(("hPhiThetaTrue"+pdgNames.at(i)).c_str(), ("#phi vs. #theta (true) "+pdgNames.at(i)+"; #theta (deg);#phi (deg)").c_str(), nBinsTheta, thetaStart, thetaEnd, nBinsPhi, phiStart, phiEnd); ;

    TH1D* hThetaReco = new TH1D(("hThetaReco"+pdgNames.at(i)).c_str(), ("#theta (reco) "+pdgNames.at(i)+";#theta (deg)").c_str() , nBinsTheta, thetaStart, thetaEnd);;
    TH1D* hPhiReco = new TH1D(("hPhiReco"+pdgNames.at(i)).c_str(), ("#phi (reco) "+pdgNames.at(i)+";#phi (deg)").c_str() , nBinsPhi, phiStart, phiEnd);
    TH2D* hPhiThetaReco = new TH2D(("hPhiThetaReco"+pdgNames.at(i)).c_str(), ("#phi vs. #theta (true) "+pdgNames.at(i)+"; #theta (deg);#phi (deg)").c_str(), nBinsTheta, thetaStart, thetaEnd, nBinsPhi, phiStart, phiEnd);

    fThetaTrueMap[pdgCode.at(i)] = hThetaTrue;
    fThetaRecoMap[pdgCode.at(i)] = hThetaReco;
    fPhiTrueMap[pdgCode.at(i)] = hPhiTrue;
    fPhiRecoMap[pdgCode.at(i)] = hPhiReco;
    fPhiThetaTrueMap[pdgCode.at(i)] = hPhiThetaTrue;
    fPhiThetaRecoMap[pdgCode.at(i)] = hPhiThetaReco;

  }

}

Efficiency::~Efficiency(){}

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

void Efficiency::fillMap2D(int pdg, map<int, TH2D*> map, double fillX, double fillY ){
  //fill the map if the pdg code of the best tParticleId
  if( map.find(pdg) != map.end() )
    map[pdg]->Fill( fillX, fillY );
  else
    map[0]->Fill( fillX, fillY );
}

void Efficiency::fill(){
  //fill the reco quantitiess

  //match with truth
  this->matchTruth();

  fFileNumber = fParticleMap[ fBestTrackID ].run.getFileNumber();
  fEvent = fParticleMap[ fBestTrackID ].eventNumber;
  fParticleId = fBestTrackID;
  fTrackId = fTrack.trackID;
  fPdg = fParticleMap[ fBestTrackID ].pdgCode;
  fTrueTheta = fParticleMap[fBestTrackID].startTheta;
  fTruePhi = fParticleMap[fBestTrackID].startPhi;
  fTrueE = fParticleMap[fBestTrackID].startE;

  fRecoTheta = fParticleMap[fBestTrackID].startTheta;
  fRecoPhi = fParticleMap[fBestTrackID].startPhi;
  fRecoE = fParticleMap[fBestTrackID].startE;

  fMcTTree->Fill();
  fRecoTTree->Fill();

  //fill first the mc quanties
  fillMap1D( abs(fPdg), fThetaTrueMap, fParticleMap[fBestTrackID].startTheta);
  fillMap1D( abs(fPdg), fPhiTrueMap, fParticleMap[fBestTrackID].startPhi );
  fillMap2D( abs(fPdg), fPhiThetaTrueMap, fParticleMap[fBestTrackID].startTheta, fParticleMap[fBestTrackID].startPhi);

  if(fCompleteness>0.5 && fPurirty>0.5){

    //fill the reco quantities
    fillMap1D( abs(fPdg), fThetaRecoMap, fParticleMap[fBestTrackID].startTheta);
    fillMap1D( abs(fPdg), fPhiRecoMap, fParticleMap[fBestTrackID].startPhi );
    fillMap2D( abs(fPdg), fPhiThetaRecoMap, fParticleMap[fBestTrackID].startTheta, fParticleMap[fBestTrackID].startPhi);
  }

}


void Efficiency::makeEfficiencyPlot(){
  //calculate the efficieny plots

  size_t arraySize = max(pdgCode.size(), pdgNames.size() );

  for( size_t i =0; i< arraySize; i++ ){

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
