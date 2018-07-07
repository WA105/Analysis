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
  fMcTTree->Branch("Pdg", &fPdg, "Pdg/D");
  fMcTTree->Branch("Theta", &fTrueTheta, "Theta/D");
  fMcTTree->Branch("Phi", &fTruePhi, "Phi/D");
  fMcTTree->Branch("Energy", &fTrueE, "Energy/D");

  fRecoTTree->Branch("Pdg", &fPdg, "Pdg/D");
  fRecoTTree->Branch("Theta", &fRecoTheta, "Theta/D");
  fRecoTTree->Branch("Phi", &fRecoPhi, "Phi/D");
  fRecoTTree->Branch("Energy", &fRecoE, "Energy/D");
  fRecoTTree->Branch("Completeness", &fCompleteness, "Completeness/D");
  fRecoTTree->Branch("Purity", &fPurirty, "Purity/D");

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

void Efficiency::fill(){
  //fill the reco quantitiess

  //match with truth
  this->matchTruth();

  fPdg = fParticleMap[ fBestTrackID ].pdgCode;
  fTrueTheta = fParticleMap[fBestTrackID].startTheta;
  fTruePhi = fParticleMap[fBestTrackID].startPhi;
  fTrueE = fParticleMap[fBestTrackID].startE;

  fRecoTheta = fParticleMap[fBestTrackID].startTheta;
  fRecoPhi = fParticleMap[fBestTrackID].startPhi;
  fRecoE = fParticleMap[fBestTrackID].startE;

  fMcTTree->Fill();
  fRecoTTree->Fill();

}

void Efficiency::write(){
  //write tree
  fMcTTree->Write();
  fRecoTTree->Write();
}

void Efficiency::clean(){

  fParticleMap.clear();
  fHits.clear();
}
