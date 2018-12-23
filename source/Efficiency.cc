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

//==============================================================================

HitEfficiency::HitEfficiency( map< int, vector<Hit> > hitsMap  )
{
  fHitsMap = hitsMap;
}

//------------------------------------------------------------------------------

HitEfficiency::~HitEfficiency()
{
  clean();
}

//------------------------------------------------------------------------------

void HitEfficiency::clean()
{
  fHitsMap.clear();
}

//------------------------------------------------------------------------------

int HitEfficiency::getNumOfHits( int id, int view )
{
  int numOfHits = 0;

  auto lambda = [ view ]( int accumulator, Hit h )
  {
     return h.view == view ? accumulator + 1 : accumulator;
  };

  numOfHits = accumulate( fHitsMap[id].begin(), fHitsMap[id].end(), 0, lambda );

  return numOfHits;
}

//------------------------------------------------------------------------------

double HitEfficiency::getAvgCompleteness( int id, int view )
{
  double sum = 0.;

  auto lambda = [ view, &sum ]( double accumulator, Hit hit )
  {
     if( hit.view == view && hit.chanMaxElectrons > 0 )
     {
       ++sum;
       return accumulator + hit.electronsFromSumADC("Montecarlo")/hit.chanMaxElectrons;
     }
     else
     {
       return accumulator;
     }
  };

  double accumulated = accumulate( fHitsMap[id].begin(), fHitsMap[id].end(), 0.0, lambda );

  return sum > 0. ? accumulated/sum : 0.;
}

//------------------------------------------------------------------------------

double HitEfficiency::getAvgPurity( int id, int view )
{
  double sum = 0.;

  auto lambda = [ view, &sum ]( double accumulator, Hit hit )
  {
     if( hit.view == view && hit.chanElectrons > 0 )
     {
       ++sum;
       return accumulator + hit.electronsFromSumADC("Montecarlo")/hit.chanElectrons;
     }
     else
     {
       return accumulator;
     }
  };

  double accumulated = accumulate( fHitsMap[id].begin(), fHitsMap[id].end(), 0.0, lambda );

  return sum > 0. ? accumulated/sum : 0.;
}

//------------------------------------------------------------------------------

double HitEfficiency::getAvgSignalToNoise( int id, int view )
{
  double sum = 0.;

  auto lambda = [ view, &sum ]( double accumulator, Hit hit )
  {
     if( hit.view == view && hit.getPeakToNoise() > 0 )
     {
       ++sum;
       return accumulator + hit.getPeakToNoise();
     }
     else
     {
       return accumulator;
     }
  };

  double accumulated = accumulate( fHitsMap[id].begin(), fHitsMap[id].end(), 0.0, lambda );

  return sum > 0. ? accumulated/sum : 0.;
}

//------------------------------------------------------------------------------

double HitEfficiency::getAvgIntegralToNoise( int id, int view )
{
  double sum = 0.;

  auto lambda = [ view, &sum ]( double accumulator, Hit hit )
  {
     if( hit.view == view && hit.getIntegralToNoise("Montecarlo") > 0 )
     {
       ++sum;
       return accumulator + hit.getIntegralToNoise("Montecarlo");
     }
     else
     {
       return accumulator;
     }
  };

  double accumulated = accumulate( fHitsMap[id].begin(), fHitsMap[id].end(), 0.0, lambda );

  return sum > 0. ? accumulated/sum : 0.;
}

//------------------------------------------------------------------------------

double HitEfficiency::getAvgIntegral( int id, int view )
{
  double sum = 0.;

  auto lambda = [ view, &sum ]( double accumulator, Hit hit )
  {
     if( hit.view == view && hit.chargeIntegral > 0 )
     {
       ++sum;
       return accumulator + hit.chargeIntegral;
     }
     else
     {
       return accumulator;
     }
  };

  double accumulated = accumulate( fHitsMap[id].begin(), fHitsMap[id].end(), 0.0, lambda );

  return sum > 0. ? accumulated/sum : 0.;
}

//------------------------------------------------------------------------------

double HitEfficiency::getAvgAmplitude( int id, int view )
{
  double sum = 0.;

  auto lambda = [ view, &sum ]( double accumulator, Hit hit )
  {
     if( hit.view == view && hit.amplitude > 0 )
     {
       ++sum;
       return accumulator + hit.amplitude;
     }
     else
     {
       return accumulator;
     }
  };

  double accumulated = accumulate( fHitsMap[id].begin(), fHitsMap[id].end(), 0.0, lambda );

  return sum > 0. ? accumulated/sum : 0.;
}

//------------------------------------------------------------------------------

double HitEfficiency::getAvgWidth( int id, int view )
{
  double sum = 0.;

  auto lambda = [ view, &sum ]( double accumulator, Hit hit )
  {
     if( hit.view == view && hit.width > 0 )
     {
       ++sum;
       return accumulator + hit.width;
     }
     else
     {
       return accumulator;
     }
  };

  double accumulated = accumulate( fHitsMap[id].begin(), fHitsMap[id].end(), 0.0, lambda );

  return sum > 0. ? accumulated/sum : 0.;
}

//==============================================================================

Efficiency::Efficiency()
{
  this->initClass();
}

Efficiency::Efficiency( int fileNumber )
{
  fFileNumber = fileNumber;
  this->initClass();
}

Efficiency::~Efficiency(){}

void Efficiency::initClass()
{
  //constructor of the class: initialize here the TTree
  fMcTTree = new TTree("mctree", "Truth info");
  fRecoTTree = new TTree("recotree", "Reco information");
  fUnmatchTTree = new TTree("unmatchtree", "True track that didn't match with a reco track");

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
  fMcTTree->Branch("TrueLength", &fTrueLength, "TrueLength/D");
  fMcTTree->Branch("NumberOfHits", &fNumberOfHits, "NumberOfHits/I");
  fMcTTree->Branch("NumberOfHitsView0", &fNumberOfHitsView0, "NumberOfHitsView0/I");
  fMcTTree->Branch("NumberOfHitsView1", &fNumberOfHitsView1, "NumberOfHitsView1/I");
  fMcTTree->Branch("AvgCompletenessView0", &fAvgCompletenessView0, "AvgCompletenessView0/D");
  fMcTTree->Branch("AvgCompletenessView1", &fAvgCompletenessView1, "AvgCompletenessView1/D");
  fMcTTree->Branch("AvgPurityView0", &fAvgPurityView0, "AvgPurityView0/D");
  fMcTTree->Branch("AvgPurityView1", &fAvgPurityView1, "AvgPurityView1/D");
  fMcTTree->Branch("AvgSignalToNoiseView0", &fAvgSignalToNoiseView0, "AvgSignalToNoiseView0/D");
  fMcTTree->Branch("AvgSignalToNoiseView1", &fAvgSignalToNoiseView1, "AvgSignalToNoiseView1/D");
  fMcTTree->Branch("AvgIntegralToNoiseView0", &fAvgIntegralToNoiseView0, "AvgIntegralToNoiseView0/D");
  fMcTTree->Branch("AvgIntegralToNoiseView1", &fAvgIntegralToNoiseView1, "AvgIntegralToNoiseView1/D");
  fMcTTree->Branch("AvgIntegralView0", &fAvgIntegralView0, "AvgIntegralView0/D");
  fMcTTree->Branch("AvgIntegralView1", &fAvgIntegralView1, "AvgIntegralView1/D");
  fMcTTree->Branch("AvgIntegralView0", &fAvgAmplitudeView0, "AvgAmplitudeView0/D");
  fMcTTree->Branch("AvgAmplitudeView1", &fAvgAmplitudeView1, "AvgAmplitudeView1/D");
  fMcTTree->Branch("AvgIntegralView0", &fAvgWidthView0, "AvgWidthView0/D");
  fMcTTree->Branch("AvgWidthView1", &fAvgWidthView1, "AvgWidthView1/D");


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
  fRecoTTree->Branch("TrueStartX", &fTrueStartX, "TrueStartX/D");
  fRecoTTree->Branch("TrueStartY", &fTrueStartY, "TrueStartY/D");
  fRecoTTree->Branch("TrueStartZ", &fTrueStartZ, "TrueStartZ/D");
  fRecoTTree->Branch("TrueEndX", &fTrueEndX, "TrueEndX/D");
  fRecoTTree->Branch("TrueEndY", &fTrueEndY, "TrueEndY/D");
  fRecoTTree->Branch("TrueEndZ", &fTrueEndZ, "TrueEndZ/D");
  fRecoTTree->Branch("TrueStartDirectionX", &fTrueStartDirectionX, "TrueStartDirectionX/D");
  fRecoTTree->Branch("TrueStartDirectionY", &fTrueStartDirectionY, "TrueStartDirectionY/D");
  fRecoTTree->Branch("TrueStartDirectionZ", &fTrueStartDirectionZ, "TrueStartDirectionZ/D");
  fRecoTTree->Branch("TrueEndDirectionX", &fTrueEndDirectionX, "TrueEndDirectionX/D");
  fRecoTTree->Branch("TrueEndDirectionY", &fTrueEndDirectionY, "TrueEndDirectionY/D");
  fRecoTTree->Branch("TrueEndDirectionZ", &fTrueEndDirectionZ, "TrueEndDirectionZ/D");
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
  fRecoTTree->Branch("NumberOfHits", &fNumberOfHits, "NumberOfHits/I");
  fRecoTTree->Branch("NumberOfHitsView0", &fNumberOfHitsView0, "NumberOfHitsView0/I");
  fRecoTTree->Branch("NumberOfHitsView1", &fNumberOfHitsView1, "NumberOfHitsView1/I");
  fRecoTTree->Branch("AvgSignalToNoiseView0", &fAvgSignalToNoiseView0, "AvgSignalToNoiseView0/D");
  fRecoTTree->Branch("AvgSignalToNoiseView1", &fAvgSignalToNoiseView1, "AvgSignalToNoiseView1/D");
  fRecoTTree->Branch("AvgIntegralToNoiseView0", &fAvgIntegralToNoiseView0, "AvgIntegralToNoiseView0/D");
  fRecoTTree->Branch("AvgIntegralToNoiseView1", &fAvgIntegralToNoiseView1, "AvgIntegralToNoiseView1/D");
  fRecoTTree->Branch("AvgIntegralView0", &fAvgIntegralView0, "AvgIntegralView0/D");
  fRecoTTree->Branch("AvgIntegralView1", &fAvgIntegralView1, "AvgIntegralView1/D");
  fRecoTTree->Branch("AvgIntegralView0", &fAvgAmplitudeView0, "AvgAmplitudeView0/D");
  fRecoTTree->Branch("AvgAmplitudeView1", &fAvgAmplitudeView1, "AvgAmplitudeView1/D");
  fRecoTTree->Branch("AvgIntegralView0", &fAvgWidthView0, "AvgWidthView0/D");
  fRecoTTree->Branch("AvgWidthView1", &fAvgWidthView1, "AvgWidthView1/D");

  fUnmatchTTree->Branch("FileNumber", &fFileNumber, "FileNumber/I");
  fUnmatchTTree->Branch("Event", &fEvent, "Event/I");
  fUnmatchTTree->Branch("UniqueEvent", &fUniqueEventLabel, "UniqueEvent/I");
  fUnmatchTTree->Branch("ParticleId", &fParticleId, "ParticleId/I");
  fUnmatchTTree->Branch("Pdg", &fPdg, "Pdg/D");
  fUnmatchTTree->Branch("Theta", &fTrueTheta, "Theta/D");
  fUnmatchTTree->Branch("Phi", &fTruePhi, "Phi/D");
  fUnmatchTTree->Branch("Energy", &fTrueE, "Energy/D");
  fUnmatchTTree->Branch("TrueStartX", &fTrueStartX, "TrueStartX/D");
  fUnmatchTTree->Branch("TrueStartY", &fTrueStartY, "TrueStartY/D");
  fUnmatchTTree->Branch("TrueStartZ", &fTrueStartZ, "TrueStartZ/D");
  fUnmatchTTree->Branch("TrueEndX", &fTrueEndX, "TrueEndX/D");
  fUnmatchTTree->Branch("TrueEndY", &fTrueEndY, "TrueEndY/D");
  fUnmatchTTree->Branch("TrueEndZ", &fTrueEndZ, "TrueEndZ/D");
  fUnmatchTTree->Branch("TrueStartDirectionX", &fTrueStartDirectionX, "TrueStartDirectionX/D");
  fUnmatchTTree->Branch("TrueStartDirectionY", &fTrueStartDirectionY, "TrueStartDirectionY/D");
  fUnmatchTTree->Branch("TrueStartDirectionZ", &fTrueStartDirectionZ, "TrueStartDirectionZ/D");
  fUnmatchTTree->Branch("TrueEndDirectionX", &fTrueEndDirectionX, "TrueEndDirectionX/D");
  fUnmatchTTree->Branch("TrueEndDirectionY", &fTrueEndDirectionY, "TrueEndDirectionY/D");
  fUnmatchTTree->Branch("TrueEndDirectionZ", &fTrueEndDirectionZ, "TrueEndDirectionZ/D");
  fUnmatchTTree->Branch("TrueLength", &fTrueLength, "TrueLength/D");
  fUnmatchTTree->Branch("NumberOfHits", &fNumberOfHits, "NumberOfHits/I");
  fUnmatchTTree->Branch("NumberOfHitsView0", &fNumberOfHitsView0, "NumberOfHitsView0/I");
  fUnmatchTTree->Branch("NumberOfHitsView1", &fNumberOfHitsView1, "NumberOfHitsView1/I");

  fHitTree = new TTree("hits", "All hits information");

  fHitTree->Branch("FileNumber", &fFileNumber, "FileNumber/I");
  fHitTree->Branch("Event" , &fEvent, "Event/I");
  fHitTree->Branch("UniqueEvent", &fUniqueEventLabel, "UniqueEvent/I");
  fHitTree->Branch("View" , &fView, "View/I");
  fHitTree->Branch("Channel" , &fWire, "Channel/D");
  fHitTree->Branch("PeakTime" , &fPeakTime, "PeakTime/D");
  fHitTree->Branch("StartTime" , &fStartTime, "StartTime/D");
  fHitTree->Branch("EndTime" , &fEndTime, "EndTime/D");
  fHitTree->Branch("Width" , &fWidth, "Width/D");
  fHitTree->Branch("Amplitude" , &fAmplitude, "Amplitude/D");
  fHitTree->Branch("SummedADC" , &fSummedADC, "SummedADC/D");
  fHitTree->Branch("Integral" , &fIntegral, "Integral/D");
  fHitTree->Branch("GoodnessOfFit" , &fGoodnessOfFit, "GoodnessOfFit/D");
  fHitTree->Branch("Multiplicity" , &fMultiplicity, "Multiplicity/D");
  fHitTree->Branch("HitParticleId", &fHitParticleId, "HitParticleId/I");
  fHitTree->Branch("HitPurity", &fHitPurity, "HitPurity/D");
  fHitTree->Branch("HitCompleteness", &fHitCompleteness, "HitCompleteness/D");
  fHitTree->Branch("fIDERatio", &fIDERatio, "IDERatio/D");
  fHitTree->Branch("fIntegralToNoise", &fIntegralToNoise, "IntegralToNoise/D");
  fHitTree->Branch("fPeakToNoise", &fPeakToNoise, "PeakToNoiseToNoise/D");
  fHitTree->Branch("Pdg", &fPdg, "Pdg/D");
  fHitTree->Branch("Theta", &fTrueTheta, "Theta/D");
  fHitTree->Branch("Phi", &fTruePhi, "Phi/D");
  fHitTree->Branch("TrueLength", &fTrueLength, "TrueLength/D");

  //histograms
  size_t arraySize = max(pdgCode.size(), pdgNames.size() );

  for( size_t i =0; i< arraySize; i++ )
  {
    TH1D* hThetaG4 = new TH1D(("hThetaG4"+pdgNames.at(i)).c_str(), ("#theta (geant) "+pdgNames.at(i)+";#theta (deg)").c_str() , nBinsTheta, thetaStart, thetaEnd);
    TH1D* hPhiG4 = new TH1D(("hPhiG4"+pdgNames.at(i)).c_str(), ("#phi (geant) "+pdgNames.at(i)+";#phi (deg)").c_str() , nBinsPhi, phiStart, phiEnd);

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

void Efficiency::setRecoHits( vector<Hit> hits )
{
  //Handles the hit object, create map and fill hits ttree
  fHits = hits;

  for(Hit hit : hits ) { fHitsMap[ hit.particleID ].push_back(hit); }
  fHitEfficiency = new HitEfficiency( fHitsMap );
}

void Efficiency::fillHitTree( int id )
{

  for( auto hit : fHitsMap[id] )
  {
    fEvent = hit.event;
    fUniqueEventLabel = fFileNumber*100 + fEvent;
    fView = hit.view;
    fWire = hit.channel;
    fPeakTime = hit.peakTime;
    fAmplitude = hit.amplitude;
    fSummedADC = hit.chargeSummedADC;
    fIntegral = hit.chargeIntegral;
    fStartTime = hit.startTime;
    fEndTime = hit.endTime;
    fWidth = hit.width;
    fGoodnessOfFit = hit.goodnessOfFit;
    fMultiplicity = hit.multiplicity;

    //backtracker quantites
    fHitParticleId = hit.particleID;
    fHitPurity = hit.electronsFromSumADC("Montecarlo")/hit.chanElectrons;
    fHitCompleteness = hit.electronsFromSumADC("Montecarlo")/hit.chanMaxElectrons;
    fIDERatio = hit.chanMaxElectrons/hit.chanElectrons;

    //signal-to-noise
    fIntegralToNoise = hit.getSummedADCToNoise();
    fPeakToNoise = hit.getPeakToNoise();

    //associated true particle info
    fPdg = fParticleMap[ id ].pdgCode;
    fTrueTheta = fParticleMap[ id ].startTheta;
    fTruePhi = fParticleMap[ id ].startPhi;
    fTrueLength = fParticleMap[ id ].length;

    fHitTree->Fill();
  }
}

void Efficiency::fillMap1D(int pdg, map<int, TH1D*> map, double fillIn )
{
  //fill the map if the pdg code of the best tParticleId
    if( map.find(pdg) != map.end() )
      map[pdg]->Fill( fillIn );
    else
      map[0]->Fill( fillIn );
}

void Efficiency::fillMap2D(int pdg, map<int, TH2D*> map, double fillZ, double fillY )
{
  //fill the map if the pdg code of the best tParticleId
  if( map.find(pdg) != map.end() )
    map[pdg]->Fill( fillZ, fillY );
  else
    map[0]->Fill( fillZ, fillY );
}

void Efficiency::setMapEntry(int id, MCTrack mctrack ){

  fParticleMap[id] = mctrack;

  //fill g4 histograms and ttree
  fillMap1D( abs( fParticleMap[id].pdgCode ), fThetaG4Map, fParticleMap[id].startTheta );
  fillMap1D( abs( fParticleMap[id].pdgCode ), fPhiG4Map, fParticleMap[id].startPhi );

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
  fTrueLength = fParticleMap[ id ].length;
  fNumberOfHitsView0 = fHitEfficiency->getNumOfHits( id, 0 );
  fNumberOfHitsView1 = fHitEfficiency->getNumOfHits( id, 1 );
  fNumberOfHits = fNumberOfHitsView0+fNumberOfHitsView1;
  fAvgCompletenessView0 = fHitEfficiency->getAvgCompleteness( id, 0 );
  fAvgCompletenessView1 = fHitEfficiency->getAvgCompleteness( id, 1 );
  fAvgPurityView0 = fHitEfficiency->getAvgPurity( id, 0 );
  fAvgPurityView1 = fHitEfficiency->getAvgPurity( id, 1 );
  fAvgSignalToNoiseView0 = fHitEfficiency->getAvgSignalToNoise( id, 0 );
  fAvgSignalToNoiseView1 = fHitEfficiency->getAvgSignalToNoise( id, 1 );
  fAvgIntegralToNoiseView0 = fHitEfficiency->getAvgIntegralToNoise( id, 0 );
  fAvgIntegralToNoiseView1 = fHitEfficiency->getAvgIntegralToNoise( id, 1 );
  fAvgIntegralView0 = fHitEfficiency->getAvgIntegral( id, 0 );
  fAvgIntegralView1 = fHitEfficiency->getAvgIntegral( id, 1 );
  fAvgAmplitudeView0 = fHitEfficiency->getAvgAmplitude( id, 0 );
  fAvgAmplitudeView1 = fHitEfficiency->getAvgAmplitude( id, 1 );
  fAvgWidthView0 = fHitEfficiency->getAvgWidth( id, 0 );
  fAvgWidthView1 = fHitEfficiency->getAvgWidth( id, 1 );

  fMcTTree->Fill();

  if(fMakeHitsTree){ this->fillHitTree(id); }

}

void Efficiency::fill()
{
  //fill the histograms

  fBestTrackID = fTrack.bestParticleID;
  fCompleteness = fTrack.completeness;
  fPurirty = fTrack.purity;

  //check if fBestTrackID belongs to fParticleMap
  if ( fParticleMap.find( fBestTrackID ) == fParticleMap.end() )
  {
    std::cout << "WARNING! Efficiency::fill(): key not found in fParticleMap" << std::endl;
    return; //don't fill the ttree
  }

  //if the key exists, fill up all the quantities
  fEvent = fParticleMap[ fBestTrackID ].eventNumber; //get to the event a progressive number across files
  fUniqueEventLabel = fFileNumber*100 + fEvent; //assumung 100 events per file

  fParticleId = fBestTrackID;
  fTrackId = fTrack.trackID;
  fPdg = fParticleMap[ fBestTrackID ].pdgCode;
  fRecoTheta = fParticleMap[fBestTrackID].startTheta;
  fRecoPhi = fParticleMap[fBestTrackID].startPhi;
  fRecoE = fParticleMap[fBestTrackID].startE;
  fTrueStartX = fParticleMap[ fBestTrackID ].startX;
  fTrueStartY = fParticleMap[ fBestTrackID ].startY;
  fTrueStartZ = fParticleMap[ fBestTrackID ].startZ;
  fTrueEndX = fParticleMap[ fBestTrackID ].endX;
  fTrueEndY = fParticleMap[ fBestTrackID ].endY;
  fTrueEndZ = fParticleMap[ fBestTrackID ].endZ;

  fTrueStartDirectionX = fParticleMap[ fBestTrackID ].startDirection.X();
  fTrueStartDirectionY = fParticleMap[ fBestTrackID ].startDirection.Y();
  fTrueStartDirectionZ = fParticleMap[ fBestTrackID ].startDirection.Z();
  fTrueEndDirectionX = fParticleMap[ fBestTrackID ].endDirection.X();
  fTrueEndDirectionY = fParticleMap[ fBestTrackID ].endDirection.Y();
  fTrueEndDirectionZ = fParticleMap[ fBestTrackID ].endDirection.Z();

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

  fNumberOfHitsView0 = fHitEfficiency->getNumOfHits( fBestTrackID, 0 );
  fNumberOfHitsView1 = fHitEfficiency->getNumOfHits( fBestTrackID, 1 );
  fNumberOfHits = fNumberOfHitsView0+fNumberOfHitsView1;
  fAvgSignalToNoiseView0 = fHitEfficiency->getAvgSignalToNoise( fBestTrackID, 0 );
  fAvgSignalToNoiseView1 = fHitEfficiency->getAvgSignalToNoise( fBestTrackID, 1 );
  fAvgIntegralToNoiseView0 = fHitEfficiency->getAvgIntegralToNoise( fBestTrackID, 0 );
  fAvgIntegralToNoiseView1 = fHitEfficiency->getAvgIntegralToNoise( fBestTrackID, 1 );
  fAvgIntegralView0 = fHitEfficiency->getAvgIntegral( fBestTrackID, 0 );
  fAvgIntegralView1 = fHitEfficiency->getAvgIntegral( fBestTrackID, 1 );
  fAvgAmplitudeView0 = fHitEfficiency->getAvgAmplitude( fBestTrackID, 0 );
  fAvgAmplitudeView1 = fHitEfficiency->getAvgAmplitude( fBestTrackID, 1 );
  fAvgWidthView0 = fHitEfficiency->getAvgWidth( fBestTrackID, 0 );
  fAvgWidthView1 = fHitEfficiency->getAvgWidth( fBestTrackID, 1 );

  fRecoTTree->Fill();

  fBestTrackIDs.push_back( fBestTrackID );

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

  if( fTrack.completeness > 0.5 & fTrack.purity > 0.5 ){

    //fill the reco quantities
    this->fillMap1D( abs(fPdg), fThetaRecoMap, fParticleMap[fBestTrackID].startTheta);
    this->fillMap1D( abs(fPdg), fPhiRecoMap, fParticleMap[fBestTrackID].startPhi );
    this->fillMap2D( abs(fPdg), fPhiThetaRecoMap, fParticleMap[fBestTrackID].startTheta, fParticleMap[fBestTrackID].startPhi);
  }

}

void Efficiency::checkUnmatch()
{
  //build up a TTree with all the particles in the PartcileList that didn't get a match with a reco track

  for( auto const & particleIt : fParticleMap )
  {
    int id = particleIt.first;

    if ( std::find( fBestTrackIDs.begin(), fBestTrackIDs.end(), id ) == fBestTrackIDs.end() )
    {
      fEvent = fParticleMap[ id ].eventNumber;
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
      fTrueLength = fParticleMap[ id ].length;

      fUnmatchTTree->Fill();
    }
  }
}

void Efficiency::makeEfficiencyPlot()
{
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

void Efficiency::write()
{

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
  fUnmatchTTree->Write();

  if( fMakeHitsTree ) { fHitTree->Write(); }

}

void Efficiency::clean()
{

  //quantities to be reset after the event
  fParticleMap.clear();
  fHits.clear();
  fBestTrackIDs.clear();
  fHitsMap.clear();
  delete fHitEfficiency;
}
