////////////////////////////////////////////////////////////////////////////////
// Reconstruction efficiency class definition
//
//mailto:andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////

#ifndef  __EFFICIENCY_H
#define __EFFICIENCY_H

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


class Efficiency
{
  public:
    Efficiency();
    Efficiency( int fileNumber );
    ~Efficiency();

    //setters
    void setMapEntry(int id, MCTrack mctrack );
    void setRecoTrack( Track track ){ fTrack = track; }
    void setRecoHits( vector<Hit> hits ){ fHits = hits; }
    void setNumberOfTracksEvent( int nTracks ){ fNtracksEvent = nTracks;} //not elegant, but need that afterwards in the analysis

    //others
    void makeEfficiencyPlot();
    void fill();
    void write();

    //cleaner
    void clean();

  private:
    void initClass();
    void matchTruth();

    void fillMap1D(int pdg, map<int, TH1D*> map, double fillIn );
    void fillMap2D(int pdg, map<int, TH2D*> map, double fillX, double fillY );


    //Particle that I consider for the efficiency
    vector<int> pdgCode = { 13, 11, 211, 2212, 0 }; // NB this can be custom set

    //NB: this can be fetched from a database (has ROOT something already)
    vector<string> pdgNames = { "Muons", "Electrons","Pions", "Protons", "Other" };

    //binning

    //dir
    int nBinsDir = 100;

    //pos
    int nBinsPos = 100;

    //theta
    int nBinsTheta = 100;
    double thetaStart = 0;
    double thetaEnd = 180;

    //phi
    int nBinsPhi = 100;
    double phiStart = -180;
    double phiEnd = 180;

    //histogram maps
    map<int, TH1D*> fThetaG4Map;
    map<int, TH1D*> fPhiG4Map;
    map<int, TH1D*> fDirMap;
    map<int, TH2D*> fDirThetaMap;
    map<int, TH2D*> fDirPhiMap;
    map<int, TH1D*> fPosXMap;
    map<int, TH2D*> fPosXThetaMap;
    map<int, TH2D*> fPosXPhiMap;
    map<int, TH1D*> fPosYMap;
    map<int, TH2D*> fPosYThetaMap;
    map<int, TH2D*> fPosYPhiMap;
    map<int, TH1D*> fPosZMap;
    map<int, TH2D*> fPosZThetaMap;
    map<int, TH2D*> fPosZPhiMap;
    map<int, TH1D*> fThetaTrueMap;
    map<int, TH1D*> fPhiTrueMap;
    map<int, TH2D*> fPhiThetaTrueMap;
    map<int, TH1D*> fThetaRecoMap;
    map<int, TH1D*> fPhiRecoMap;
    map<int, TH2D*> fPhiThetaRecoMap;
    map<int, TEfficiency*> fPhiEfficiency;
    map<int, TEfficiency*> fThetaEfficiency;
    map<int, TEfficiency*> fPhiThetaEfficiency;

    //more quantities
    TTree *fMcTTree;
    TTree *fRecoTTree;

    map<int, MCTrack> fParticleMap; //particleID mctrack association
    map<int, double> fEnergyMap; //particleID energy association
    Track fTrack;
    vector<Hit> fHits;

    int fFileNumber;
    int fEvent;
    int fUniqueEventLabel;
    int fParticleId;
    int fTrackId;
    int fNtracksEvent;
    double fPurirty;
    double fCompleteness;
    double fPdg;
    double fBestTrackID;
    double fTruePhi;
    double fTrueTheta;
    double fTrueE;
    double fTrueStartX;
    double fTrueStartY;
    double fTrueStartZ;
    double fTrueEndX;
    double fTrueEndY;
    double fTrueEndZ;
    double fTrueStartDirectionX;
    double fTrueStartDirectionY;
    double fTrueStartDirectionZ;
    double fTrueEndDirectionX;
    double fTrueEndDirectionY;
    double fTrueEndDirectionZ;
    double fRecoPhi;
    double fRecoTheta;
    double fRecoE;

    double fDirection;
    double fDiffStartX;
    double fDiffStartY;
    double fDiffStartZ;

    double fStartX;
    double fStartY;
    double fStartZ;

    double fEndX;
    double fEndY;
    double fEndZ;

    double fStartDirectionX;
    double fStartDirectionY;
    double fStartDirectionZ;

    double fEndDirectionX;
    double fEndDirectionY;
    double fEndDirectionZ;

    double fRecoLength;
    double fTrueLength;
};

#endif // __EFFICIENCY_H
