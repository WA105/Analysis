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
    ~Efficiency();

    //setters
    void setMapEntry(int id, MCTrack mctrack ){ fParticleMap[id] = mctrack; }
    void setRecoTrack( Track track ){ fTrack = track; }
    void setRecoHits( vector<Hit> hits ){ fHits = hits; }

    //getters

    //others
    void fill();
    void makeEfficiencyPlot();
    void write(TFile *ofile);

    //cleaner
    void clean();

  private:
    void matchTruth();
    void fillMap1D(int pdg, map<int, TH1D*> map, double fillIn );
    void fillMap2D(int pdg, map<int, TH2D*> map, double fillX, double fillY );

    //Particle that I consider for the efficiency
    vector<int> pdgCode = { 13, -13, 11, -11, 211, -211, 2212, 0 }; // NB this can be custom set

    //NB: this can be fetched from a database (has ROOT something already)
    vector<string> pdgNames = { "Muons", "Antimuons", "Electrons", "Positron", "PiPlus", "PiMinus", "Proton", "Other" };

    //binning
    //theta
    int nBinsTheta = 100;
    double thetaStart = 0;
    double thetaEnd = 180;

    //phi
    int nBinsPhi = 100;
    double phiStart = -180;
    double phiEnd = 180;

    //histogram maps
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
    map<int, MCTrack> fParticleMap; //particleID mctrack association
    map<int, double> fEnergyMap; //particleID energy association
    Track fTrack;
    vector<Hit> fHits;
    double fPurirty;
    double fCompleteness;
    double fPdg;
    double fBestTrackID;

};

#endif // __EFFICIENCY_H
