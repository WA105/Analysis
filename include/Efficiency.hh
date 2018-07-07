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
    void write();

    //cleaner
    void clean();

  private:
    void matchTruth();

    //more quantities
    TTree *fMcTTree;
    TTree *fRecoTTree;

    map<int, MCTrack> fParticleMap; //particleID mctrack association
    map<int, double> fEnergyMap; //particleID energy association
    Track fTrack;
    vector<Hit> fHits;
    double fPurirty;
    double fCompleteness;
    double fPdg;
    double fBestTrackID;
    double fTruePhi;
    double fTrueTheta;
    double fTrueE;
    double fRecoPhi;
    double fRecoTheta;
    double fRecoE;
};

#endif // __EFFICIENCY_H
