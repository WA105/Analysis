////////////////////////////////////////////////////////////////////////////////
// Cuts class definition
//
//This class holds the cuts the user may want to apply on data
//
//mailto:andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////

#ifndef  __CUTS_H
#define __CUTS_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "DataStructure.hh"
#include "Geometry.hh"

using namespace std;

//Cut implement Track, therefore all the public classes and variables that are
//defined into Track are available also in Cut

class Cut
{
  public:
    Cut();
    Cut( vector<int> volumeCut, double length );
    ~Cut();

    //setters
    void setActiveVolumeCut( vector<int> volumeCut );
    void setTrackLengthCut(double length){fLengthCut = length; };
    void setTrack( Track track ){ fTrack = track; }

    //getters
    bool isPassingCut(); //crossing TPC from anode to cathode

  private:

    //here one can implement its own cut functions
    bool isCrossingTPC();
    void initActiveVolumeBounds();

    vector<int> fVolumeCut = {0,0,0,0,0,0}; //in cm
    double fLengthCut=100; //in cm

    bool fIsPassingCut; //flag retuned by isPassingCut() into the track object

    Track fTrack;

    int minx;
    int maxx;
    int miny;
    int maxy;
    int minz;
    int maxz;

};

#endif //__CUTS_H
