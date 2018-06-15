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

class Cut
{
  public:
    Cut();
    Cut( vector<int> volumeCut, double length );
    ~Cut();

    //setters
    void setActiveVolumeCut( vector<int> volumeCut );
    void setTrackLengthCut(double length){fLengthCut = length; };

    //cuts
    bool isCrossingTPC( Track track ); //crossing TPC from anode to cathode
    bool passBoxCut( Track track ); //box cut for charge

  private:

    void initActiveVolumeBounds();

    vector<int> fVolumeCut = {0,0,0,0,0,0}; //in cm
    double fLengthCut=100; //in cm

    int minx;
    int maxx;
    int miny;
    int maxy;
    int minz;
    int maxz;

};

#endif //__CUTS_H
