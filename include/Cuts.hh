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
    void setTrackLengthCut(double length){fLengthCut = length; };
    void setTrack( Track track ){ fTrack = track; }

    //getters
    bool isPassingCut( Track track ); //crossing TPC from anode to cathode
    bool isCrossingTPC( Track track );
    bool cutEndPoint( Track track );

    //cut type
    void setActiveVolumeCut( vector<int> volumeCut );
    bool cutsFromFile( string filename );
    bool appyNoCuts(){ fSelNoCuts = true; }
    bool cutsFromHighway( string filename, double allCBR, double sumCBR );

  private:

    //here one can implement its own cut functions
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

    vector<Track> fCutsList;
    vector<Track> fHighwayCutsList;
    bool fCutsFromFile=false;
    bool fCutsFromHighway=false;
    bool fSelNoCuts=false;

};

#endif //__CUTS_H
