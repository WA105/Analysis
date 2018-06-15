////////////////////////////////////////////////////////////////////////////////
// Cuts class methods
//
//This class holds the cuts the user may want to apply on data
//
//mailto:andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "Geometry.hh"
#include "Cuts.hh"

Cut::Cut(){
  //standard constructor
  this->initActiveVolumeBounds();
}

Cut::Cut( vector<int> volumeCut, double length ){
  //alternative constructor
  this->setActiveVolumeCut( volumeCut );
  this->setTrackLengthCut( length );
}

Cut::~Cut(){};

void Cut::initActiveVolumeBounds(){
  //initialize the active volume bounds
  minx = tpc_boundaries[0] + fVolumeCut[0];
  maxx = tpc_boundaries[1] - fVolumeCut[1];
  miny = tpc_boundaries[2] + fVolumeCut[2];
  maxy = tpc_boundaries[3] - fVolumeCut[3];
  minz = tpc_boundaries[4] + fVolumeCut[4];
  maxz = tpc_boundaries[5] - fVolumeCut[5];

  return;
}

void Cut::setActiveVolumeCut( vector<int> volumeCut ){
  //set a cut over the active volume

  if( volumeCut.size() != 6 ){
    cout << "Cut::setActiveVolumeCut() ERROR: active volume array must contain 6 elements " <<endl;
  }

  for(int i=0; i<6; i++ ){
    fVolumeCut.at(i) = volumeCut.at(i);
  }

  //make it effective
  this->initActiveVolumeBounds();

  return;
}

bool Cut::isCrossingTPC( Track track ){
  //Check if the selected track cross the tpc

  bool isCrossing = false;

  //first select based on the track length:
  if( track.length > fLengthCut){
    if( max(track.endPointX, track.startPointX) > maxx
        && min(track.endPointX, track.startPointX) < minx ){
          isCrossing = true;
    }
  }

  return isCrossing;
}

/*
bool Cut::passBoxCut( Track track ){

}
*/
