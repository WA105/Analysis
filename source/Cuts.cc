////////////////////////////////////////////////////////////////////////////////
// Cuts class methods
//
// This class holds the cuts the user may want to apply on data
//
//mailto:andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "Geometry.hh"
#include "DataStructure.hh"
#include "Cuts.hh"

Cut::Cut(){
  //standard constructor

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

bool Cut::isPassingCut(){
  return fIsPassingCut;
}

bool Cut::isCrossingTPC(){
  //Check if the selected track cross the tpc

  fIsPassingCut = false;

  //first select based on the track length:
  if( fTrack.length > fLengthCut){
    if( max(fTrack.endPointX, fTrack.startPointX) > maxx
        && min(fTrack.endPointX, fTrack.startPointX) < minx ){
          fIsPassingCut = true;
    }
  }

  return fIsPassingCut;
}

/*
void Cut::cutsFromFile( string filename ){
  //import an external list of cuts from a file, make a map with a list of cuts

}

bool Cut::checkCutFile( Track *track ){
  //check if the track is among the list of cuts imported

}
*/
