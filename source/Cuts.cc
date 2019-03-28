////////////////////////////////////////////////////////////////////////////////
// Cuts class methods
//
// This class holds the cuts the user may want to apply on data
//
// mailto:andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>

#include <bits/stdc++.h>
#include <boost/algorithm/string.hpp>

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

bool Cut::isPassingCut( Track testTrack ){

  fIsPassingCut=false;

  //----------------------------------------------------------------------------
  if( fCutsFromFile )
  {
    //go inside the list and check if this track matches the correct criteria
    auto lambda = [ testTrack ]( Track cutTrack )
    {
       bool runCheck = (testTrack.run == cutTrack.run);
       bool subrunCheck = (testTrack.subRun == cutTrack.subRun);
       bool eventCheck = (testTrack.event == cutTrack.event);
       bool idCheck = (testTrack.trackID == cutTrack.trackID);

       return runCheck && subrunCheck && eventCheck && idCheck;
    };
    auto pos = find_if( fCutsList.begin(),fCutsList.end(), lambda );

    if( pos != fCutsList.end() ){ fIsPassingCut = true; }
    else{ fIsPassingCut = false; }
  }

  //----------------------------------------------------------------------------
  if ( fCutsFromHighway )
  {
    //go inside the list and check if this track matches the correct criteria
    auto lambda = [ testTrack ]( Track cutTrack )
    {
       bool runCheck = (testTrack.run == cutTrack.run);
       bool subrunCheck = (testTrack.subRun == cutTrack.subRun);
       bool eventCheck = (testTrack.event == cutTrack.event);
       bool idCheck = (testTrack.trackID == cutTrack.trackID);

       return runCheck && subrunCheck && eventCheck && idCheck;
    };
    auto pos = find_if( fHighwayCutsList.begin(),fHighwayCutsList.end(), lambda );

    if( pos != fHighwayCutsList.end() ){ fIsPassingCut = true; }
    else{ fIsPassingCut = false; }
  }

  //----------------------------------------------------------------------------
  if( fSelNoCuts ){ fIsPassingCut=true; } //I am not applying any cut

  return fIsPassingCut;
}

bool Cut::isCrossingTPC( Track track ){
  //Check if the selected track cross the tpc

  bool isCrossing = false;

  //first select based on the track length:
  if( max(track.endPointX, track.startPointX) > maxx
        && min(track.endPointX, track.startPointX) < minx ){
          isCrossing = true;
  }

  return isCrossing;
}

bool Cut::cutEndPoint( Track track ){
  //Check if the selected track cross the tpc

  bool isCrossing = false;

  //first select based on the track length:
  if( track.endPointX <= minx ){
          isCrossing = true;
  }

  return isCrossing;
}

bool Cut::cutsFromHighway( string filename, double allCBR, double sumCBR ){
  //read tracks selected by the highway algoritm, by applying cuts on the charge
  //box ratios. 'allCBR' will apply the cut on each box, 'sumCBR' will apply the
  //cut on the sum of all boxes

  //clean the vector
  fHighwayCutsList.clear();

  //cut-type flag
  fCutsFromHighway = false;

  //read file
  ifstream infile;

  infile.open(filename.c_str());
  if(infile.fail()) // checks to see if file opended
  {
      cout << "Cut::cutsFromFile: File doesn't exist!" << endl;
      return fCutsFromHighway;
  }

  string firstline, line;

  //skip the very first line
  getline(infile, firstline);

  //now read all the other lines
  while(getline(infile, line))
  {
    Track dummyTrack;
    std::stringstream ss(line);

    int run, subrun, event, id;
    double csb0[10];
    double clb0[10];
    double csb1[10];
    double clb1[10];

    ss >> run >> subrun >> event >> id;
    dummyTrack.run = run;
    dummyTrack.subRun = subrun;
    dummyTrack.event = subrun*335+event;
    dummyTrack.trackID = id;

    //read cb small v0
    for( int i=0; i < 10; i++ ) { ss >> csb0[i]; }

    //read cb large v0
    for( int i=0; i < 10; i++ ){ ss >> clb0[i]; }

    //read cb small v1
    for( int i=0; i < 10; i++ ){ ss >> csb1[i]; }

    //read cb large v1
    for( int i=0; i < 10; i++ ){ ss >> clb1[i]; }

    //calculate the cbr for all the boxes and check if match with the cuts
    bool fisGood = true;
    for( int i=0; i<10; i++ )
    {
      if ( clb0[i] == -9999 || clb1[i] ==  -9999 ||
           csb0[i] == -9999 || csb1[i] ==  -9999 ) { continue; }

      double cbr0 = ( clb0[i] - csb0[i] )/csb0[i];
      double cbr1 = ( clb1[i] - csb1[i] )/csb1[i];

      if( cbr0 > allCBR || cbr1 > allCBR || cbr0 < -0.1 || cbr1 < -0.1 )
      {
        fisGood=false;
        break;
      }
    }

    //TODO: calculate the summed charge box ratio

    if (fisGood){ fHighwayCutsList.push_back(dummyTrack); }
  }

  infile.close();
  fCutsFromHighway=true;

  cout << "File read " << endl;

  return fCutsFromHighway;
}

bool Cut::cutsFromFile( string filename ){
  //import an external list of cuts from a file, make a map with a list of cuts

  //clean the vector
  fCutsList.clear();

  //cut-type flag
  fCutsFromFile = false;

  //read file
  ifstream infile;

  infile.open(filename.c_str());
  if(infile.fail()) // checks to see if file opended
  {
      cout << "Cut::cutsFromFile: File doesn't exist!" << endl;
      return fCutsFromFile;
  }
  string line, firstline;

  //skip the very first line
  getline(infile, firstline);

  while(getline(infile, line))
  {

    Track dummyTrack; // define a dummyTrack object with minimal inport quantites
    std::stringstream ss(line);
    int evtid, run, subrun, event, id;

    if ( ss >> evtid >> run >> subrun >> event >> id )
    {
        dummyTrack.run = run;
        dummyTrack.subRun = subrun;
        dummyTrack.event = subrun*335+event;
        dummyTrack.trackID = id;

        fCutsList.push_back(dummyTrack);
      }
  }

  infile.close();
  fCutsFromFile=true;

  cout << "File read " << endl;

  return fCutsFromFile;
}
