////////////////////////////////////////////////////////////////////////////////
// Data structures class definition
//
//mailto:andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////

#ifndef  __DATASTRUCTURE_H
#define __DATASTRUCTURE_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "TTree.h"
#include "Geometry.hh"
#include "Run.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//Hit data structure
class Hit
{

  public:
    Hit();
    ~Hit();

    Run run;
    int subRun;
    int event;
    int trackID;
    double X;
    double Y;
    double Z;
    double dxLocalTrackDirection;
    double dx3DPosition;
    double TPC;
    double view;
    double channel;
    double peakTime;
    double chargeSummedADC;
    double chargeIntegral;
    double startTime;
    double endTime;
    double width;
    double goodnessOfFit;
    double multiplicity;

    //find lem
    //project view
    //deprecated channels
};

////////////////////////////////////////////////////////////////////////////////
//Track data Strucutre
class Track
{

  public:
    Track();
    ~Track();

    Run run;
    int subRun;
    int event;
    int trackID;
    int numberOfHits;

    double length;
    double startPointX;
    double startPointY;
    double startPointZ;
    double startPointDistanceToBoundary;
    double endPointX;
    double endPointY;
    double endPointZ;
    double endPointDistanceToBoundary;
    double startDirectionTheta;
    double startDirectionPhi;
    double startDirectionX;
    double startDirectionY;
    double startDirectionZ;
    double endDirectionTheta;
    double endDirectionPhi;
    double endDirectionX;
    double endDirectionY;
    double endDirectionZ;

    vector<Hit> hitsTrk; //Array holding the hits associated with this track
};

////////////////////////////////////////////////////////////////////////////////
//Geant4 Tracks data Strucutre
class MCTrack : public Track
{

  public:
    MCTrack();
    ~MCTrack();


  private:

};

////////////////////////////////////////////////////////////////////////////////

//similary here a new class can be created for inputs from qScan

////////////////////////////////////////////////////////////////////////////////

//load the LArSoft Parser files into the data structures
class LArParser
{
  public:
    LArParser();
    LArParser( TTree *tree, Run *run );
    ~LArParser();

    //setters
    void setTTree( TTree *tree ){ fTree = tree; }
    void setRun(Run *run){fRun = run;}

    //getters
    void getRecoHitsEvent( vector<Hit> & hits, int event  );
    void getRecoTracksEvent( vector<Track> & tracks, int event  );
    void getRecoHits( vector<Hit> & hits );
    void getRecoTracks( vector<Track> & tracks );
    //void readG4Tree( &vector<MCTracks> tracks );

    bool isTreeGood( );
  private:

    void setRecoBranches();
    void fillRecoHits( vector<Hit> & hits );
    void fillRecoTrack( vector<Track> & tracks );

    TTree *fTree;
    Run *fRun;

    static const int NMaxHitsPerEvent=10000;;
    static const int NMaxClustersPerEvent=10000;
    static const int NMaxTracksPerEvent=1000;
    static const int NMaxTracksPerEventTimesNViews=NMaxTracksPerEvent*NUM_OF_VIEWS;
    static const int NEventsPerRun=335;

    //Define variables to store the data of the ROOT file
    //Metadata
    int tRun;
    int tSubrun;
    int tEventNumberInRun;
    //  int  tEventTimeSeconds;
    //  int  tEventTimeNanoseconds;
    //  char tIsData;

    //Hit variables

    int tNumberOfHits;
    short tHit_TPC[NMaxHitsPerEvent];
    short tHit_View[NMaxHitsPerEvent];
    short tHit_Channel[NMaxHitsPerEvent];
    float tHit_PeakTime[NMaxHitsPerEvent];
    float tHit_ChargeSummedADC[NMaxHitsPerEvent];
    float tHit_ChargeIntegral[NMaxHitsPerEvent];
    float tHit_StartTime[NMaxHitsPerEvent];
    float tHit_EndTime[NMaxHitsPerEvent];
    float tHit_Width[NMaxHitsPerEvent];
    float tHit_GoodnessOfFit[NMaxHitsPerEvent];
    short tHit_Multiplicity[NMaxHitsPerEvent];
    short tHit_TrackID[NMaxHitsPerEvent];
    short tHit_ClusterID[NMaxHitsPerEvent];


    //Cluster variables
    //  short tNumberOfClusters;
    //  short tClusterID[NMaxClustersPerEvent];
    //  short tCluster_NumberOfHits[NMaxClustersPerEvent];
    //  short tCluster_View[NMaxClustersPerEvent];
    //  float tCluster_ChargeIntegral[NMaxHitsPerEvent];
    //  short tCluster_StartChannel[NMaxHitsPerEvent];
    //  short tCluster_StartTick[NMaxHitsPerEvent];
    //  short tCluster_EndChannel[NMaxHitsPerEvent];
    //  short tCluster_EndTick[NMaxHitsPerEvent];
    //  float tCluster_StartCharge[NMaxHitsPerEvent];
    //  float tCluster_StartAngle[NMaxHitsPerEvent];
    //  float tCluster_EndCharge[NMaxHitsPerEvent];
    //  float tCluster_EndAngle[NMaxHitsPerEvent];

    //Track variables
    short tNumberOfTracks;
    short tTrackID[NMaxTracksPerEventTimesNViews];
    short tTrack_NumberOfHits[NMaxTracksPerEventTimesNViews];
    float tTrack_Length[NMaxTracksPerEventTimesNViews];
    float tTrack_StartPoint_X[NMaxTracksPerEvent];
    float tTrack_StartPoint_Y[NMaxTracksPerEvent];
    float tTrack_StartPoint_Z[NMaxTracksPerEvent];
    float tTrack_StartPoint_DistanceToBoundary[NMaxTracksPerEvent];
    float tTrack_EndPoint_X[NMaxTracksPerEvent];
    float tTrack_EndPoint_Y[NMaxTracksPerEvent];
    float tTrack_EndPoint_Z[NMaxTracksPerEvent];
    float tTrack_EndPoint_DistanceToBoundary[NMaxTracksPerEvent];
    float tTrack_StartDirection_Theta[NMaxTracksPerEvent];
    float tTrack_StartDirection_Phi[NMaxTracksPerEvent];
    float tTrack_StartDirection_X[NMaxTracksPerEvent];
    float tTrack_StartDirection_Y[NMaxTracksPerEvent];
    float tTrack_StartDirection_Z[NMaxTracksPerEvent];
    float tTrack_EndDirection_Theta[NMaxTracksPerEvent];
    float tTrack_EndDirection_Phi[NMaxTracksPerEvent];
    float tTrack_EndDirection_X[NMaxTracksPerEvent];
    float tTrack_EndDirection_Y[NMaxTracksPerEvent];
    float tTrack_EndDirection_Z[NMaxTracksPerEvent];
    float tTrack_PitchInViews[NMaxTracksPerEvent];
    short tTrack_NumberOfHitsPerView[NMaxTracksPerEvent][2];

    //Track hit variables
    float tTrack_Hit_X[NMaxHitsPerEvent];
    float tTrack_Hit_Y[NMaxHitsPerEvent];
    float tTrack_Hit_Z[NMaxHitsPerEvent];
    float tTrack_dx_LocalTrackDirection[NMaxHitsPerEvent];
    float tTrack_dx_3DPosition[NMaxHitsPerEvent];
    short tTrack_Hit_TPC[NMaxHitsPerEvent];
    short tTrack_Hit_View[NMaxHitsPerEvent];
    short tTrack_Hit_Channel[NMaxHitsPerEvent];
    float tTrack_Hit_PeakTime[NMaxHitsPerEvent];
    float tTrack_Hit_ChargeSummedADC[NMaxHitsPerEvent];
    float tTrack_Hit_ChargeIntegral[NMaxHitsPerEvent];
    float tTrack_Hit_StartTime[NMaxHitsPerEvent];
    float tTrack_Hit_EndTime[NMaxHitsPerEvent];
    float tTrack_Hit_Width[NMaxHitsPerEvent];
    float tTrack_Hit_GoodnessOfFit[NMaxHitsPerEvent];
    short tTrack_Hit_Multiplicity[NMaxHitsPerEvent];

};

#endif // __DATASTRUCTURE_H
