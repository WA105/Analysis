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
#include "TVector3.h"

#include "Geometry.hh"
#include "Run.hh"
//#include "Cuts.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//Channel data structure: holds raw adc vectors deposited on each wire

class Channel
{
  public:
    Channel();
    ~Channel();

    int run;
    int subRun;
    int event;
    int timeSeconds;
    int timeNanoSeconds;
    int channel;
    int nTicks;
    int view;
    std::vector<float> signal;

    bool isDead();
    bool isBad();
    void subtractPedestal( bool subtractPedestal );
    double sumAdcInROI( int startTime, int endTime );

};

////////////////////////////////////////////////////////////////////////////////
//Hit data structure
class Hit
{

  public:
    Hit();
    ~Hit();

    int run;
    int subRun;
    int event;
    int trackID;
    double X;
    double Y;
    double Z;
    double dxLocalTrackDirection;
    double dx3DPosition;
    double TPC;
    int view;
    double channel;
    double peakTime;
    double amplitude;
    double chargeSummedADC;
    double chargeIntegral;
    double startTime;
    double endTime;
    double width;
    double goodnessOfFit;
    double multiplicity;
    double theta;
    double phi;

    //if available
    double channelNoise;

    //backtracker quantites
    int particleID;
    double trueEnergy;
    double trueEnergyFraction;
    double chanElectrons;
    double chanMaxElectrons;
    int chanBestID;

    int lem;
    //project view
    //deprecated channels
    void calibrateCharge( string type );
    double getdQdx(string type);
    void findLem();
    bool isGoodLem( vector<int> lems );
    double electronsFromSumADC( string type );
    double electronsFromIntegral( string type );
    double getPeakToNoise(){ return  channelNoise > 0 ? amplitude/channelNoise : 0; };
    double getIntegralToNoise( string type ){ return  channelNoise > 0 ? electronsFromIntegral( type )/channelNoise*900 : 0; };
    double getSummedADCToNoise(){ return  channelNoise > 0 ? chargeSummedADC/channelNoise : 0; };

};

////////////////////////////////////////////////////////////////////////////////
//Track data Strucutre
class Track
{

  public:
    Track();
    ~Track();

    int run;
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
    TVector3 startDirection;
    TVector3 endDirection;
    double endDirectionTheta;
    double endDirectionPhi;
    double endDirectionX;
    double endDirectionY;
    double endDirectionZ;

    int bestParticleID;
    double completeness;
    double purity;

    vector<Hit> hitsTrk; //Array holding the hits associated with this track

};

////////////////////////////////////////////////////////////////////////////////
//Geant4 Tracks data Strucutre
class MCTrack
{

  public:
    MCTrack();
    ~MCTrack();

    int run;
    int subrun;
    int eventNumber;
    int mcTracksPerEvent;
    int pdgCode;
    int particleID;
    int isInTPCAV;
    double startE;
    double endE;
    double mom;
    double momX;
    double momY;
    double momZ;
    double startTime;
    double startX;
    double startY;
    double startZ;
    double startTheta;
    double startPhi;
    double lengthAV;
    double endX;
    double endY;
    double endZ;
    double length;
    TVector3 startDirection;
    TVector3 endDirection;

  };

////////////////////////////////////////////////////////////////////////////////

// similary here a new class can be created to convert inputs from qScan into
// the data structure holds by this framework

////////////////////////////////////////////////////////////////////////////////

//load the LArSoft Parser files into the data structures
class LArParser
{
  public:
    LArParser();
    ~LArParser();

    LArParser(bool activeBranch);

    //setters
    //void setTTree( TTree *tree ){ fTree = tree; }
    //void setRun(Run *run){fRun = run;}
    void setChannelNoise( string model );

    //getters
    void getRawChannelsEvent(TTree *tree, vector<Channel> & channels, int event  );
    void getRecoChannelsEvent(TTree *tree, vector<Channel> & channels, int event  );
    void getMCTracksEvent(TTree *tree, vector<MCTrack>  & mctracks, int event );
    void getRecoHitsEvent(TTree *tree, vector<Hit> & hits, int event );
    void getRecoTracksEvent(TTree *tree, vector<Track> & tracks, int event );
    void getMCTracks(TTree *tree, vector<MCTrack>  & mctracks );
    void getRecoHits(TTree *tree, vector<Hit> & hits );
    void getRecoTracks(TTree *tree, vector<Track> & tracks );
    //void readG4Tree( &vector<MCTracks> tracks );

    bool isTreeGood(TTree *tree);

  private:
    void setRawBranches(TTree *tree);
    void setRecoChannelBranches(TTree *tree);
    void setRecoHitBranches(TTree *tree);
    void setRecoClusterBranches(TTree *tree);
    void setRecoTrackBranches(TTree *tree);
    void setMCBranches(TTree *tree);
    void fillRawChannels( vector<Channel> & channels );
    void fillRecoChannels( vector<Channel> & channels );
    void fillMCTrack( vector<MCTrack> & tracks );
    void fillRecoHits( vector<Hit> & hits );
    void fillRecoTrack( vector<Track> & tracks );
    //void clean();

    static const int NMaxGeantTrackPerEvent=10000;
    static const int NMaxHitsPerEvent=10000;
    static const int NMaxClustersPerEvent=10000;
    static const int NMaxTracksPerEvent=1000;
    static const int NMaxTracksPerEventTimesNViews=NMaxTracksPerEvent*NUM_OF_VIEWS;
    static const int NEventsPerRun=335;
    static const int nMaxNumChannels= maxNumChannels;
    static const int nMaxNumTdc = 2133760;

    //Common ///////////////////////////////////////////////////////////////////
    int tRun;
    int tSubrun=0;
    int tEventNumberInRun=0;
    int tEventTimeSeconds=0;
    int tEventTimeNanoseconds=0;

    //Raw //////////////////////////////////////////////////////////////////////
    int tRawWaveform_NumberOfChannels=0;
    int tRawWaveform_NumberOfTicks=0;
    int tRawWaveform_Channel[nMaxNumChannels]={0};
    int tRawWaveform_NumberOfTicksInAllChannels[nMaxNumChannels] = {0};
    short tRawWaveform_ADC[nMaxNumTdc] = {0};

    //G4 //////////////////////////////////////////////////////////////////////
    int tNGeantTrackPerEvent=0;
    int tPdg[NMaxGeantTrackPerEvent]={0};
    int tParticleId[NMaxGeantTrackPerEvent]={0};
    int tIsInTPCAV[NMaxGeantTrackPerEvent]={0};
    float tStartE[NMaxGeantTrackPerEvent]={0};
    float tEndE[NMaxGeantTrackPerEvent]={0};
    float tMom[NMaxGeantTrackPerEvent]={0};
    float tMomX[NMaxGeantTrackPerEvent]={0};
    float tMomY[NMaxGeantTrackPerEvent]={0};
    float tMomZ[NMaxGeantTrackPerEvent]={0};
    float tStartTime[NMaxGeantTrackPerEvent]={0};
    float tStartX[NMaxGeantTrackPerEvent]={0};
    float tStartY[NMaxGeantTrackPerEvent]={0};
    float tStartZ[NMaxGeantTrackPerEvent]={0};
    float tStartPhi[NMaxGeantTrackPerEvent]={0};
    float tStartTheta[NMaxGeantTrackPerEvent]={0};
    float tLengthAV[NMaxGeantTrackPerEvent]={0};
    float tEndX[NMaxGeantTrackPerEvent]={0};
    float tEndY[NMaxGeantTrackPerEvent]={0};
    float tEndZ[NMaxGeantTrackPerEvent]={0};

    //Reco /////////////////////////////////////////////////////////////////////

    //reco wires
    int tRecoWaveform_NumberOfChannels=0;
    int tRecoWaveform_NumberOfTicks[nMaxNumChannels]={0};
    int tRecoWaveform_Channel[nMaxNumChannels]={0};
    int tRecoWaveform_NumberOfTicksInAllChannels=0;
    float tRecoWaveform_ADC[nMaxNumTdc]={0};

    //Hit variables
    int tNumberOfHits=0;
    short tHit_TPC[NMaxHitsPerEvent]={0};
    short tHit_View[NMaxHitsPerEvent]={0};
    short tHit_Channel[NMaxHitsPerEvent]={0};
    float tHit_PeakTime[NMaxHitsPerEvent]={0};
    float tHit_Amplitude[NMaxHitsPerEvent]={0};
    float tHit_ChargeSummedADC[NMaxHitsPerEvent]={0};
    float tHit_ChargeIntegral[NMaxHitsPerEvent]={0};
    float tHit_StartTime[NMaxHitsPerEvent]={0};
    float tHit_EndTime[NMaxHitsPerEvent]={0};
    float tHit_Width[NMaxHitsPerEvent]={0};
    float tHit_GoodnessOfFit[NMaxHitsPerEvent]={0};
    short tHit_Multiplicity[NMaxHitsPerEvent]={0};
    short tHit_TrackID[NMaxHitsPerEvent]={0};
    short tHit_ClusterID[NMaxHitsPerEvent]={0};
    int tHit_particleID[NMaxHitsPerEvent]={0};
    float tHit_TrueEnergy[NMaxHitsPerEvent]={0};
    float tHit_TrueEnergyFraction[NMaxHitsPerEvent]={0};
    float tHit_chanElectrons[NMaxHitsPerEvent] = {0};
    float tHit_chanMaxElectrons[NMaxHitsPerEvent] = {0};
    int tHit_bestChanIDE[NMaxHitsPerEvent] = {0};

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
    short tNumberOfTracks=0;
    short tTrackID[NMaxTracksPerEventTimesNViews]={0};
    short tTrack_NumberOfHits[NMaxTracksPerEventTimesNViews]={0};
    int tTrack_BestParticleID[NMaxTracksPerEventTimesNViews]={0};
    double tTrack_Completeness[NMaxTracksPerEventTimesNViews]={0};
    double tTrack_Purity[NMaxTracksPerEventTimesNViews]={0};
    float tTrack_Length[NMaxTracksPerEventTimesNViews]={0};
    float tTrack_StartPoint_X[NMaxTracksPerEvent]={0};
    float tTrack_StartPoint_Y[NMaxTracksPerEvent]={0};
    float tTrack_StartPoint_Z[NMaxTracksPerEvent]={0};
    float tTrack_StartPoint_DistanceToBoundary[NMaxTracksPerEvent]={0};
    float tTrack_EndPoint_X[NMaxTracksPerEvent]={0};
    float tTrack_EndPoint_Y[NMaxTracksPerEvent]={0};
    float tTrack_EndPoint_Z[NMaxTracksPerEvent]={0};
    float tTrack_EndPoint_DistanceToBoundary[NMaxTracksPerEvent]={0};
    float tTrack_StartDirection_Theta[NMaxTracksPerEvent]={0};
    float tTrack_StartDirection_Phi[NMaxTracksPerEvent]={0};
    float tTrack_StartDirection_X[NMaxTracksPerEvent]={0};
    float tTrack_StartDirection_Y[NMaxTracksPerEvent]={0};
    float tTrack_StartDirection_Z[NMaxTracksPerEvent]={0};
    float tTrack_EndDirection_Theta[NMaxTracksPerEvent]={0};
    float tTrack_EndDirection_Phi[NMaxTracksPerEvent]={0};
    float tTrack_EndDirection_X[NMaxTracksPerEvent]={0};
    float tTrack_EndDirection_Y[NMaxTracksPerEvent]={0};
    float tTrack_EndDirection_Z[NMaxTracksPerEvent]={0};
    float tTrack_PitchInViews[NMaxTracksPerEvent]={0};
    short tTrack_NumberOfHitsPerView[NMaxTracksPerEvent][2]={};

    //Track hit variables
    float tTrack_Hit_X[NMaxHitsPerEvent]={0};
    float tTrack_Hit_Y[NMaxHitsPerEvent]={0};
    float tTrack_Hit_Z[NMaxHitsPerEvent]={0};
    float tTrack_dx_LocalTrackDirection[NMaxHitsPerEvent]={0};
    float tTrack_Hit_LocalTrackDirection_Theta[NMaxHitsPerEvent]={0};
    float tTrack_Hit_LocalTrackDirection_Phi[NMaxHitsPerEvent]={0};
    float tTrack_dx_3DPosition[NMaxHitsPerEvent]={0};
    short tTrack_Hit_TPC[NMaxHitsPerEvent]={0};
    short tTrack_Hit_View[NMaxHitsPerEvent]={0};
    short tTrack_Hit_Channel[NMaxHitsPerEvent]={0};
    float tTrack_Hit_PeakTime[NMaxHitsPerEvent]={0};
    float tTrack_Hit_Amplitude[NMaxHitsPerEvent]={0};
    float tTrack_Hit_ChargeSummedADC[NMaxHitsPerEvent]={0};
    float tTrack_Hit_ChargeIntegral[NMaxHitsPerEvent]={0};
    float tTrack_Hit_StartTime[NMaxHitsPerEvent]={0};
    float tTrack_Hit_EndTime[NMaxHitsPerEvent]={0};
    float tTrack_Hit_Width[NMaxHitsPerEvent]={0};
    float tTrack_Hit_GoodnessOfFit[NMaxHitsPerEvent]={0};
    short tTrack_Hit_Multiplicity[NMaxHitsPerEvent]={0};
    int tTrack_Hit_particleID[NMaxHitsPerEvent]={0};
    float tTrack_Hit_TrueEnergy[NMaxHitsPerEvent]={0};
    float tTrack_Hit_TrueEnergyFraction[NMaxHitsPerEvent]={0};
    float tTrack_Hit_ChannelElectrons[NMaxHitsPerEvent]={0};
    float tTrack_Hit_MaxChannelElectrons[NMaxHitsPerEvent]={0};
    int tTrack_Hit_BestChanIDE[NMaxHitsPerEvent]={0};

    bool fUseChannelNoise = false;
    std::map< int, double > fChannelNoise;

    bool fActiveBranch = false; //a flag to use if the root file has been made with a specific Parser version which lives on my personal featue branch ( andrea.scarpelli@cern.ch 04/03/19 )
};

#endif // __DATASTRUCTURE_H
