////////////////////////////////////////////////////////////////////////////////
//  Routine to study a run signal-to-noise ratio
//
//  Usage: ./build/signalToNoise -s /path/to/data.root -c /path/to/cuts.txt -n /path/to/noise_model.root -o path/to/outputFile.root
//  -s FastReco data runs
//  -c List of tracks which pass a given cut calculated externally (e.g. highway cut)
//  -r Type of runs: Data or Montecarlo
//  -n File .root with the reference noise model
//  -o outputfile.root
//
//  mailto:andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////

//c++ includes
#include <glob.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <numeric>

//root includes
#include "TFile.h"
#include "TTree.h"

//projects include
#include "Utils.hh"
#include "DataStructure.hh"
#include "Cuts.hh"
#include "Efficiency.hh"
#include "Geometry.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////////////

class MakeFile
{
  public:
    MakeFile();
    ~MakeFile();

    void setRunType( string runType ){ fRunType = runType; };
    void fillTrack( Track t );
    void fillHit( Hit h );
    void Write();

  private:

    //general
    string fRunType = "Data";

    //root
    TTree *fTrackTree;
    TTree *fHitTree;

    //event specific
    int fRun;
    int fSubrun;
    int fEvent;
    int fNumerOfHits;

    //track specific
    int fEventID=0;
    int fTrackId;
    double fSDirX;
    double fSDirY;
    double fSDirZ;
    double fEDirX;
    double fEDirY;
    double fEDirZ;
    double fSPointX;
    double fSPointY;
    double fSPointZ;
    double fEPointX;
    double fEPointY;
    double fEPointZ;
    double fStartTheta;
    double fStartPhi;
    double fLength;
    int fNumberOfHitsView0;
    int fNumberOfHitsView1;
    double fAvgWidthView0;
    double fAvgWidthView1;
    double fAvgAmplitudeView0;
    double fAvgAmplitudeView1;
    double fAvgSummedADCView0;
    double fAvgSummedADCView1;

    //hit specific
    int fView;
    int fChannel;
    double fPeakTime;
    double fIntegral;
    double fSummedADC;
    double fAmplitude;
    double fWidth;
    double fChi2;
    double fPeakToNoise;
    double fAreaToNoise;
    double fLocalGamma; //Angle between wire and time direction
    double fTheta;
    double fPhi;
    double fX;
    double fY;
    double fZ;
    double fLocaldS;
    double f3DPositiondS;
    double fdQds;

    //functions
    int getNumOfHits( vector<Hit> hits, int view );
    double getAvgWidth( vector<Hit> hits, int view );
    double getAvgAmplitude( vector<Hit> hits, int view );
    double getAvgSummedADC( vector<Hit> hits, int view );
};

MakeFile::MakeFile()
{
  fTrackTree = new TTree("track", "Track information");
  fHitTree = new TTree("hits", "Hits information");

  //track tree
  fTrackTree->Branch("Run", &fRun, "Run/I");
  fTrackTree->Branch("Subrun", &fSubrun, "Subrun/I");
  fTrackTree->Branch("Event", &fEvent, "Event/I");
  fTrackTree->Branch("Track", &fTrackId, "Track/I");
  fTrackTree->Branch("SDirX", &fSDirX, "SDirX/D");
  fTrackTree->Branch("SDirY", &fSDirY, "SDirY/D");
  fTrackTree->Branch("SDirZ", &fSDirZ, "SDirZ/D");
  fTrackTree->Branch("EDirX", &fEDirX, "EDirX/D");
  fTrackTree->Branch("EDirY", &fEDirY, "EDirY/D");
  fTrackTree->Branch("EDirZ", &fEDirZ, "EDirZ/D");
  fTrackTree->Branch("Length", &fLength, "Length/D");
  fTrackTree->Branch("SPointX", &fSPointX, "SPointX/D");
  fTrackTree->Branch("SPointY", &fSPointY, "SPointY/D");
  fTrackTree->Branch("SPointZ", &fSPointZ, "SPointZ/D");
  fTrackTree->Branch("EPointX", &fEPointX, "EPointX/D");
  fTrackTree->Branch("EPointY", &fEPointY, "EPointY/D");
  fTrackTree->Branch("EPointZ", &fEPointZ, "EPointZ/D");
  fTrackTree->Branch("StartTheta", &fStartTheta, "StartTheta/D");
  fTrackTree->Branch("StartPhi", &fStartPhi, "StartPhi/D");
  fTrackTree->Branch("NumberOfHitsView0", &fNumberOfHitsView0, "NumberOfHitsView0/I");
  fTrackTree->Branch("NumberOfHitsView1", &fNumberOfHitsView1, "NumberOfHitsView1/I");
  fTrackTree->Branch("AvgAmplitudeView0", &fAvgAmplitudeView0, "AvgAmplitudeView0/D");
  fTrackTree->Branch("AvgAmplitudeView1", &fAvgAmplitudeView1, "AvgAmplitudeView1/D");
  fTrackTree->Branch("AvgWidthView0", &fAvgWidthView0, "AvgWidthView0/D");
  fTrackTree->Branch("AvgWidthView1", &fAvgWidthView1, "AvgWidthView1/D");
  fTrackTree->Branch("AvgSummedADCView0", &fAvgSummedADCView0, "AvgSummedADCView0/D");
  fTrackTree->Branch("AvgSummedADCView1", &fAvgSummedADCView1, "AvgSummedADCView1/D");

  //hittree
  fHitTree->Branch("Run", &fRun, "Run/I");
  fHitTree->Branch("Subrun", &fSubrun, "Subrun/I");
  fHitTree->Branch("Event", &fEvent, "Event/I");
  fHitTree->Branch("Track", &fTrackId, "Track/I");
  fHitTree->Branch("View", &fView, "View/I");
  fHitTree->Branch("Channel", &fChannel, "Channel/I");
  fHitTree->Branch("PeakTime", &fPeakTime, "PeakTime/D");
  fHitTree->Branch("Integral", &fIntegral, "Integral/D");
  fHitTree->Branch("SummedADC", &fSummedADC, "SummedADC/D");
  fHitTree->Branch("fLocaldS", &fLocaldS, "fLocaldS/D");
  fHitTree->Branch("f3DPositiondS", &f3DPositiondS, "f3DPositiondS/D");
  fHitTree->Branch("dQds", &fdQds, "dQds/D" );
  fHitTree->Branch("Amplitude", &fAmplitude, "Amplitude/D");
  fHitTree->Branch("Width", &fWidth, "Width/D");
  fHitTree->Branch("Chi2", &fChi2, "Chi2/D");
  fHitTree->Branch("PeakToNoise", &fPeakToNoise, "PeakToNoise/D");
  fHitTree->Branch("AreaToNoise", &fAreaToNoise, "AreaToNoise/D");
  fHitTree->Branch("LocalTheta", &fTheta, "LocalTheta/D");
  fHitTree->Branch("LocalPhi", &fPhi, "LocalPhi/D");
  fHitTree->Branch("X", &fX, "X/D" );
  fHitTree->Branch("Y", &fY, "Y/D" );
  fHitTree->Branch("Z", &fZ, "Z/D" );

  //add information on the signle hit of the track he belong to
  fHitTree->Branch("SDirX", &fSDirX, "SDirX/D");
  fHitTree->Branch("SDirY", &fSDirY, "SDirY/D");
  fHitTree->Branch("SDirZ", &fSDirZ, "SDirZ/D");
  fHitTree->Branch("EDirX", &fEDirX, "EDirX/D");
  fHitTree->Branch("EDirY", &fEDirY, "EDirY/D");
  fHitTree->Branch("EDirZ", &fEDirZ, "EDirZ/D");
  fHitTree->Branch("Length", &fLength, "Length/D");
  fHitTree->Branch("SPointX", &fSPointX, "SPointX/D");
  fHitTree->Branch("SPointY", &fSPointY, "SPointY/D");
  fHitTree->Branch("SPointZ", &fSPointZ, "SPointZ/D");
  fHitTree->Branch("EPointX", &fEPointX, "EPointX/D");
  fHitTree->Branch("EPointY", &fEPointY, "EPointY/D");
  fHitTree->Branch("EPointZ", &fEPointZ, "EPointZ/D");
  fHitTree->Branch("StartTheta", &fStartTheta, "StartTheta/D");
  fHitTree->Branch("StartPhi", &fStartPhi, "StartPhi/D");
  fHitTree->Branch("NumberOfHitsView0", &fNumberOfHitsView0, "NumberOfHitsView0/I");
  fHitTree->Branch("NumberOfHitsView1", &fNumberOfHitsView1, "NumberOfHitsView1/I");
  fHitTree->Branch("AvgAmplitudeView0", &fAvgAmplitudeView0, "AvgAmplitudeView0/D");
  fHitTree->Branch("AvgAmplitudeView1", &fAvgAmplitudeView1, "AvgAmplitudeView1/D");
  fHitTree->Branch("AvgWidthView0", &fAvgWidthView0, "AvgWidthView0/D");
  fHitTree->Branch("AvgWidthView1", &fAvgWidthView1, "AvgWidthView1/D");
  fHitTree->Branch("AvgSummedADCView0", &fAvgSummedADCView0, "AvgSummedADCView0/D");
  fHitTree->Branch("AvgSummedADCView1", &fAvgSummedADCView1, "AvgSummedADCView1/D");

}

MakeFile::~MakeFile(){}

int MakeFile::getNumOfHits( vector<Hit> hits, int view )
{
  int numOfHits = 0;
  auto lambda = [ view ]( int accumulator, Hit h )
  {
     return h.view == view ? accumulator + 1 : accumulator;
  };
  numOfHits = accumulate( hits.begin(), hits.end(), 0, lambda );

  return numOfHits;
}

double MakeFile::getAvgWidth( vector<Hit> hits, int view )
{
  double sum = 0.;
  auto lambda = [ view, &sum ]( double accumulator, Hit hit )
  {
     if( hit.view == view && hit.width > 0 )
     {
       ++sum;
       return accumulator + hit.width;
     }
     else
     {
       return accumulator;
     }
  };
  double accumulated = accumulate( hits.begin(), hits.end(), 0.0, lambda );
  return sum > 0. ? accumulated/sum : 0.;
}

double MakeFile::getAvgAmplitude( vector<Hit> hits, int view )
{
  double sum = 0.;
  auto lambda = [ view, &sum ]( double accumulator, Hit hit )
  {
     if( hit.view == view && hit.amplitude > 0 )
     {
       ++sum;
       return accumulator + hit.amplitude;
     }
     else
     {
       return accumulator;
     }
  };
  double accumulated = accumulate( hits.begin(), hits.end(), 0.0, lambda );
  return sum > 0. ? accumulated/sum : 0.;
}

double MakeFile::getAvgSummedADC( vector<Hit> hits, int view )
{
  double sum = 0.;
  auto lambda = [ view, &sum ]( double accumulator, Hit hit )
  {
     if( hit.view == view && hit.chargeSummedADC > 0 )
     {
       ++sum;
       return accumulator + hit.chargeSummedADC;
     }
     else
     {
       return accumulator;
     }
  };
  double accumulated = accumulate( hits.begin(), hits.end(), 0.0, lambda );
  return sum > 0. ? accumulated/sum : 0.;
}

void MakeFile::fillTrack( Track t )
{
  //Wrapper to fill the track ttree
  fEventID += 1;
  fRun = t.run;
  fSubrun = t.subRun;
  fEvent = t.event;
  fTrackId = t.trackID;

  if ( max( t.startPointX, t.endPointX ) == t.startPointX )
  {
    fSDirX = t.startDirectionX;
    fSDirY = t.startDirectionY;
    fSDirZ = t.startDirectionZ;
    fEDirX = t.endDirectionX;
    fEDirY = t.endDirectionY;
    fEDirZ = t.endDirectionZ;
    fSPointX = t.startPointX;
    fSPointY = t.startPointY;
    fSPointZ = t.startPointZ;
    fEPointX = t.endPointX;
    fEPointY = t.endPointY;
    fEPointZ = t.endPointZ;
  }
  else if ( max( t.startPointX, t.endPointX ) == t.endPointX ) //flip start-end
  {
    fSDirX = t.endDirectionX;
    fSDirY = t.endDirectionY;
    fSDirZ = t.endDirectionZ;
    fEDirX = t.startDirectionX;
    fEDirY = t.startDirectionY;
    fEDirZ = t.startDirectionZ;
    fSPointX = t.endPointX;
    fSPointY = t.endPointY;
    fSPointZ = t.endPointZ;
    fEPointX = t.startPointX;
    fEPointY = t.startPointY;
    fEPointZ = t.startPointZ;
  }
  else{ cout << "fillTrack: invalid maximum!" << endl; }

  fLength = t.length;
  fStartTheta = t.startDirectionTheta;
  fStartPhi =t.startDirectionPhi;
  fNumberOfHitsView0 = this->getNumOfHits( t.hitsTrk, 0 );
  fNumberOfHitsView1 = this->getNumOfHits( t.hitsTrk, 1 );
  fAvgWidthView0 = this->getAvgWidth( t.hitsTrk, 0 );
  fAvgWidthView1 = this->getAvgWidth( t.hitsTrk, 1 );
  fAvgAmplitudeView0 = this->getAvgAmplitude( t.hitsTrk, 0 );
  fAvgAmplitudeView1 = this->getAvgAmplitude( t.hitsTrk, 1 );
  fAvgSummedADCView0 = this->getAvgSummedADC( t.hitsTrk, 0 );
  fAvgSummedADCView1 = this->getAvgSummedADC( t.hitsTrk, 1 );
  fTrackTree->Fill();
}

void MakeFile::fillHit( Hit h )
{
  //Wrapper to fill the hit ttree
  double adc2fc = getCalAmpConstant(h.view, fRunType);

  fRun = h.run;
  fSubrun = h.subRun;
  fEvent = h.event;
  fView = h.view;
  fTrackId = h.trackID;
  fChannel = h.channel;
  fPeakTime = h.peakTime;
  fAmplitude = h.amplitude;
  fWidth = h.width;
  fIntegral = h.chargeIntegral*adc2fc; //in fC
  fSummedADC = h.chargeSummedADC*adc2fc; //in fC
  fPeakToNoise = h.getPeakToNoise();
  fAreaToNoise = h.getIntegralToNoise(fRunType);
  fLocaldS = h.dxLocalTrackDirection ;
  f3DPositiondS = h.dx3DPosition;
  if ( fLocaldS>0 ){
    fdQds = fSummedADC/f3DPositiondS;
  }
  else{
    fdQds = 0;
  }
  fPhi = h.phi;
  fTheta = h.theta;
  fX = h.X;
  fY = h.Y;
  fZ = h.Z;
  fChi2 = h.goodnessOfFit;

  fHitTree->Fill();
}

void MakeFile::Write()
{
  fTrackTree->Write();
  fHitTree->Write();
}

////////////////////////////////////////////////////////////////////////////////
// Main macro

int main( int argc, char* argv[] ){

  //gROOT->ProcessLine( "gErrorIgnoreLevel = 3001;")

  //parse the input argument ===================================================

  string cutsFile;
  string recoFile;
  string noiseFile;
  string outputFile;
  string runType;

  for ( int i=1; i<argc; i=i+2 ) {
   if      ( string(argv[i]) == "-s" ) recoFile = argv[i+1];
   else if ( string(argv[i]) == "-c" ) cutsFile = argv[i+1];
   else if ( string(argv[i]) == "-r" ) runType = argv[i+1];
   else if ( string(argv[i]) == "-n" ) noiseFile = argv[i+1];
   else if ( string(argv[i]) == "-o" ) outputFile = argv[i+1];
   else {
     cout << "Unknown option" << endl;
     return 1;
   }
  }

  cout << runType << endl;

  if( runType != "Montecarlo" && runType != "Data" )
  {
    cout << "main(): invalid runType options."
        << "Possible options are 'Data' or 'Montecarlo'" << endl;
    return 1;
  }

  cout << "Processing file: " << recoFile << " cuts: " << cutsFile << endl;

  //define LArParser  ==========================================================

  LArParser *recoParser = new LArParser();
  recoParser->setChannelNoise( noiseFile );

  //define output file and TTree address =======================================

  TFile *ofile = new TFile(outputFile.c_str(), "RECREATE");
  if(!ofile->IsOpen())
  {
    cout << "File: " << outputFile << " cannot be open!" << endl;
    return 1;
  }

  MakeFile *myFile = new MakeFile();
  myFile->setRunType(runType);

  //set the TTree  =============================================================

  TTree *recoTree = getTTree( recoFile );

  if( !recoTree ){
   cout << "Error! Trees doesn't exist! " << endl;
   return 1;
  }

  //Here define the required cuts on the track =================================

  Cut *myCuts = new Cut();
  vector<int> activeVolumeCut = { 10, 10, 10, 10, 10, 10 };
  myCuts->setActiveVolumeCut( activeVolumeCut );

  if( runType == "Montecarlo" ){ myCuts->appyNoCuts(); }
  if( runType == "Data" )
  {
    //if ( !myCuts->cutsFromHighway(cutsFile, 0.1, 9999) ){
    if(!myCuts->cutsFromFile(cutsFile)){
      cout << "ERROR! " << cutsFile << "does not exist! No cut applied" << endl;
    }
  }

  //Event looper ===============================================================
  for(int evt=0; evt<recoTree->GetEntries(); evt++)
  {

      cout << "Processing event: " << evt << endl;

      //data structures array
      vector<Track> recoTracks;
      vector<Hit> recoHits; //hits not associated to a track

      recoParser->getRecoTracksEvent(recoTree, recoTracks, evt);

      //loop over tracks
      for( Track track : recoTracks )
      {

        if( myCuts->isPassingCut( track ) && track.length > 50) {

          if( track.startDirectionTheta > 92 && track.startDirectionTheta < 160
                && track.startDirectionPhi > 25 && track.startDirectionPhi < 65 )
          {

            //loop over hits and sort the hits per views
            map< int, vector<Hit> > fViewHitMap;
            for( Hit hit : track.hitsTrk ){ fViewHitMap[hit.view].push_back(hit); }

            //ensure to have enough hits per view to remove the first and last 3
            if ( fViewHitMap[0].size() < 7 || fViewHitMap[1].size() <7  ) { continue; };

            //fill the track ttree
            myFile->fillTrack( track );

            //now loop over the map entries and remove some hits:
            for( auto viewList : fViewHitMap )
            {
              //loop only between the second and the second to last hit of the collection
              for( vector<Hit>::iterator hitList = viewList.second.begin()+3;
                                     hitList != viewList.second.end()-3; hitList++ )
              {

                Hit hit = *hitList;

                // 1. remove hits were loose wires are and under troublesome lems
                bool condRegion = ( hit.Z > 52 && hit.Z < 248 );

                // 2. remove hits with unphysical ds
                bool condPitch = ( hit.dx3DPosition < 0.3125 );

                // 3. remove hits with multiplicity and chi2 too large
                bool condFit = ( hit.multiplicity > 1 && ( hit.goodnessOfFit ==0 || hit.goodnessOfFit > 5)  );

                  // 4. remove problematic channles
                 bool condChannel = ( hit.channel > 576 && hit.channel < 607 );

                 // 5. remove unphisical LocalTheta ( horizontal )
                 bool condLocalTheta = ( hit.theta < 92 );

                if ( condRegion && !condPitch && !condFit && !condChannel && !condLocalTheta )
                {
                  myFile->fillHit( hit );
                }
              }//end for hits
            } //end view loops

            fViewHitMap.clear();

          } //end if angle
        } //end if myCuts
      }//end tracks
  }//end event loop

  //wite in ofile
  ofile->cd();
  myFile->Write();
  ofile->Close();

}//end macro
