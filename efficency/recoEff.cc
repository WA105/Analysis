////////////////////////////////////////////////////////////////////////////////
//  Routine to study the tracks reconstruction efficiency
//
//  Usage: ./build/recoEff -s /path/to/detsimg4.root -r /path/to/fastreco.root -o path/to/outputfile.root
//  -s detsimg4 parser root file
//  -r fastreco parser root file
//  -o outputfile.root
//
//mailto:andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////

//c++ includes
#include <glob.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

//root includes
#include "TFile.h"
#include "TEfficiency.h"

//projects include
#include "Utils.hh"
#include "DataStructure.hh"
#include "Cuts.hh"
#include "Efficiency.hh"
#include "Geometry.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Custom functions

bool inBufferBox( MCTrack track , int buffer)
{
  //check if the mctrack starting-ending point are both too close to the tpc edges
  //Ideally we'd prefer not to count those events because they don't containt much information

  double buffMinX = tpc_boundaries[0]+buffer;
  double buffMaxX = tpc_boundaries[1]-buffer;
  double buffMinY = tpc_boundaries[2]+buffer;
  double buffMaxY = tpc_boundaries[3]-buffer;
  double buffMinZ = tpc_boundaries[4]+buffer;
  double buffMaxZ = tpc_boundaries[5]-buffer;

  bool inBufferX = ( ((track.startX < buffMinX) && ( track.endX < buffMinX)) || ((track.startX > buffMaxX) && ( track.endX > buffMaxX)) );
  bool inBufferY = ( ((track.startY < buffMinY) && ( track.endY < buffMinY)) || ((track.startY > buffMaxY) && ( track.endY > buffMaxY)) );
  bool inBufferZ = ( ((track.startZ < buffMinZ) && ( track.endZ < buffMinZ)) || ((track.startZ > buffMaxZ) && ( track.endZ > buffMaxZ)) );

  return ( inBufferX || inBufferY || inBufferZ );

}

bool isGoodMcTrack( MCTrack track )
{
  //check if the mctrack is inside the active volume for a good portion
  bool isInAv = ( track.isInTPCAV == 1 ); // first condition: be inside the tpc
  bool inBuffer = inBufferBox( track, 10 ); //is start and end contained into a buffer box?
  bool isLongEnough = ( track.length > 35 ); //third condition: have a decent length

  return ( isInAv && ( !inBuffer || isLongEnough ) );

}

////////////////////////////////////////////////////////////////////////////////
// Main macro

int main( int argc, char* argv[] ){

  //parse the input argument ===================================================
  string simFile;
  string recoFile;
  string outputFile;

  for ( int i=1; i<argc; i=i+2 ) {
   if      ( string(argv[i]) == "-s" ) simFile = argv[i+1];
   else if ( string(argv[i]) == "-r" ) recoFile = argv[i+1];
   else if ( string(argv[i]) == "-o" ) outputFile = argv[i+1];
   else {
     cout << "Unknown option" << endl;
     return 1;
   }
  }

  //define and variables =======================================================

  //extract the filenumber from the simFile
  //int fileNumber = getFileNumber( simFile );

  //define here the output file
  Track recoTrack;
  MCTrack mcTrack;

  //define output file

  TFile *ofile = new TFile(outputFile.c_str(), "RECREATE");
  if(!ofile->IsOpen()){
    cout << "File: " << outputFile << " cannot be open!" << endl;
    return 1;
  }

  //Prepare inputs =============================================================

  LArParser *mcParser = new LArParser();
  LArParser *recoParser = new LArParser();

  int fileNum = stoi( getFileNumber( simFile ) );

  //and here i define the class efficiency
  Efficiency *recoEfficiency = new Efficiency( fileNum );

  TTree *mcTree = getTTree( simFile );
  TTree *recoTree = getTTree( recoFile );

  //check if the tree has been correctly set ===================================

  if( !mcTree || !recoTree ){
   cout << "Error! Trees doesn't exist! " << endl;
   return 1;
  }

  //check if the two trees have the same number of entries.
  //I ideally want to loop over one of them
  if( mcTree->GetEntries() != recoTree->GetEntries() ){
    cout << "Not the same number of events" << endl;
    return 1;
  }

  //Event looper ===============================================================
  for(int evt=0; evt<mcTree->GetEntries(); evt++){

      cout << "Processing event: " << evt << endl;

      //data structures array
      vector<MCTrack> mcTracks;
      vector<Track> recoTracks;
      vector<Hit> recoHits; //hits not associated to a track

      mcParser->getMCTracksEvent(mcTree, mcTracks, evt);
      recoParser->getRecoTracksEvent(recoTree, recoTracks, evt);
      recoParser->getRecoHitsEvent( recoTree, recoHits, evt );

      //set recoHits
      recoEfficiency->setRecoHits( recoHits );

      for( auto track : mcTracks ){
        //order truth track into a map sorted by their id so is easier to make
        //a match between true and reco

        if ( !isGoodMcTrack( track ) ) { continue; } //jump tracks which aren't inside the active volume

        recoEfficiency->setMapEntry( track.particleID, track );
      }

      //insert the reconstruced hits inside the event
      recoEfficiency->setNumberOfTracksEvent( (int)recoTracks.size() );

      for( auto track : recoTracks ){
        //here we do the real efficiency analysis trying to match reconstructed
        //calculate efficiency for the track and fill the historams
        recoEfficiency->setRecoTrack( track );
        recoEfficiency->fill();
      }

      recoEfficiency->checkUnmatch(); //build up a ttree with all the unmatched true particles in the event
      recoEfficiency->clean();

  }//end event loop

  //calculate the root TEfficency plots
  recoEfficiency->makeEfficiencyPlot();

  ofile->cd();
  recoEfficiency->write();

  ofile->Close();

}//end macro
