////////////////////////////////////////////////////////////////////////////////
// Macro studying the reconstruction efficiency of one MC sample
//
//mailto:andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////

//TODO dynamic path to input and output file

#include <glob.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h> //to use struct stat

#include "TChain.h"
#include "Run.hh"
#include "DataStructure.hh"
#include "Cuts.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//initalize some variables as global variables (easy to edit)

vector<int> fileList = {}; //runs to process
int startFile = 0;
int endFile = 500;

int mockRun = 840; //query the metadata of this run from db

////////////////////////////////////////////////////////////////////////////////
// Functions

inline bool ExistTest (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

inline vector<string> glob(const string& pat){
  glob_t glob_result;
  glob(pat.c_str(),GLOB_TILDE,NULL,&glob_result);
  vector<string> ret;
  for(unsigned int i=0;i<glob_result.gl_pathc;++i){
      string full_path = string(glob_result.gl_pathv[i]);
      string basename = full_path.substr(full_path.find_last_of("/")+1);
      ret.push_back(basename.data());
  }
  globfree(&glob_result);
  return ret;
}

void fillFileList(string path){
  //fill up the runList with all the file if empty

  for(i=startFile; i<endFile; i++){
    fileList.push_back(i);
  }

  /*
  if (fileList.size() == 0){
    #if verbose
    cout << "Processing all runs in " << path << "*..." << endl;
    #endif
    string wildcard_path = path + "*";
    for( auto irun : glob(wildcard_path) ){
      fileList.push_back(atoi(irun.data()));
    }
  }
  */

  return;
}

TTree *getTTree( string path, string prefix, string suffix, int runNumber ){

  TTree *ttree;

  //no need to run over subruns (are mc). Just jump to filename definition
  string file=prefix+to_string(runNumber)+suffix;
  string filename=path+file;

  if( ExistTest(filename) ){

      cout << "Processing file: " << filename << endl;
      TFile *file = new TFile(filename.c_str(), "READ");

      if( file->IsOpen() )
        ttree = (TTree*)file->Get("analysistree/anatree");
      else
        cout << "getTTree::Error: Invalid TTree name" << endl;

    }else{
      cout << "getTTree::Error: Invalid file name " << filename << endl;
    }

    return ttree;
}

////////////////////////////////////////////////////////////////////////////////
// Reco efficiency class definition and associated methods
// Handles efficiency calculation and the relevant plots

class Efficiency
{
  public:
    Efficiency();
    ~Efficiency();

    //setters
    void setMapEntry(int id, MCTrack mctrack ){ fParticleMap[id] = mctrack; }
    void setRecoTrack( Track track ){ fTrack = track; }
    void setRecoHits( vector<Hit> hits ){ fHits = hits; }

    //getters

    //others
    void fill();
    void write(TFile *ofile);

    //cleaner
    void clean();

  private:
    void matchTruth();
    void fillMap1D(int pdg, map<int, TH1D*> map, double fillIn );
    void fillMap2D(int pdg, map<int, TH2D*> map, double fillX, double fillY );

    //Particle that I consider for the efficiency
    vector<int> pdgCode = { 13, -13, 11, -11, 211, -211, 2212, 0 };
    vector<string> pdgNames = { "Muons", "Antimuons", "Electrons", "Positron", "PiPlus", "PiMinus", "Proton", "Other" };

    //binning
    //theta
    int nBinsTheta = 100;
    double thetaStart = 0;
    double thetaEnd = 180;

    //phi
    int nBinsPhi = 100;
    double phiStart = -180;
    double phiEnd = 180;

    //histogram maps
    map<int, TH1D*> fThetaTrueMap;
    map<int, TH1D*> fPhiTrueMap;
    map<int, TH2D*> fPhiThetaTrueMap;
    map<int, TH1D*> fThetaRecoMap;
    map<int, TH1D*> fPhiRecoMap;
    map<int, TH2D*> fPhiThetaRecoMap;

    //more quantities
    map<int, MCTrack> fParticleMap; //particleID mctrack association
    map<int, double> fEnergyMap; //particleID energy association
    Track fTrack;
    vector<Hit> fHits;
    double fPurirty;
    double fCompleteness;
    double fPdg;
    double fBestTrackID;

};

Efficiency::Efficiency(){
  //constructor of the class: initialize here all the histograms

  size_t arraySize = max(pdgCode.size(), pdgNames.size() );

  for( size_t i =0; i< arraySize; i++ ){

    TH1D* hThetaTrue = new TH1D(("hThetaTrue"+pdgNames.at(i)).c_str(), ("#theta (true) "+pdgNames.at(i)+";#theta (deg)").c_str() , nBinsTheta, thetaStart, thetaEnd);
    TH1D* hPhiTrue = new TH1D(("hPhiTrue"+pdgNames.at(i)).c_str(), ("#phi (true) "+pdgNames.at(i)+";#phi (deg)").c_str() , nBinsPhi, phiStart, phiEnd);
    TH2D* hPhiThetaTrue = new TH2D(("hPhiThetaTrue"+pdgNames.at(i)).c_str(), ("#phi vs. #theta (true) "+pdgNames.at(i)+"; #theta (deg);#phi (deg)").c_str(), nBinsTheta, thetaStart, thetaEnd, nBinsPhi, phiStart, phiEnd); ;

    TH1D* hThetaReco = new TH1D(("hThetaReco"+pdgNames.at(i)).c_str(), ("#theta (reco) "+pdgNames.at(i)+";#theta (deg)").c_str() , nBinsTheta, thetaStart, thetaEnd);;
    TH1D* hPhiReco = new TH1D(("hPhiReco"+pdgNames.at(i)).c_str(), ("#phi (reco) "+pdgNames.at(i)+";#phi (deg)").c_str() , nBinsPhi, phiStart, phiEnd);
    TH2D* hPhiThetaReco = new TH2D(("hPhiThetaReco"+pdgNames.at(i)).c_str(), ("#phi vs. #theta (true) "+pdgNames.at(i)+"; #theta (deg);#phi (deg)").c_str(), nBinsTheta, thetaStart, thetaEnd, nBinsPhi, phiStart, phiEnd);

    fThetaTrueMap[pdgCode.at(i)] = hThetaTrue;
    fThetaRecoMap[pdgCode.at(i)] = hThetaReco;
    fPhiTrueMap[pdgCode.at(i)] = hPhiTrue;
    fPhiRecoMap[pdgCode.at(i)] = hPhiReco;
    fPhiThetaTrueMap[pdgCode.at(i)] = hPhiThetaTrue;
    fPhiThetaRecoMap[pdgCode.at(i)] = hPhiThetaReco;

  }
}

Efficiency::~Efficiency(){}

void Efficiency::matchTruth(){
  //match reco and truth, calculate purity and completeness of a reco track

  //first of all reset all the variables
  fEnergyMap.clear();
  fPurirty=0;
  fCompleteness=0;
  fPdg=0;

  //loop over the hits in tracks and associate every energy deposit to the correct particleID
  double energyTrk = 0.;
  for( auto hit : fTrack.hitsTrk ){
    fEnergyMap[ hit.particleID ] += hit.trueEnergy;
    energyTrk += (hit.trueEnergy/hit.trueEnergyFraction);
  }

  //find the best particle ID (the one that contribute the most in the track energy account )
  fBestTrackID =0.;
  double maxe = 0.;
  for( auto const & val : fEnergyMap  ){
    if( maxe < val.second ){
      maxe = val.second;
      fBestTrackID= val.first;
    }
  }

  //calculate the total energy of the best track in hits
  double totalEnergy = 0;
  for( auto hit : fHits ){
    if( fBestTrackID == hit.particleID )
      totalEnergy+=hit.trueEnergy;
  }

  fPurirty = (fEnergyMap[ fBestTrackID ])/totalEnergy;
  fCompleteness = fEnergyMap[ fBestTrackID ]/energyTrk;
  fPdg = fParticleMap[ fBestTrackID ].pdgCode;

}

void Efficiency::fillMap1D(int pdg, map<int, TH1D*> map, double fillIn ){
  //fill the map if the pdg code of the best tParticleId
    if( map.find(pdg) != map.end() )
      map[pdg]->Fill( fillIn );
    else
      map[0]->Fill( fillIn );
}

void Efficiency::fillMap2D(int pdg, map<int, TH2D*> map, double fillX, double fillY ){
  //fill the map if the pdg code of the best tParticleId
  if( map.find(pdg) != map.end() )
    map[pdg]->Fill( fillX, fillY );
  else
    map[0]->Fill( fillX, fillY );
}

void Efficiency::fill(){
  //fill up the correct histogram for the track

  this->matchTruth();

  //fill first the mc quanties
  fillMap1D( fPdg, fThetaTrueMap, fParticleMap[fBestTrackID].startTheta);
  fillMap1D( fPdg, fPhiTrueMap, fParticleMap[fBestTrackID].startPhi );
  fillMap2D( fPdg, fPhiThetaTrueMap, fParticleMap[fBestTrackID].startTheta, fParticleMap[fBestTrackID].startPhi);

  if(fCompleteness>0.5 && fPurirty>0.5){

    //fill the reco quantities
    fillMap1D( fPdg, fThetaRecoMap, fParticleMap[fBestTrackID].startTheta);
    fillMap1D( fPdg, fPhiRecoMap, fParticleMap[fBestTrackID].startPhi );
    fillMap2D( fPdg, fPhiThetaRecoMap, fParticleMap[fBestTrackID].startTheta, fParticleMap[fBestTrackID].startPhi);
  }
}

void Efficiency::write(TFile *ofile){
  //write histogram in root file

  ofile->cd();
  size_t arraySize = max(pdgCode.size(), pdgNames.size() );

  for( size_t i =0; i< arraySize; i++ ){

    fThetaTrueMap[pdgCode.at(i)]->Write();
    fThetaRecoMap[pdgCode.at(i)]->Write();
    fPhiTrueMap[pdgCode.at(i)]->Write();
    fPhiRecoMap[pdgCode.at(i)]->Write();
    fPhiThetaTrueMap[pdgCode.at(i)]->Write();
    fPhiThetaRecoMap[pdgCode.at(i)]->Write();

  }

}

void Efficiency::clean(){

  fParticleMap.clear();
  fHits.clear();
}

////////////////////////////////////////////////////////////////////////////////
// Main macro

void recoEff(){
  //Example of a macro to study the reconsturction eff.

  //define here the output file
  Track recoTrack;
  MCTrack mcTrack;

  TFile *ofile = new TFile("simOutput.root", "RECREATE");

  //mc tree
  //TTree *mcOutputTree = new TTree("mcTracks", "contains geant tracks");
  //mcOutputTree->Branch("mcTrack", &mcTrack );
  //reco tree
  //TTree *recoOutputTree = new TTree("recoTracks", "contains reco tracks");
  //recoOutputTree->Branch("recoTrack", &recoTrack );

  //here I define the parser object
  LArParser *mcParser = new LArParser();
  LArParser *recoParser = new LArParser();

  //define the Run object using the mockRun flag (always the same in this case)
  Run *run = new Run(mockRun, "metadata/test.db");

  //fill the filelist if empty
  //fillFileList();

  //and here i define the class efficiency
  Efficiency *recoEfficiency = new Efficiency();

  //loop over all the runs
  for(int fileNumber : fileList){

     TTree *mcTree = getTTree("/Users/scarpell/cernbox/311/simulation/g4detsim/", "", "-G4Detsim-Parser.root", fileNumber ); //TODO dynamic path
     TTree *recoTree = getTTree("/Users/scarpell/cernbox/311/simulation/ana/", "", "-RecoFull-Parser.root", fileNumber ); //TODO: dyynamic path

     mcParser->setTTree(mcTree);
     mcParser->setRun(run);

     recoParser->setTTree(recoTree);
     recoParser->setRun(run);

     //check if the tree exists and has been correctly set
     if( !mcParser->isTreeGood() || !recoParser->isTreeGood() ){
       cout << "Invalid ttree" << endl;
       continue;
     }

    //check if the two trees have the same number of entries. I ideally want to loop over one of them
    if( mcTree->GetEntries() != recoTree->GetEntries() ){
      cout << "Not the same number of events" << endl;
      continue;
    }

    //loop over the events in mcTree. Should be the same for also recoTree
    for(int evt=0; evt<mcTree->GetEntries(); evt++){

      //data structures array
      vector<MCTrack> mcTracks;
      vector<Track> recoTracks;
      vector<Hit> recoHits; //hits not associated to a track

      mcParser->getMCTracksEvent(mcTracks, evt);
      recoParser->getRecoTracksEvent(recoTracks, evt);
      recoParser->getRecoHitsEvent( recoHits, evt );

      for( auto track : mcTracks ){

        //order truth track into a map sorted by their id so is easier to make
        //a match between true and reco
        recoEfficiency->setMapEntry( track.particleID, track );

        //mcTrack = track;
        //mcOutputTree->Fill();
      }

      //insert the reconstruced hits inside the event
      recoEfficiency->setRecoHits( recoHits );

      for( auto track : recoTracks ){

        //here we do the real efficiency analysis trying to match reconstructed
        //calculate efficiency for the track and fill the historams
        recoEfficiency->setRecoTrack( track );
        recoEfficiency->fill();

        //recoTrack = track;
        //recoOutputTree->Fill();
      }

      //recoEfficiency->clean(); //clean the trueParticle map inside the class
    }//end event loop
  }//end filelist run


  ofile->cd();
  //write efficiency histograms
  recoEfficiency->write( ofile );

  //write ttree
  //mcOutputTree->Write();
  //recoOutputTree->Write();
  ofile->Close();

}//end macro
