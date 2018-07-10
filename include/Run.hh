////////////////////////////////////////////////////////////////////////////////
// Run class definition
//
// mailto:andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////

#ifndef  __RUN_H
#define __RUN_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

class Run
{

  public:
    Run();
    Run(int run, string dbase);
    ~Run();

    //setters
    void setDBName( string dbName ){ fDbName = dbName; };
    void setRunNumber( int run ){ fRunNumber = run; };
    void setFileNumber( int fileNumber ){ fFileNumber = fileNumber; }; //montecarlo doesn't have an unique run number
    void setNumberOfSubrun( int nSubruns ){ fNumberOfSubruns = nSubruns; };

    //getters
    int getRunNumber(){ return fRunNumber; }
    int getNumberOfSubruns(){ return fNumberOfSubruns; }
    int getFileNumber(){ return fFileNumber; };

  private:

    void importRunMetadata();

    string fDbName;
    int fRunNumber;
    int fNumberOfSubruns;
    int fFileNumber;

};

#endif //__RUN_H
