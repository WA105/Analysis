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

class Run
{

  public:
    Run();
    Run(int run, string dbase);
    ~Run();

    //setters
    //void setDBName( string dbName ){ fDbName = dbName; };
    void setRunNumber( int run ){ fRunNumber = run; };
    void setNumberOfSubrun( int nSubruns ){ fNumberOfSubruns = nSubruns; };

    //getters
    int getRunNumber(){ return fRunNumber; }
    int getNumberOfSubruns(){ return fNumberOfSubruns; }

  private:

    //void importRunMetadata();

    string fDbName;
    int fRunNumber;
    int fNumberOfSubruns;

};

#endif //__RUN_H
