////////////////////////////////////////////////////////////////////////////////
// Run class methods
//
// mailto:andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <sqlite3.h>
#include "Run.hh"

using namespace std;

Run::Run(){
  fDbName = "../metadata/test.db";
  fRunNumber = 840;
  fNumberOfSubruns = 118;
}

Run::Run(int run, string dbase){
  fDbName = dbase;
  fRunNumber = run;
  //this->importRunMetadata();
}

Run::~Run(){}

void Run::importRunMetadata(){
  //read run metadata from database
  sqlite3* fdb;
  sqlite3_stmt *statement;

  //build the sql statement query
  string kStatement = "SELECT * FROM runs WHERE run="+to_string(fRunNumber)+";";

  int res=sqlite3_open(fDbName.c_str(), &fdb);

  if (res!= SQLITE_OK){
     cout << "Run::readFromDB(): impossible to open database: " << fDbName << " error: " << sqlite3_errmsg(fdb) << endl;
     return;
   }

   if ( sqlite3_prepare(fdb, kStatement.c_str(), -1, &statement, 0 ) == SQLITE_OK ){

   int res=0;
   res = sqlite3_step(statement);
   if ( res == SQLITE_ROW ){
     fNumberOfSubruns = sqlite3_column_double(statement,1);

     cout << " Importing metadata for run: "<< fRunNumber << " from database " << fDbName << endl;
   }else{
     cout << "Unexpected sqlite3_step return value: (" <<res<<"); "<<"ERROR:"<<sqlite3_errmsg(fdb )<< endl;
   }
 }else{
   cout << "Error preparing statement: (" <<kStatement<<"); "<<"ERROR:"<<sqlite3_errmsg(fdb) << endl;
 }

 sqlite3_close(fdb);

 return;
}
