#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "Geometry.hh"

using namespace std;

unsigned int ViewToDAQChan(unsigned int ViewChan){
  //convert View number to daq channel number

  size_t crate = ViewChan / 320;
  size_t Chan311;

  ViewChan = 8*(ViewChan/8+1)-ViewChan%8 -1;

  if(crate == 0)
  {
    ViewChan = 32*(ViewChan/32+1)-ViewChan%32 -1;
    size_t card = 4 - ((ViewChan / 32) % 5);
    if(ViewChan > 159)
    {
        size_t shift = 31 - (ViewChan % 32);
        Chan311 = (2*card)*32 + shift;
    }
    else
    {
       size_t shift = 31 - (ViewChan % 32);
       Chan311 = (2*card + 1)*32 + shift;
    }
  }
  else
  {
     size_t new_ViewChan = ViewChan - crate*320;
     size_t card = ((new_ViewChan / 32) % 5);
     if(new_ViewChan > 159)
     {
        size_t shift = new_ViewChan % 32;
        Chan311 = (2*card)*32 + shift;
     }
     else
     {
       size_t shift = new_ViewChan % 32;
       Chan311 = (2*card + 1)*32 + shift;
     }
     Chan311 = Chan311 + crate*320;
  } // end of if/else statementi

  return Chan311;
}


/*
int find_lem(double y, double z){
  //for the 311 geometry retun the lem nubmer associated to a certain geometry
  int lem =-99999;
  int lem_in_module = 4; //num of minimal modules of 4 lems

 //fist of all check if the set of coodinate is valid (tolerance on boundaries)
 if( fabs(z) > 300 || fabs(y) > 50){
   //cout << "z " << z << " and y " << y <<" are unknown coordinates" << endl;
   return lem;
 }

  int module = (int)z/100;
  int pos = (int)z % 100; //I base the classification on the first module and then rescale

  if(pos < 50 && y >= 0 )       { lem = (module*lem_in_module) + 1; }
  else if( pos < 50 && y < 0 )  { lem = (module*lem_in_module) + 3; }
  else if( pos >= 50 && y >= 0 ){ lem = (module*lem_in_module) + 2; }
  else if( pos >= 50 && y < 0 ) { lem = (module*lem_in_module) + 4; }
  else {
    cout << "z " << pos << " y " << y <<" are unknown coordinate: LEM not found " << endl;
    return lem;
  }

  return lem;
}

bool isGood_lem( int lem ){
  //select only good lems among a list

  vector<int> lems = { 2, 4, 5, 6, 7, 8, 9, 11}; //active lems

   if (find(lems.begin(),lems.end(), lem) != lems.end())
     return true;
   else
      return false;
}

bool isGood_lem(vector<int> lems, int lem){
  //select only good lems among a list
  if (find(lems.begin(),lems.end(), lem) != lems.end())
     return true;
  else
      return false;
}

double find_projection(hit h){
  //return the correct variable projection on the view
  double coor = -99999;
  switch (h.view) {
    case 0:
      coor = h.sp_y;
      break;
    case 1:
      coor = h.sp_z;
      break;
    default:
      cout << " Unknown projection " << endl;
  }
  return coor;
}
*/
