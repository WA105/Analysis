////////////////////////////////////////////////////////////////////////////////
// Geometry class definition
//
//mailto:andrea.scarpelli@cern.ch
////////////////////////////////////////////////////////////////////////////////

#ifndef  __GEOMETRY_H
#define __GEOMETRY_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

const int NUM_OF_VIEWS = 2;
const int NUM_OF_LEMS = 12;
const int tpc_boundaries[6] = {-50, 50, -50, 50, 0, 300}; //minx,maxx,miny,maxy,minz,max
const int Ch_0 = 320;
const int Ch_1 = 960;
const double pitch = 0.3;
const int tdc = 1667;
const int  maxNumChannels = 1280;
const int maxNumTdc = maxNumChannels*tdc;
const float sampling_freq = 2.5;    //sampling feq in MHz
const double ADC2CHARGE = 45.31875; //ADC*ticks (from qScan)
const double adc2fc0=1./56.5; //ADC*ticks to fC view 0 (from pulsing)
const double adc2fc1=1./67.4; //ADC*ticks to fC view 1 (from pulsing)

const int lem_in_module = 4; //num of minimal modules of 4 lems
const double lem_size = 50; //in cm

unsigned int ViewToDAQChan(unsigned int ViewChan);

//int find_lem(double y, double z);

//bool isGood_lem( int lem );

//bool isGood_lem(vector<int> lems, int lem);

//double find_projection(hit h);

#endif // __GEOMETRY_H
