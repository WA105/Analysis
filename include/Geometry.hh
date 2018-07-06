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
const double Ch_0 = 320;
const double Ch_1 = 960;
const double pitch = 0.3;
const int tdc = 1667;
const double ADC2CHARGE = 45.31875; //ADC*ticks (from qScan)

//int find_lem(double y, double z);

//bool isGood_lem( int lem );

//bool isGood_lem(vector<int> lems, int lem);

//double find_projection(hit h);

#endif // __GEOMETRY_H
