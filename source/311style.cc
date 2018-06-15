////////////////////////////////////////////////////////////////////////////////
//
// Set of libraries to handle styles for plot uniformity
//
////////////////////////////////////////////////////////////////////////////////


//project libraries
#include <311Lib.h>
#include <311style.h>


Color_t getColor(int view){

  Color_t color = kBlack;

  switch (view) {
    case 1:
      color = kBlue;
      break;
    case 0:
      color = kRed;
      break;
    default:
      break;
  }

  return color;
}
