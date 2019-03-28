//Calculate the integral of the calorimeter response function of DUNE based on
//the pulsing measurements

#include "TF1.h"

double getIntegral( TF1 *f, const int nticks,
                            double samplingPeriod, double mv2fc, double ADC2mV){
  double integral =0;
  double element[nticks];
  double max = 0.;

  for(size_t t=0; t<nticks; t++)
  {
    double time = double(t)*samplingPeriod; //convert from ticks to us
    element[t] = f->Eval(time);
    cout << time << " " << f->Eval(time) << endl;
    if( element[t] > max ){ max = element[t]; }
  }

  for(size_t t=0; t<nticks; t++)
  {
    element[t] /= max;
    element[t] *= mv2fc;
    element[t] *= ADC2mV;
    integral += element[t];
  }

  return integral;
}

void caloAmp()
{
  //define response function----------------------------------------------------
  //View 0 (3 m)
  TF1 *fAmp0 = new TF1("amp0", "[0]*exp((x-[1])/[2]) / ( 1 + exp((x-[1])/[3]))", 0., 667);
  fAmp0->SetParameters(1.21581, 8.64873e-1, 2.88648e-1, 2.7413e-1 );
  fAmp0->SetNpx(100000);
  //View 1 (1 m)
  TF1 *fAmp1 = new TF1("amp1", "[0]*exp((x-[1])/[2]) / ( 1 + exp((x-[1])/[3]))", 0., 667);
  fAmp1->SetParameters(1.46135, 1.0228, 3.36029e-1, 2.9685e-1 );
  fAmp0->SetNpx(100000);

  //Get integral----------------------------------------------------------------
  double int0 = getIntegral(fAmp0, 200, 0.4, 2.53, 1.0);
  double int1 = getIntegral(fAmp1, 200, 0.4, 6.38, 1.0);
  cout << "ADCxus to fC calibration" << endl;
  cout << "View 0: " << int0  << "\nView 1: " << int1 << endl;

  //draw------------------------------------------------------------------------
  fAmp0->SetLineColor(2);
  fAmp0->Draw("");

  fAmp1->SetLineColor(4);
  fAmp1->Draw("SAME");

  fAmp0->GetXaxis()->SetTitle("Time [us]");
  fAmp0->GetXaxis()->SetRangeUser(0,20);

  fAmp0->GetYaxis()->SetTitle("Voltage [mV]");
  fAmp0->GetYaxis()->SetRangeUser(0,1.2);

  gPad->BuildLegend();
}
