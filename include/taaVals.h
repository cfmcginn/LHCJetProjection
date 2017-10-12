#ifndef TAAVALS_H
#define TAAVALS_H

double getTAA(const int centLow, const int centHi)
{
  double taaVal = 0.;

  if(centLow == 0 && centHi == 5) taaVal = 25.98;
  else if(centLow == 5 && centHi == 10) taaVal = 20.46;
  else if(centLow == 10 && centHi == 30) taaVal = 11.51;
  else if(centLow == 30 && centHi == 50) taaVal = 3.819;
  else if(centLow == 50 && centHi == 70) taaVal = 0.9345;
  else{
    std::cout << "Input centLow-centHi \'" << centLow << "-" << centHi << "\' is invalid" << std::endl;
  }

  return taaVal;
}

#endif


