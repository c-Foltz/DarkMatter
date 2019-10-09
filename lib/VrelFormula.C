

#include <iostream>
#include <fstream>
#include <string>
#include "stdio.h"
#include "stdlib.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "math.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"

void VrelFormula(){
  
  double rad=180.0/3.141592;
  double theta= 135./rad;
  double Vsun=220.0;
  double Vstar=220.0;
  double Rsun=8.5;  //unit is parsec

  int ic=0;
  double xx[1000],yy[1000];
  double xxerror[1000],yyerror[1000];

  for(double rr=0.5; rr<30.0; rr=rr+0.25){
    double Sint=sin(theta);
    double Cost=cos(theta);
    double RR=sqrt(rr*rr + Rsun*Rsun -2.0*rr*Rsun*Cost);
    double isign=-1.0;
    if((rr*rr) > (RR*RR+Rsun*Rsun)) isign= 1.0;
    double Vrel=0.0;
    Vrel = Vsun*Sint - Vstar*rr/RR*Sint*Cost + isign*Vstar*sqrt(1.0-(rr*rr/RR/RR*Sint*Sint))*Sint;
    xx[ic]=rr;
    yy[ic]=Vrel;
    xxerror[ic]=0.01;
    yyerror[ic]=0.01;
    ic=ic+1;
    //cout<<rr<<" "<<Vrel<< endl;
  }
    TGraphErrors *gtdata = new TGraphErrors(ic-1,xx,yy,xxerror,yyerror);
    gtdata ->Draw("ALP");
}
