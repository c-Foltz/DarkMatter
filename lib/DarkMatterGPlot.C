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

Double_t Pol2(Double_t *x, Double_t *par)
//Polynomial 
{ 
  Double_t arg1= 0; 
  arg1=par[0]+ par[1]*x[0]+ par[2]*x[0]*x[0] ; 
  Double_t fitval = arg1 ; 
  /* 
  cout <<par[0]<<" "<<par[1]<<" "<<par[2]<<"\n"; 
  cout <<par[3]<<" "<<par[4]<<" "<<par[5]<<"\n"; 
  cout <<x[0]<<" "<<arg1<<" "<<arg2<<" "<<arg3<<" "<<fitval<<" "<<"\n"; 
  */ 
  return fitval; 
}

Double_t Gauss(Double_t *x, Double_t *par)  
//Gauss fitting 
{
    Double_t xnew=(x[0]-par[1]);
    Double_t fitval=par[0]*exp(-xnew*xnew/2.0/par[2]/par[2]);
    return fitval;
}
Double_t DoubleGauss(Double_t *x, Double_t *par)  
//Two-Gaussian fitting with an extra constant term.
{
    Double_t xnew1=(x[0]-par[1]);
    Double_t xnew2=(x[0]-par[4]);
    Double_t fitval=par[0]*exp(-xnew1*xnew1/2.0/par[2]/par[2])+par[3]*exp(-xnew2*xnew2/2.0/par[5]/par[5])+par[6];
    return fitval;
}
Double_t GaussPoly(Double_t *x, Double_t *par)  
//Gauss+Pol3 fitting 
{
    Double_t xx=x[0];
    Double_t xmid=(xx-par[1]);
    Double_t fitval=par[0]*exp(-xmid*xmid/2.0/par[2]/par[2]);
    fitval=fitval+par[3]+par[4]*xmid+par[5]*xmid*xmid+par[6]*xmid*xmid*xmid;
    return fitval;
}

Double_t GaussPoly2(Double_t *x, Double_t *par)  
//Gauss+Pol2 fitting 
{
    Double_t xnew=(x[0]-par[1]);
    Double_t fitval=par[0]*exp(-xnew*xnew/2.0/par[2]/par[2]);
    fitval=fitval+par[3]+par[4]*xnew+par[5]*xnew*xnew;
    return fitval;
}

void printToFile(double input_array[],const char* name){
  const double size = sizeof(input_array);
  // double* ptr = input_array;

  cout<<size<<"Length"<<endl;

  ofstream myfile (name);
  if (myfile.is_open())
  {
    // myfile << "This is a line.\n";
    // myfile << "This is another line.\n";
    for(int count = 0; count < 100; count ++){
        myfile << input_array[count] <<endl ;
        // ptr++;
    }
    myfile.close();
  }
  else cout << "Unable to open file";
  return 0;
}
////////////////////////////////////////////////////////////////

void DarkMatterGPlot(){
  // read in data
  vector<Double_t> vx;
  vector<Double_t> vy;
  Double_t xdat,ydat;
  double xx[1000],yV[1000],yB[1000],yP[1000],yD[1000],yT[1000],xv[1000],xmv[1000];
  double xxerror[1000],yVerror[1000],yDerror[1000],yPerror[1000],yTerror[1000];

  // Read in Data

  fstream infile1;
  char filename[50]={"GS4_Trace72"}; //************************

  cout<<" "<< endl;
  cout<<"File Name to Process = "<<filename<< endl;
  infile1.open(filename, ios_base::in);
  char c1dat[20],c2dat[20],c3dat[20],c4dat[20];
  infile1>>c1dat>>c2dat>>c3dat>>c4dat;  //skip first line;
  cout <<c1dat<<" "<<c3dat<< endl;

  int ncounter=0;
  int Type=1; //0=voltage, 1=dBm
  while (infile1>>xdat>>ydat){
    //cout <<ncounter<<" "<<xdat<<" "<<ydat<< endl;
    ncounter=ncounter+1;
    xdat=xdat/1.0e+6; //MHz
    if(fabs(ydat)<1.0) Type=0; //0: voltage, 1: dBm    
    if(fabs(ydat)<1.0) ydat=ydat*1e+6;  //change to micro-volt
    vx.push_back(xdat);
    vy.push_back(ydat) ;
  }
  //infile1.close();
  Int_t vsize = vy.size();
  cout<<"data size ="<<vsize<<" "<<vx[0]<<" "<<vx[vsize-1]<< endl;  
  cout<<"data type ="<<Type<<" 0:V, 1:dBm "<< endl;

  //This is to check the data is what is expected
  double binw=0.0000033333333;
  int ntbin=(vx[vsize-1]-vx[0])/binw;
  ntbin=ntbin+1.;
  int nbin=ntbin;
  cout<<"bin size ="<<ntbin<<" "<<nbin<< endl;

  // book histogram
   for (Int_t i=0; i<vsize; i++) {
    xx[i]=vx[i];
    xv[i] = 3e5*(1.0-1420.40/xx[i]); // conver to speed in km
    xmv[i]= xv[i];
    if(Type==0) {
         yV[i]=vy[i];
	 double tmp=vy[i]/1000000.0;  //to volt
	 yD[i]=10.*log10(tmp*tmp*1000./50); //conver to dBm. 1000 for mili, 50 for impedance
         yP[i]=1e12*tmp*tmp/50.0;  //convert to power/flu
                 }
    if(Type==1) {
         yD[i]=vy[i];
	 yV[i]=1000000.*sqrt( TMath::Power(10,vy[i]/10.0)*50.0/1000.0);
	 yP[i]=1e9*TMath::Power(10,vy[i]/10.0);
                 }
      xxerror[i]=0.0;
      yVerror[i]=0.03;  //set the errors
      yDerror[i]=0.03;      
      yPerror[i]=0.03;
      yTerror[i]=2.0;
   }
    
    TGraphErrors *gtdata = new TGraphErrors(ncounter,xx,yV,xxerror,yVerror);
    TGraphErrors *gwdata = new TGraphErrors(ncounter,xx,yD,xxerror,yDerror);
    TGraphErrors *gxdata = new TGraphErrors(ncounter,xx,yP,xxerror,yPerror);

    TGraphErrors *gfdata = new TGraphErrors(ncounter,xv,yV,xxerror,yVerror);
    TGraphErrors *ghdata = new TGraphErrors(ncounter,xv,yD,xxerror,yDerror);
    TGraphErrors *gkdata = new TGraphErrors(ncounter,xv,yP,xxerror,yPerror);   

  gtdata->GetXaxis()->SetTitle(" MHz      ");
  gtdata->GetYaxis()->SetTitle(" micro-V   ");  
  gtdata ->SetMarkerStyle(20);
  gtdata ->SetMarkerSize(0.5);
  TCanvas *myc1 =new TCanvas("myc1","Sig (V)");
  gtdata->Draw("ALP");

  gwdata->GetXaxis()->SetTitle(" MHz      ");
  gwdata->GetYaxis()->SetTitle(" dBm         ");  
  gwdata ->SetMarkerStyle(20);
  gwdata ->SetMarkerSize(0.5);
  TCanvas *myc2 =new TCanvas("myc2","Sig (dBm)");
  gwdata->Draw("ALP");

  gxdata->GetXaxis()->SetTitle(" MHz      ");
  gxdata->GetYaxis()->SetTitle(" Power ( x10E-12)        ");    
  gxdata ->SetMarkerStyle(20);
  gxdata ->SetMarkerSize(0.5);
  TCanvas *myc3 =new TCanvas("myc3","Sig (Power)");
  gxdata->Draw("ALP");  

  gkdata->GetXaxis()->SetTitle(" Km/sec      ");
  gkdata->GetYaxis()->SetTitle(" Power ( x10E-12)        ");  
  gkdata ->SetMarkerStyle(20);
  gkdata ->SetMarkerSize(0.5);
  TCanvas *myc6 =new TCanvas("myc6","Sig (Power)");
  gkdata->Draw("ALP");

  double ymax=5.25;
  double mult=1.10;
  //First delete sharp peaks
   for (Int_t i=3; i<vsize-3; i++) {
    double ytmp=yP[i];
    //if(ytmp < ymax) continue;
    double ytmpM=yP[i-2];
    double ytmpP=yP[i+2];
    double yavg=(ytmpM+ytmpP)/2.0;
    if(ytmp < mult*yavg) continue;
    cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    yP[i]=yavg;
    yP[i-1]=yavg;
    yP[i+1]=yavg;}

  
  //delete the signal region and fit the spectrum with second degree polynomial.
  double x1del=-50.0;
  double x2del= 50.0;
  int npoint=0;
  double XB[1000],YB[1000];
   for (Int_t i=0; i<vsize; i++) {
    double xtmp=xv[i];
    if(xtmp > x1del && xtmp < x2del) continue;
    XB[npoint]=xtmp;
    YB[npoint]=yP[i];
    npoint=npoint+1;}
   
    TGraphErrors *gndata = new TGraphErrors(npoint,XB,YB,xxerror,yPerror);
   gndata->GetXaxis()->SetTitle(" Km/sec      ");
   gndata->GetYaxis()->SetTitle(" Power ( x10E-12)        ");  
   gndata ->SetMarkerStyle(20);
   gndata ->SetMarkerSize(0.5);
   TCanvas *mye1 =new TCanvas("mye1","Bag (Power)");
   gndata->Draw("ALP");

  //Fitting background
  Double_t par[6];  //
  //define POL2 fit
  TF1 *mfit1 = new TF1("mfit1",Pol2, -300., 300.,3); // range & number of parameters
  mfit1 ->SetParameters(4.5,0.0,0.0);  
  gndata ->Fit("mfit1","R"," ",-230.0,230.0);

   //subtrack the fitted background
   double p0 = mfit1 ->GetParameter(0);
   double p1 = mfit1 ->GetParameter(1);
   double p2 = mfit1 ->GetParameter(2);   
   npoint=0;
   double P2T =50.0;  //comvert power to T
   double XV[1000],YP[1000],YT[1000];
   double xtmp;
   for (Int_t i=0; i<vsize; i++) {
    xtmp=xv[i];
    XV[i]=xtmp;
    double ybag=p0+p1*xtmp+p2*xtmp*xtmp;
    //double ybag=4.2;  //constant background
    YP[i]=yP[i]-ybag;
    YT[i]=P2T*YP[i];  //convert power to Temperature
    npoint=npoint+1;}
   
    TGraphErrors *gmdata = new TGraphErrors(npoint,XV,YT,xxerror,yTerror);
   gmdata->SetTitle("Sig-Bag (T)");
   gmdata->GetXaxis()->SetTitle(" Km/sec      ");
   gmdata->GetYaxis()->SetTitle(" T (Kelvin)   ");  
   gmdata ->SetMarkerStyle(20);
   gmdata ->SetMarkerSize(0.5);
   TCanvas *mye2 =new TCanvas("mye2",filename);
   gmdata->Draw("ALP");
   TF1 *fa1 = new TF1("fa1","0.0",-230.,230.);
   fa1->Draw("same");  
  
   //Rebin to 100
   double rXV[1000],rYT[1000];  //rebin
   int NC=0;
   for (Int_t i=0; i<vsize-6; i=i+6) {
     double xsum=0.0;
     double ysum=0.0;
     double counter=0;
     for(int j=i; j<i+6; j++){
       if(YT[j]==0) continue;
       counter=counter+1;
       xsum=xsum+XV[j];
       ysum=ysum+YT[j];
     }
     rXV[NC]=xsum/counter;
     rYT[NC]=ysum/counter;
     NC=NC+1;
   }

    TGraphErrors *gpdata = new TGraphErrors(NC,rXV,rYT,xxerror,yTerror);
    gpdata->SetTitle("Sig-Bag (T) vel");
   gpdata->GetXaxis()->SetTitle(" Km/sec      ");
   gpdata->GetYaxis()->SetTitle("  T (Kelvin)      ");  
   gpdata ->SetMarkerStyle(20);
   gpdata ->SetMarkerSize(0.5);
   TCanvas *mye3 =new TCanvas("mye3",filename);
   gpdata->Draw("ALP");
   fa1->Draw("same");  

  int len_NC = sizeof(NC);
  int len_rXV = sizeof(rXV);
  int len_rYT = sizeof(rYT);
  int len_xxerror=sizeof(xxerror); 
  int len_yTerror=sizeof(yTerror);

  cout<<len_NC<<endl;
  cout<<len_rXV<<"realLength"<<endl;
  cout<<len_rYT<<endl;
  cout<<len_xxerror<<endl;
  cout<<len_yTerror<<endl;

  cout<<"File Name to Processed = "<<filename<< endl;
  cout<<" "<< endl;
  
  printToFile(rXV,"rXV.txt");
  printToFile(rYT,"rYT.txt");
  printToFile(xxerror,"xxerror.txt");
  printToFile(yTerror,"yTerror.txt");

  // for(int i=0;i<1000;i++){
  //   cout<<yTerror[i]<<endl;
  // }


   }
   
