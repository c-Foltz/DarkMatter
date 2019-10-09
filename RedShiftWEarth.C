
#include <iostream>
#include <fstream>
#include <vector>

double aatan(double py, double px) {
  double R=6.2831853;
  double pxx=px;
  double pyy=py;
  if(px == 0.0) pxx=0.000001;
  if(py == 0.0) pyy=0.000001;
  double tmp = atan(fabs(pyy/pxx));
  if(pxx >= 0.0  && pyy >= 0.0) return tmp;
  if(pxx >= 0.0  && pyy <  0.0) return R-tmp;
  if(pxx <  0.0  && pyy >= 0.0) return 0.5*R-tmp;
  if(pxx <  0.0  && pyy <  0.0) return 0.5*R+tmp;
  return 0;
}

void RedShiftWEarth(){

  gStyle->SetOptStat(0);
  double rad=180.0/3.141592;
  double theta_r=180.0;
  int itheta=2.0*theta_r+0.5;
  double theta_sun=180.0/rad;
  double R_MW=55.0;  //diameter milkyway : 100 K light years
  double Rinner=2.0;

  double R_sun=25.0; //from the center of galaxy
  double X_sun=R_sun*cos(theta_sun);
  double Y_sun=R_sun*sin(theta_sun);
  double Z_sun=0.0;

  double VR=220.;   //km/sec plateau rotational speed
  double VR_sun=VR;
  double VX_sun=-VR*sin(theta_sun);
  double VY_sun= VR*cos(theta_sun);
  double VZ_sun= 0.0;
  int iwrite=0;

  // book histogram
  TH1F *hist1 = new TH1F("hist1","Relative Speed",    100,-VR,VR);
  TH1F *hist2 = new TH1F("hist2","Relative Max Speed",100,-VR,VR);
  TH2F *hist11 = new TH2F("hist11","VR rr .vs. angle", itheta,-theta_r,theta_r, 200, 0.0, 200.0);
  TH2F *hist12 = new TH2F("hist12","VR_max .vs. angle",itheta,-theta_r,theta_r, 100,-3*R_MW,3*R_MW);
  TH2F *hist13 = new TH2F("hist13","VR_max position rr .vs. angle",itheta,-theta_r,theta_r, 100,-3*R_MW,3*R_MW);
  TH1F *hist5 = new TH1F("hist5","Correction to the Vrel due to earth orbiting",360,-theta_r,theta_r);

  char *histname = new char[100];
  char *histtitle = new char[100];
  TH1F *Vrel[50];
  for(int ii=0;ii<36;ii++) {
    double theta1=-180+ii*10.;
    double theta2=theta1+10.;
      sprintf(histname, "Vrel_%dA",ii);
      sprintf(histtitle,"Vrel : function of rr for theta from %g to %g",theta1,theta2);
      Vrel[ii] =  new TH1F(histname,histtitle,100,0.0,100.);}

  int IDate;
  cout<<" "<<endl;
  cout<<"Enter the nuber of days counted from Jan 1 : "<< endl;
  cin>>IDate;
  //int IDate=240;
  int IShift=14;
      IDate=IDate-IShift;  //earth orbit wrt the sun in days. 0 is February 1  ********!!!!!!!!!!!!!!!!!!!********
  double Date=IDate/365.*2.0*3.141592; //convert to angle
  double R_earth=1.58e-8; // earth distance (kly) from the sun;
  double V_earth=30.0;    //earth around the sun in km/sec
  double alpha=28.0/rad;  //earth orbit tilt in GC

  double X_earth=X_sun + R_earth*cos(Date);
  double Y_earth=Y_sun + sin(alpha)*R_earth*sin(Date);
  double Z_earth=Z_sun + cos(alpha)*R_earth*sin(Date);

  double VX_earth=VX_sun - V_earth*sin(Date);
  double VY_earth=VY_sun + sin(alpha)*V_earth*cos(Date);
  double VZ_earth=VZ_sun + cos(alpha)*V_earth*cos(Date);

  double dtheta=2.0;  
  for (double theta=-theta_r; theta <theta_r; theta=theta+dtheta) { 
    //  { double theta =15.0;
    double Vmax=0.0;
    double rmax=0.0;
  for (double rr=1.5; rr <2.0*R_MW; rr =rr+1.0) {
    double Xwrt_sun= rr*cos(theta/rad); //now the star
    double Ywrt_sun= rr*sin(theta/rad);
    double Zwrt_sun=0.0;
    double Rwrt_sun=sqrt(Xwrt_sun*Xwrt_sun + Ywrt_sun*Ywrt_sun);
    
    double Xwrt_gc=Xwrt_sun - R_sun;  //the star wrt gc
    double Ywrt_gc=Ywrt_sun;
    double Zwrt_gc=0.0;
    double Rwrt_gc=sqrt(Xwrt_gc*Xwrt_gc + Ywrt_gc*Ywrt_gc);

    if(iwrite == 1){
    cout <<" "<< endl;
    cout <<"loop :  "<<theta<<" "<<rr<<endl;
    cout <<"rr_sun :  "<<Xwrt_sun<<" "<<Ywrt_sun<< endl;
    cout <<"rr_gc :  "<<Xwrt_gc<<" "<<Ywrt_gc<< endl;}
    if(Rwrt_gc > R_MW) continue;
    if(Rwrt_gc < Rinner) continue;

    double VRwrt_gc=VR; //star rotational speed. 220 km/sec
    //double VRwrt_gc=14.0*Rwrt_gc/3.0+110.0;       //a test
    //if(Rwrt_gc < 10.) VRwrt_gc=VR*Rwrt_gc/10.0;   //a test
    if(Rwrt_gc < 10.) VRwrt_gc=VR*sqrt(Rwrt_gc)/sqrt(10.0);   
    double theta_gc=aatan(Ywrt_gc,Xwrt_gc);
    double VXwrt_gc=-VRwrt_gc*sin(theta_gc);    //star speed wrt gc
    double VYwrt_gc= VRwrt_gc*cos(theta_gc);

    double RR=sqrt(R_sun*R_sun+rr*rr-2.0*R_sun*rr*cos(theta/rad)); //distance to rr from GC
    //double Vdot=VRwrt_gc*R_sun/RR*sin(theta/rad)-VR*sin(theta/rad);
    double Vdot1=(VX_sun*Xwrt_sun   +  VY_sun*Ywrt_sun)/rr;
    double Vdot2=(VXwrt_gc*Xwrt_sun +  VYwrt_gc*Ywrt_sun)/rr;
    double Vdot=Vdot1-Vdot2;

    double VEdot1=(VX_earth*Xwrt_sun +VY_earth*Ywrt_sun + VZ_earth*Zwrt_sun)/rr;
    double VEdot2=(VXwrt_gc*Xwrt_sun +  VYwrt_gc*Ywrt_sun)/rr;
    double VEdot=VEdot1-VEdot2;
    if(rr == 1.5) {hist5 ->Fill(theta,-Vdot+VEdot);
      int binN=theta+181.5;
      hist5 ->SetBinError(binN,0.0);} 

    if(iwrite == 1) {
    cout <<"v_sun :  "<<VX_sun<<" "<<VY_sun<<"    Vwrt_gc:  "<<VXwrt_gc<<" "<<VYwrt_gc<< endl;
    cout <<"vdot1,  vdot2  ="<<Vdot1<<" "<<Vdot2<< endl;
    cout <<"vEdot1, vEdot2 ="<<VEdot1<<" "<<VEdot2<< endl;
    cout <<"vdv => "<<Vdot<< endl;
    cout <<"ved => "<<VEdot<< endl;
    cout <<"Correction to the plot (km/s) from earth orbit V => "<<-Vdot+VEdot<< endl;}
    
    hist1 ->Fill(-Vdot);
    hist11 ->Fill(theta,rr,fabs(Vdot));
    int ihist=(180.0+theta)/10.;
    Vrel[ihist] ->Fill(rr,-Vdot/10.*dtheta);  //10 degree plot step, 2 degree step
    if(fabs(Vdot) > fabs(Vmax)) {
            Vmax=Vdot;
	    rmax=rr;
                                }

  } //end of theta loop
  hist2 ->Fill(Vmax);
  hist12 ->Fill(theta,Vmax,1.0);
  hist13 ->Fill(theta,rmax,Vmax);
  }// end of rr loop

  hist1->GetXaxis()->SetTitle("Vrel (Km/sec)       ");
  hist1->GetYaxis()->SetTitle("events              ");  

  hist5->GetXaxis()->SetTitle("Angle wrt GC (degree)     ");
  hist5->GetYaxis()->SetTitle("V correction (km)            ");
  hist5 ->SetMarkerStyle(20);  

  hist11->GetXaxis()->SetTitle("angle (degree)       ");
  hist11->GetYaxis()->SetTitle("rr (x1000 ly)         ");  

  hist12->GetXaxis()->SetTitle("angle (degree)       ");
  hist12->GetYaxis()->SetTitle(" Vrel Max  (Km/sec)     "); 
  hist12->SetMarkerStyle(20); 

  TFile *f = new TFile("VRELplots.root","recreate");
  for(int j=0 ; j<36 ; j++){
   Vrel[j] ->GetXaxis()->SetTitle(" rr (x 1000 light years)           ");
   Vrel[j] ->GetYaxis()->SetTitle("Vrel  (Km/sec)                ");  
   Vrel[j] ->Write();
  }
  hist5 ->Draw();
  hist1 ->Write();
  hist5 ->Write();
  hist11 ->Write();
  hist12 ->Write();
  hist13 ->Write();
  f ->Close();

  cout<<" "<< endl;
  cout<<"Run finished for the earth Day = "<< IDate+IShift <<" days from Jan 1 "<<endl;
  cout<<" "<< endl;
  cout<<" "<< endl;
}


Double_t Exp1(Double_t *x, Double_t *par)
{
  Double_t arg1= 0;
  arg1 = par[0]*exp(par[1]*x[0]);
  Double_t fitval = arg1  ;
  return fitval;
}

Double_t Ruby(Double_t *x, Double_t *par)
{
  double par0=par[0];  //amplitude
  double par1=par[1];  //g factor
  double par2=par[2]/57.295;  //crystal axis
  double par3=par[3]/57.295;  //zero point
  double xx=x[0]/57.295;

  double a1=par1*par1*sin(par2)*sin(par2)*cos(xx+par3)*cos(xx+par3);
  double a2=4.0*(1.0-sin(par2)*sin(par2)*cos(xx+par3)*cos(xx+par3));
  double fitval=par0/sqrt(a1+a2);
  return fitval;
}


