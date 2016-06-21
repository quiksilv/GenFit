/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "HelixTrackModel.h"
#include <FieldManager.h>

#include <MaterialEffects.h>
#include <RKTrackRep.h>
#include <TGeoMaterialInterface.h>

#include <MeasuredStateOnPlane.h>
#include <MeasurementOnPlane.h>


#include <assert.h>
#include <math.h>
#include <TMath.h>
namespace genfit {

HelixTrackModel::HelixTrackModel(const TVector3& pos, const TVector3& mom, double charge) {

  mom_ = mom.Mag();

  TVector3 B = genfit::FieldManager::getInstance()->getFieldVal(pos);
  B.Print();

  // B must point in Z direction
//  assert(B.X() == 0);
//  assert(B.Y() == 0);

//  double Bz = B.Z();
  double Bz = 10000;

  // calc helix parameters
  TVector3 dir2D(mom);
  dir2D.SetZ(0);
  dir2D.SetMag(1.);
  R_ = 100.*mom.Perp()/(0.0299792458*Bz) / fabs(charge);
  sgn_ = 1;
  if (charge<0) sgn_=-1.;
  center_ = pos + sgn_ * R_ * dir2D.Orthogonal();
  alpha0_ = (pos-center_).Phi();

  theta_ = mom.Theta();

  //std::cout<<"radius " << R_ << "  center ";
  //center_.Print();

}


TVector3 HelixTrackModel::getPos(double tracklength) const {

  TVector3 pos;

  double angle = alpha0_ - sgn_ * tracklength / R_ * sin(theta_);

  TVector3 radius(R_,0,0);
  radius.SetPhi(angle);
  pos = center_ + radius;
  pos.SetZ(center_.Z() - sgn_ * ((alpha0_-angle)*R_ * tan(theta_-M_PI/2.)) );

  return pos;
}



void HelixTrackModel::getPosMom(double tracklength, TVector3& pos, TVector3& mom) const {


  // wilfrid
  // trackrep
  const int pdg = 11 ;     
  genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);
  //rep-> setDebugLvl(10);


  double tracklengthNomateff =0.;
  TVector3 posNomateff ;
  TVector3 thiscenter = this->center_;
  TVector3 thismom ;

  this->getPosMomNomateff(tracklengthNomateff, posNomateff, thismom);

  const double resolution = 0.02;
  const unsigned int nMeasurements =60;

  TMatrixDSym covM(6);
  for (int i = 0; i < 3; ++i)
    covM(i,i) = resolution*resolution;
  for (int i = 3; i < 6; ++i)
    covM(i,i) = pow(resolution / nMeasurements / sqrt(3), 2);



// smeared start state
      genfit::MeasuredStateOnPlane stateSmeared(rep);
      //  rep->setPosMomCov(stateSmeared, posNomateff, thismom, covM);
  rep->setPosMom(stateSmeared, posNomateff, thismom);

  TVector3 momstateSmeared =  rep->getMom(stateSmeared);


  //std::cout << "HelixTrackModel::getPosMom mom init  = " <<  thismom.X() << " Y :  "<<thismom.Y() <<" Z :  "<<  thismom.Z() << std::endl;

  // std::cout << "HelixTrackModel::getPosMom pos  = " << " end mom "<< std::endl;

  double angle = alpha0_ - sgn_ * tracklength / R_ * sin(theta_);

  TVector3 radius(R_,0,0);
  radius.SetPhi(angle);
  pos = center_ + radius;
  pos.SetZ(center_.Z() - sgn_ * ((alpha0_-angle)*R_ * tan(theta_-M_PI/2.)) );

  mom.SetXYZ(1,1,1);
  mom.SetTheta(theta_);
  mom.SetPhi(angle - sgn_*M_PI/2.);
  mom.SetMag(mom_);

  /*std::cout<<"tracklength " << tracklength << "\n";
  std::cout<<"angle " << angle << "\n";
  std::cout<<"radius vector "; radius.Print();
  std::cout<<"pos "; pos.Print();
  std::cout<<"mom "; mom.Print();*/




  double extrapLenTest(0);
 
  extrapLenTest = rep->extrapolateToPoint(stateSmeared, pos);
 
  TVector3 posMateff =  rep->getPos(stateSmeared);
  TVector3 momMateff =  rep->getMom(stateSmeared);


  //std::cout << "HelixTrackModel::getPosMom pos no mateff = " <<  pos.X() << " Y :  "<<pos.Y() <<" Z :  "<<  pos.Z() << std::endl;

  //  std::cout << "HelixTrackModel::getPosMom pos mateff = " <<  posMateff.X() << " Y :  "<<posMateff.Y() <<" Z :  "<<  posMateff.Z() << std::endl;

  //  std::cout << "HelixTrackModel::getPosMom pos  = " << " end pos "<< std::endl;

  //std::cout << "HelixTrackModel::getPosMom mom no mateff = " <<  mom.X() << " Y :  "<<mom.Y() <<" Z :  "<<  mom.Z() << std::endl;

  //std::cout << "HelixTrackModel::getPosMom mom mateff  = " <<  momMateff.X() << " Y :  "<<momMateff.Y() <<" Z :  "<<  momMateff.Z() << std::endl;



  // std::cout << "HelixTrackModel::getPosMom pos  = " << " end mom "<< std::endl;


      pos= posMateff;
      mom=momMateff;

  delete rep;
 

}

void HelixTrackModel::getPosMomNomateff(double tracklength, TVector3& pos, TVector3& mom) const {


  // wilfrid
    // trackrep
  //  const int pdg = 11 ;     
  //  genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);
    //rep-> setDebugLvl(10);

// smeared start state
  // genfit::MeasuredStateOnPlane stateSmeared(rep);
  // rep->setPosMomCov(stateSmeared, posM, momM, covM);


  double angle = alpha0_ - sgn_ * tracklength / R_ * sin(theta_);

  TVector3 radius(R_,0,0);
  radius.SetPhi(angle);
  pos = center_ + radius;
  pos.SetZ(center_.Z() - sgn_ * ((alpha0_-angle)*R_ * tan(theta_-M_PI/2.)) );

  mom.SetXYZ(1,1,1);
  mom.SetTheta(theta_);
  mom.SetPhi(angle - sgn_*M_PI/2.);
  mom.SetMag(mom_);

  /*std::cout<<"tracklength " << tracklength << "\n";
  std::cout<<"angle " << angle << "\n";
  std::cout<<"radius vector "; radius.Print();
  std::cout<<"pos "; pos.Print();
  std::cout<<"mom "; mom.Print();*/

}



std::vector<double> HelixTrackModel::GetGenTrack(){
  //std::vector<double> HelixTrackModel::GetGenTrack( TVector3 pos , TVector3 center, TVector3 mom ){
 

  double tracklength =0.;
  TVector3 pos=this->getPos(tracklength);
  TVector3 center = this->center_;


  std::vector<double> vecZ1;
  std::vector<double> vecZ2;

  std::vector<double> vecPhihelp;
  std::vector<double> vecPhihelm;

 std::vector<double> vecPhihel;
std::vector<double> vecPhi;
 std::vector<double> tracklength_hel;

  std::vector<TVector3> vecPoslp;
  std::vector<TVector3> vecPoslm;

  std::vector<TVector3> vecPos;


  double twopi = 8.*atan(1.);

  //  double rho=60.; // 3.;

 
  double ZhelMin = -2000.;
  double ZhelMax = 2000.;


  //   double A=0.1;
  //  double B=0.4;
  //  double C=0.3;
  //  double tlam=0.3;

  //    double Xc = 60.; //0.5; //0.5;
  //  double Yc = 0.2; //-0.6;


  double A= pos.X();
    double B=pos.Y();
    double C= pos.Z();

    //    double tlam= - TMath::Tan(mom.Theta()-TMath::Pi()*0.5);
   double tlam= - TMath::Tan(theta_ -TMath::Pi()*0.5);




    double Xc = center.X();
    double Yc = center.Y();

    std::cout<< " radius =" << TMath::Sqrt((A-Xc)* (A-Xc)+(B-Yc)*(B-Yc))  << " center X Y = " << Xc << " " <<  Yc<< " pos X Y Z  = " << A << " " << B <<  " " << C<<std::endl;


    //     double rmax=   1000.;// 79.0+1.6;
    //  double rmin= 41.;// 0.5;// 51.8;
        double rmax=   79.0+1.6;
	double rmin=   42.0; // 51.8; // 42.0; //42.0; // 45.0; // 51.8;
      //  double rmax=   84.;
      // double rmin= 54.;



  double CellSize =0.5; 
  int nbwires=(int)((rmax-rmin)/CellSize);

  std::cout<<"rmax = "<<rmax<<" rmin= "<<rmin<<" CellSize = "<<CellSize<<" nbwires =" <<nbwires<<std::endl;


  int nhits=0;

  for(;;){
    nhits= nhits+1;
    if(nhits>nbwires)break;

    double rd=rmin+CellSize*nhits;
    // rd = 2.;

    double delta =  -(pow (B*Xc - A*Yc, 2)*(pow (A, 4) + pow (B, 4) - 4*pow (A, 3)*Xc - 4*pow (B, 3)*Yc + 4*B*pow (rd, 2)*Yc + 
           2*pow (A, 
             2)*(pow (B, 2) - pow (rd, 2) + 2*pow (Xc, 2) - 
              2*B*Yc) + 
           4*A*Xc*(-pow (B, 2) + pow (rd, 2) + 2*B*Yc) - 
           2*pow (B, 2)*(pow (rd, 2) - 2*pow (Yc, 2)) + 
           pow (rd, 
             2)*(pow (rd, 2) - 
		 4*(pow (Xc, 2) + pow (Yc, 2))) ) ) ; 

    std::cout<< " rd = "<< rd << " delta = "<<delta<< std::endl;

    if(delta<0.) break;


    // (2.*(pow(A - Xc,2) + pow(B - Yc,2))*(pow(Xc,2) + pow(Yc,2)));



 double cosw11 =(-(((A - Xc)*Xc + (B - Yc)*Yc)*
        (pow(A,2) + pow(B,2) - pow(rd,2) - 2*A*Xc - 2*B*Yc + 
          2*(pow(Xc,2) + pow(Yc,2)))) - 
     TMath::Sqrt(-(pow(B*Xc - A*Yc,2)*
         (pow(A,4) + pow(B,4) - 4*pow(A,3)*Xc - 4*pow(B,3)*Yc + 
           4*B*pow(rd,2)*Yc + 
           2*pow(A,2)*(pow(B,2) - pow(rd,2) + 2*pow(Xc,2) - 
              2*B*Yc) + 4*A*Xc*(-pow(B,2) + pow(rd,2) + 2*B*Yc) - 
           2*pow(B,2)*(pow(rd,2) - 2*pow(Yc,2)) + 
           pow(rd,2)*(pow(rd,2) - 4*(pow(Xc,2) + pow(Yc,2)))))))/
   (2.*(pow(A - Xc,2) + pow(B - Yc,2))*(pow(Xc,2) + pow(Yc,2)));



double delta2sin = 4*pow(B*Xc - A*Yc,2) - 
        4*(pow(A - Xc,2) + pow(B - Yc,2))*
         (pow(A,2)*pow(cosw11,2) + pow(B,2)*pow(cosw11,2) - 
           pow(rd,2) - 2*A*(-1 + cosw11)*cosw11*Xc - 
           2*B*(-1 + cosw11)*cosw11*Yc + 
	  pow(-1 + cosw11,2)*(pow(Xc,2) + pow(Yc,2)));

 if(delta2sin<0.) break;

double sinw11 = (B*Xc - A*Yc - TMath::Sqrt(4*pow(B*Xc - A*Yc,2) - 
        4*(pow(A - Xc,2) + pow(B - Yc,2))*
         (pow(A,2)*pow(cosw11,2) + pow(B,2)*pow(cosw11,2) - 
           pow(rd,2) - 2*A*(-1 + cosw11)*cosw11*Xc - 
           2*B*(-1 + cosw11)*cosw11*Yc + 
           pow(-1 + cosw11,2)*(pow(Xc,2) + pow(Yc,2))))/2.)/
   (pow(A - Xc,2) + pow(B - Yc,2));




double sinw21 = (B*Xc - A*Yc + TMath::Sqrt(4*pow(B*Xc - A*Yc,2) - 
        4*(pow(A - Xc,2) + pow(B - Yc,2))*
         (pow(A,2)*pow(cosw11,2) + pow(B,2)*pow(cosw11,2) - 
           pow(rd,2) - 2*A*(-1 + cosw11)*cosw11*Xc - 
           2*B*(-1 + cosw11)*cosw11*Yc + 
           pow(-1 + cosw11,2)*(pow(Xc,2) + pow(Yc,2))))/2.)/
   (pow(A - Xc,2) + pow(B - Yc,2));


double cosw12 =((-((A - Xc)*Xc) - (B - Yc)*Yc)*
      (pow(A,2) + pow(B,2) - pow(rd,2) - 2*A*Xc - 2*B*Yc + 
        2*(pow(Xc,2) + pow(Yc,2))) + 
     TMath::Sqrt(-(pow(B*Xc - A*Yc,2)*
         (pow(A,4) + pow(B,4) - 4*pow(A,3)*Xc - 4*pow(B,3)*Yc + 
           4*B*pow(rd,2)*Yc + 
           2*pow(A,2)*(pow(B,2) - pow(rd,2) + 2*pow(Xc,2) - 
              2*B*Yc) + 4*A*Xc*(-pow(B,2) + pow(rd,2) + 2*B*Yc) - 
           2*pow(B,2)*(pow(rd,2) - 2*pow(Yc,2)) + 
           pow(rd,2)*(pow(rd,2) - 4*(pow(Xc,2) + pow(Yc,2)))))))/
   (2.*(pow(A - Xc,2) + pow(B - Yc,2))*(pow(Xc,2) + pow(Yc,2)));


double sinw12 = (B*Xc - A*Yc - TMath::Sqrt(4*pow(B*Xc - A*Yc,2) - 
        4*(pow(A - Xc,2) + pow(B - Yc,2))*
         (pow(A,2)*pow(cosw12,2) + pow(B,2)*pow(cosw12,2) - 
           pow(rd,2) - 2*A*(-1 + cosw12)*cosw12*Xc - 
           2*B*(-1 + cosw12)*cosw12*Yc + 
           pow(-1 + cosw12,2)*(pow(Xc,2) + pow(Yc,2))))/2.)/
   (pow(A - Xc,2) + pow(B - Yc,2));

double sinw22 = (B*Xc - A*Yc + TMath::Sqrt(4*pow(B*Xc - A*Yc,2) - 
        4*(pow(A - Xc,2) + pow(B - Yc,2))*
         (pow(A,2)*pow(cosw12,2) + pow(B,2)*pow(cosw12,2) - 
           pow(rd,2) - 2*A*(-1 + cosw12)*cosw12*Xc - 
           2*B*(-1 + cosw12)*cosw12*Yc + 
           pow(-1 + cosw12,2)*(pow(Xc,2) + pow(Yc,2))))/2.)/
   (pow(A - Xc,2) + pow(B - Yc,2));





 double iden11 = pow(cosw11,2)+ pow(sinw11,2);
double iden12 = pow(cosw12,2)+ pow(sinw12,2);
 double iden21 = pow(cosw11,2)+ pow(sinw21,2);
double iden22 = pow(cosw12,2)+ pow(sinw22,2);

 double eps = 0.000001;
 double cosphip;
 double sinphip;
 double cosphim;
 double sinphim;

 if(TMath::Abs(1.-iden11)<eps){cosphip=cosw11;sinphip=sinw11;}
 if(TMath::Abs(1.-iden12)<eps){cosphip=cosw12;sinphip=sinw12;}
 if(TMath::Abs(1.-iden21)<eps){cosphim=cosw11;sinphim=sinw21;}
 if(TMath::Abs(1.-iden22)<eps){cosphim=cosw12;sinphim=sinw22;}


 std::cout<< " **************cosw11 = "<<cosw11<< " delta2sin = " << delta2sin  <<  " sinw11 = " <<sinw11  << std::endl;
 std::cout<< " **************iden11 = "  <<iden11<< "  iden12 = "  <<iden12<< "  iden21 = "  <<iden21<< "  iden22 = "  <<iden22<<std::endl;

 


     double phip = atan2(sinphip,cosphip);
      if(phip<0.) phip=phip+twopi;

     double phim = atan2(sinphim,cosphim);
      if(phim<0.) phim=phim+twopi;

      vecPhihelp.push_back(phip);
      vecPhihelm.push_back(phim);
    
     double Xdp = A*cosphip - B*sinphip + Xc - cosphip*Xc + sinphip*Yc; 
     double Ydp =  B*cosphip + A*sinphip - sinphip*Xc + Yc - cosphip*Yc;
     double Zdp = C +  TMath::Sqrt((A-Xc)* (A-Xc)+(B-Yc)*(B-Yc))* phip*tlam;

     TVector3 posd(Xdp, Ydp, Zdp);
     vecPoslp.push_back(posd);

     double Xdm = A*cosphim - B*sinphim + Xc - cosphim*Xc + sinphim*Yc;  
     double Ydm = B*cosphim + A*sinphim - sinphim*Xc + Yc - cosphim*Yc;
     double Zdm = C +  TMath::Sqrt((A-Xc)* (A-Xc)+(B-Yc)*(B-Yc))* phim*tlam;



     std::cout<< " Xdp = " << Xdp<< " "<<  GetHelPar(1 , A, B, C, phip, Xc, Yc , tlam )<<  std::endl;
     std::cout<< " Ydp = " << Ydp<< " " << GetHelPar(2 , A, B, C, phip, Xc, Yc , tlam )<<  std::endl;
     std::cout<< " Zdp = " << Zdp<< " " << GetHelPar(3 , A, B, C, phip, Xc, Yc , tlam )<<  std::endl;



     TVector3 posm(Xdm, Ydm, Zdm);
     vecPoslm.push_back(posm);

     double cosphidp =  Xdp/ rd; 
     double sinphidp =Ydp/ rd;

     double cosphidm = Xdm/ rd; 
     double sinphidm =  Ydm/rd;



     std::cout<<"nhits = "<<nhits<<" rd = "<<rd << " cosphip = " <<cosphip<< " sinphip =  "<< sinphip << " cosphip*  cosphip+sinphip* sinphip = " <<cosphip*  cosphip+sinphip* sinphip  << " cosphpim = " <<cosphim<< " sinphim =  "<<sinphim<< " cosphim*  cosphim+sinphim* sinphim = " <<cosphim*  cosphim+sinphim* sinphim  <<" phip = " << phip<<" phim = " << phim <<  " Xdp = " << Xdp <<" Ydp = " <<Ydp    << " Zdp = " <<Zdp   <<   " Xdm = " << Xdm <<" Ydm = " <<Ydm << " Zdm = " <<Zdm <<   std::endl;


     // if(abs(cosphi)>1)break;
     // if(abs(sinphi)>1)break;

  ;
 }

 nhits =nhits-1;



  for (unsigned int i=0; i<vecPhihelp.size(); ++i){
    std::cout<< " i= " <<i <<" vecPhihelp[i] = " <<vecPhihelp[i]<<   std::endl;
  }

  for (unsigned int i=0; i<vecPhihelm.size(); ++i){
    std::cout<< " i= " <<i <<" vecPhihelm[i] = " <<vecPhihelm[i]<<   std::endl;
  }

  std::cout<< " test phase *****"  <<vecPhihelm[vecPhihelm.size()-1] <<  "  " <<vecPhihelp[vecPhihelp.size()-1]<<   std::endl;

  if( vecPhihelm[vecPhihelm.size()-1]< vecPhihelp[vecPhihelp.size()-1]){
  std::cout<< "*********************  phase mp"<<   std::endl;
  for (unsigned int i=0; i<vecPhihelm.size(); ++i){
  vecPhihel.push_back(vecPhihelm[i]);
  }
  for (unsigned int i=0; i<vecPhihelp.size(); ++i){
 vecPhihel.push_back(vecPhihelp[vecPhihelp.size()-1-i]);
  }
 }else{
 std::cout<< "*********************  phase pm"<<   std::endl;
  for (unsigned int i=0; i<vecPhihelp.size(); ++i){
  vecPhihel.push_back(vecPhihelp[i]);
  }
  for (unsigned int i=0; i<vecPhihelm.size(); ++i){
 vecPhihel.push_back(vecPhihelm[vecPhihelm.size()-1-i]);
  }
  }



  // put decay vertex position
   tracklength_hel.push_back(0); //

 int nbHelTurn = 0;
 for(;;){
   double Zhel;
 
 for (unsigned int i=0; i<vecPhihel.size(); ++i){

   double phi = vecPhihel[i]+ nbHelTurn*2.*TMath::Pi();

   Zhel = GetHelPar(3,A,B,C,phi,Xc,Yc,tlam);

   std::cout<< " Zhel = " <<  Zhel << " ZhelMin = " << ZhelMin<< " ZhelMax = "<<ZhelMax << std::endl;

   if( Zhel< ZhelMin ||  Zhel > ZhelMax) break;

   vecPhi.push_back(phi);

   // tracklength / R_ * sin(theta_) = phi
   tracklength_hel.push_back( phi* R_ / sin(theta_));

   TVector3 posd(GetHelPar(1,A,B,C,phi,Xc,Yc,tlam),GetHelPar(2,A,B,C,phi,Xc,Yc,tlam),Zhel);
     vecPos.push_back(posd);

     std::cout<< " i= " <<i <<" vecPhihel[i] = " <<vecPhihel[i]<< " vecPhi = " <<vecPhi[vecPhi.size()-1]<< " pos Z =" <<vecPos[vecPhi.size()-1].Z() << " tracklength_hel  = " <<   tracklength_hel[i]  << std::endl;
 // vecPos[vecPos.size()-1].Print();
  }

 if( Zhel< ZhelMin ||  Zhel > ZhelMax) break;
   nbHelTurn =nbHelTurn +1;

 }

 //  return vecPhi;
 return tracklength_hel;
}


double  HelixTrackModel::GetHelPar(int i , double A, double B, double C, double phi, double Xc, double Yc , double tlam ){
  double cosphim = cos(phi);
  double sinphim = sin(phi);


  if(i==1) return  A*cosphim - B*sinphim + Xc - cosphim*Xc + sinphim*Yc;  
  if(i==2) return   B*cosphim + A*sinphim - sinphim*Xc + Yc - cosphim*Yc;
  if(i==3)return   C + TMath::Sqrt((A-Xc)* (A-Xc)+(B-Yc)*(B-Yc))* phi*tlam;
  return 99999999999.;
  }




} /* End of namespace genfit */
