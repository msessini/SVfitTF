#include "../interface/CrystalBall.h"
#include <iostream>
#include <string>
#include <assert.h> 
#include <math.h>

CrystalBall::CrystalBall(){ 

  decayMode_=0;
  nDecayMode_=5;
  cbParSize_=7;
  xx_=new double[1];
  pp_=new double[cbParSize_];

  std::vector<std::string> parDir{"../L2L3Corr/par0/", "../L2L3Corr/par1/", "../L2L3Corr/par2/", "../L2L3Corr/par3/", "../L2L3Corr/par4/", "../L2L3Corr/par5/", "../L2L3Corr/par6/"};

  std::vector<std::string> fileL2 {
      "TauJec11V1_L2Relative_AK5tauHPSlooseCombDBcorrAll.txt",
      "TauJec11V1_L2Relative_AK5tauHPSlooseCombDBcorrOneProng0Pi0.txt",
      "TauJec11V1_L2Relative_AK5tauHPSlooseCombDBcorrOneProng1Pi0.txt",
      "TauJec11V1_L2Relative_AK5tauHPSlooseCombDBcorrOneProng2Pi0.txt",
      "TauJec11V1_L2Relative_AK5tauHPSlooseCombDBcorrThreeProng0Pi0.txt",
      "TauJec11V1_L2Relative_AK5tauHPSmediumCombDBcorrAll.txt",
      "TauJec11V1_L2Relative_AK5tauHPSmediumCombDBcorrOneProng0Pi0.txt",
      "TauJec11V1_L2Relative_AK5tauHPSmediumCombDBcorrOneProng1Pi0.txt",
      "TauJec11V1_L2Relative_AK5tauHPSmediumCombDBcorrOneProng2Pi0.txt",
      "TauJec11V1_L2Relative_AK5tauHPSmediumCombDBcorrThreeProng0Pi0.txt",
      "TauJec11V1_L2Relative_AK5tauHPStightCombDBcorrAll.txt",
      "TauJec11V1_L2Relative_AK5tauHPStightCombDBcorrOneProng0Pi0.txt",
      "TauJec11V1_L2Relative_AK5tauHPStightCombDBcorrOneProng1Pi0.txt",
      "TauJec11V1_L2Relative_AK5tauHPStightCombDBcorrOneProng2Pi0.txt",
      "TauJec11V1_L2Relative_AK5tauHPStightCombDBcorrThreeProng0Pi0.txt"
  };

  std::vector<std::string> fileL3 {
      "TauJec11V1_L3Absolute_AK5tauHPSlooseCombDBcorrAll.txt",
      "TauJec11V1_L3Absolute_AK5tauHPSlooseCombDBcorrOneProng0Pi0.txt",
      "TauJec11V1_L3Absolute_AK5tauHPSlooseCombDBcorrOneProng1Pi0.txt",
      "TauJec11V1_L3Absolute_AK5tauHPSlooseCombDBcorrOneProng2Pi0.txt",
      "TauJec11V1_L3Absolute_AK5tauHPSlooseCombDBcorrThreeProng0Pi0.txt",
      "TauJec11V1_L3Absolute_AK5tauHPSmediumCombDBcorrAll.txt",
      "TauJec11V1_L3Absolute_AK5tauHPSmediumCombDBcorrOneProng0Pi0.txt",
      "TauJec11V1_L3Absolute_AK5tauHPSmediumCombDBcorrOneProng1Pi0.txt",
      "TauJec11V1_L3Absolute_AK5tauHPSmediumCombDBcorrOneProng2Pi0.txt",
      "TauJec11V1_L3Absolute_AK5tauHPSmediumCombDBcorrThreeProng0Pi0.txt",
      "TauJec11V1_L3Absolute_AK5tauHPStightCombDBcorrAll.txt",
      "TauJec11V1_L3Absolute_AK5tauHPStightCombDBcorrOneProng0Pi0.txt",
      "TauJec11V1_L3Absolute_AK5tauHPStightCombDBcorrOneProng1Pi0.txt",
      "TauJec11V1_L3Absolute_AK5tauHPStightCombDBcorrOneProng2Pi0.txt",
      "TauJec11V1_L3Absolute_AK5tauHPStightCombDBcorrThreeProng0Pi0.txt"
  };


  for(int decay=0; decay<nDecayMode_; decay++){
    for(int par=0; par<cbParSize_; par++){
      mapPar_[decay][par]=new CrystalBallPar(decay, par, parDir[par]+fileL2[decay], parDir[par]+fileL3[decay]); // each parameter has it own directory
    }
  }

}

CrystalBall::~CrystalBall(){} 

void CrystalBall::setDecayMode( int decayMode ){ 

  decayMode_=decayMode; 

  if ( mapPar_.find(decayMode_) == mapPar_.end() ) { std::cerr << "Error: Invalid decay mode = " << decayMode_ << " !!"; assert(0); }

  for ( int iPar = 0; iPar < cbParSize_; ++iPar ) thePar_.push_back( mapPar_[decayMode_][iPar] );

}

//double CrystalBall::operator()(double recPt, double genPt, double genEta) const{ 
double CrystalBall::operator()(double recPt, double genPt, double genEta){ 

  *xx_=recPt/genPt; 

  for ( int iPar = 0; iPar < cbParSize_; ++iPar ) pp_[iPar] = (*thePar_[iPar])(genPt, genEta);

  return fnc_dscb( xx_, pp_);
}


//double CrystalBall::fnc_dscb(double *xx,double *pp) const{
double CrystalBall::fnc_dscb(double *xx,double *pp) {
  double x   = xx[0];
  // gaussian core
  double N   = pp[0];//norm
  double mu  = pp[1];//mean
  double sig = pp[2];//variance
  // transition parameters
  double a1  = pp[3];
  double p1  = pp[4];
  double a2  = pp[5];
  double p2  = pp[6];
  //
  double u   = (x-mu)/sig;
  double A1  = pow(p1/fabs(a1),p1)*exp(-a1*a1/2);
  double A2  = pow(p2/fabs(a2),p2)*exp(-a2*a2/2);
  double B1  = p1/fabs(a1) - fabs(a1);
  double B2  = p2/fabs(a2) - fabs(a2);
  //
  double result(N);
  if      (u<-a1) result *= A1*pow(B1-u,-p1);
  else if (u<a2)  result *= exp(-u*u/2);
  else            result *= A2*pow(B2+u,-p2);
  return result;
}

