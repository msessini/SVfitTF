#include "../interface/CrystalBall.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include <iostream>
#include <string>
#include <assert.h> 
#include <math.h>

namespace
{
  std::string findFile(const std::string& fileName)
  {
    edm::FileInPath inputFile(fileName);
    if ( inputFile.fullPath() == "" ) {
      std::cerr << "Error: Cannot find file = " << fileName << " !!" << std::endl;
      assert(0);
    }
    return inputFile.fullPath().data();
  }
}

CrystalBall::CrystalBall(const std::string& inputFilePath)
{ 
  std::cout << "<CrystalBall::CrystalBall>:" << std::endl;

  decayMode_ = 0;
  nDecayMode_ = 5;
  cbParSize_ = 7;
  xx_ = new double[1];
  pp_ = new double[cbParSize_];

  std::vector<std::string> parDir {
    inputFilePath + "/par0/", 
    inputFilePath + "/par1/", 
    inputFilePath + "/par2/", 
    inputFilePath + "/par3/", 
    inputFilePath + "/par4/", 
    inputFilePath + "/par5/", 
    inputFilePath + "/par6/", 
    inputFilePath + "/par7/"
  };

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

  for ( int decayMode = 0; decayMode < nDecayMode_; ++decayMode ){
    for ( int par = 0; par < cbParSize_; ++par ){
      std::string l2JetParFileName = findFile(parDir[par] + fileL2[decayMode]); // each parameter has it own directory
      std::string l3JetParFileName = findFile(parDir[par] + fileL3[decayMode]);
      std::cout << "reading tauES correction parameters fpr decayMode = " << decayMode << ", par #" << par << " from files:" << std::endl;
      std::cout << " L2 = " << l2JetParFileName << std::endl;
      std::cout << " L3 = " << l3JetParFileName << std::endl;
      mapPar_[decayMode][par] = new CrystalBallPar(decayMode, par, l2JetParFileName, l3JetParFileName);
    }
  }
}

CrystalBall::~CrystalBall()
{
  for ( int decayMode = 0; decayMode < nDecayMode_; ++decayMode ){
    for ( int par = 0; par < cbParSize_; ++par ){
      delete mapPar_[decayMode][par];
    }
  }
} 

void CrystalBall::setDecayMode(int decayMode)
{ 
  decayMode_ = decayMode; 

  if ( mapPar_.find(decayMode_) == mapPar_.end() ) { 
    std::cerr << "Error: Invalid decay mode = " << decayMode_ << " !!"; 
    assert(0); 
  }

  for ( int iPar = 0; iPar < cbParSize_; ++iPar ) {
    thePar_.push_back(mapPar_[decayMode_][iPar]);
  }
}

namespace
{
  double fnc_dscb(double* xx, double* pp) 
  { 
    double x   = xx[0];
    // gaussian core
    double N   = pp[0]; // norm
    double mu  = pp[1]; // mean
    double sig = pp[2]; // variance
    // transition parameters
    double a1  = pp[3];
    double p1  = pp[4];
    double a2  = pp[5];
    double p2  = pp[6];
    //
    double u   = (x - mu)/sig;
    double A1  = pow(p1/fabs(a1), p1)*exp(-a1*a1/2);
    double A2  = pow(p2/fabs(a2), p2)*exp(-a2*a2/2);
    double B1  = p1/fabs(a1) - fabs(a1);
    double B2  = p2/fabs(a2) - fabs(a2);
    //
    double result = N;
    if      ( u < -a1 ) result *= A1*pow(B1 - u, -p1);
    else if ( u <  a2 ) result *= exp(-u*u/2);
    else                result *= A2*pow(B2 + u, -p2);
    return result;
  }
}

double CrystalBall::operator()(double recPt, double genPt, double genEta) const
{ 
  (*xx_) = recPt/genPt; 

  for ( int iPar = 0; iPar < cbParSize_; ++iPar ) {
    pp_[iPar] = (*thePar_[iPar])(genPt, genEta);
  }

  return fnc_dscb(xx_, pp_);
}
