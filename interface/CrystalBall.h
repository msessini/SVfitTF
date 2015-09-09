// This class returns the CrystalBall probability 
// value for parameters produced with the JetMET packadge. 
// Author: Betty Calpas, Christian Veelken
// Email: betty.calpas@cern.ch, christian.veelken@cern.ch

#ifndef CRYSTAL_BALL_H
#define CRYSTAL_BALL_H

#include "CrystalBallPar.h"
#include <string>
#include <map>
#include <vector>

class CrystalBall
{
 public:
  CrystalBall(const std::string& = "TauAnalysis/SVfitTF/data/L2L3Corr/"); // default constructor
  virtual ~CrystalBall(); // default destructor
  double operator()(double recPt, double genPt, double genEta) const; // call the function "fnc_dscb" which return the CrystalBall probalility value

  void setDecayMode(int decayMode);

  const CrystalBallPar* getPar(int par); 

 private:
  int cbParSize_;  // nb of CrystalBall parameter
  int nDecayMode_; // nb of different supported tau decay modes
  enum { kAll, kOneProng0Pi0, kOneProng1Pi0, kOneProng2Pi0, kThreeProng0Pi0 };
  int decayMode_;  // reconstructed decay mode of given tau (to be set for each event by calling "setDecayMode" function)
  mutable double* xx_;
  mutable double* pp_; // point to the 7 parameter function value 
  std::vector<CrystalBallPar*> thePar_;
  std::map<int, std::vector<CrystalBallPar*> > mapPar_; // <decayMode, parametersVec>
  mutable double genPt_cache_;
  mutable double genEta_cache_;
};
#endif
