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
    CrystalBall(); // default constructor
    virtual ~CrystalBall(); // default destructor
    //double operator()(double recPt, double genPt, double genEta) const; // call the function "fnc_dscb" which return the CrystalBall probalility value
    double operator()(double recPt, double genPt, double genEta); // call the function "fnc_dscb" which return the CrystalBall probalility value

    void setDecayMode(int decayMode);

    const CrystalBallPar * getPar(int par); 

  private:
    int cbParSize_; // nb of CrystalBall parameter
    int decayMode_; // tau decay mode recontruction
    int nDecayMode_; 
    double *xx_;
    double *pp_; // point to the 7 parameter function value 
    //double fnc_dscb(double *xx, double *pp) const; // return the CrystalBall probalility value
    double fnc_dscb(double *xx, double *pp); // return the CrystalBall probalility value
    std::vector<CrystalBallPar*> thePar_;
    std::map<int, std::vector<CrystalBallPar*> > mapPar_; // <decayMode, parametersVec>
};
#endif
