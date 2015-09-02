#ifndef CRYSTAL_BALL_PAR_H
#define CRYSTAL_BALL_PAR_H

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include <string>
#include <vector>

class CrystalBallPar
{
  public:
    CrystalBallPar(int parNumber, int decayMode, const std::string &jetParFileName); // implement 1 correction level at a time
    CrystalBallPar(int parNumber, int decayMode, const std::string &l2JetParFileName, const std::string &l3JetParFileName); // 2 corr level
    virtual ~CrystalBallPar();
    double operator()(double genPt, double genEta) const;

  private:
    int parNumber_;
    int decayMode_;
    std::string jetParFileName_;
    std::string l2JetParFileName_;
    std::string l3JetParFileName_;
    JetCorrectorParameters *jetPar_;
    JetCorrectorParameters *l2JetPar_;
    JetCorrectorParameters *l3JetPar_;
    FactorizedJetCorrector *jetCorrector_;
    std::vector<JetCorrectorParameters> vPar_;
};

#endif
