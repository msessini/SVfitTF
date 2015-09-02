#include "../interface/CrystalBallPar.h"
#include <iostream>

using namespace std;

CrystalBallPar::CrystalBallPar( int parNumber, int decayMode, const std::string &jetParFileName ): parNumber_{parNumber}, decayMode_{decayMode}, jetParFileName_{jetParFileName}{
 
  cout << "\nStarting class instantiation for 1 correction only!!" << endl;
  
  // Create the JetCorrectorParameter objects, the order does not matter.
  cout << "Creating the JetCorrectorParameter objects..." << endl; 
  jetPar_  = new JetCorrectorParameters(jetParFileName_.data());

  // Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
  cout << "Load the JetCorrectorParameter objects into a vector..." << endl;
  vPar_.push_back(*jetPar_);

  // Construct a FactorizedJetCorrector object with the vector
  cout << "Construct a FactorizedJetCorrector object with the vector..." << endl;
  jetCorrector_ = new FactorizedJetCorrector(vPar_);

  cout << "Finished class instantiation for 1 correction level only!!\n" << endl;
}


CrystalBallPar::CrystalBallPar( int parNumber, int decayMode, const std::string &l2JetParFileName, const std::string &l3JetParFileName ): parNumber_{parNumber}, decayMode_{decayMode}, l2JetParFileName_{l2JetParFileName}, l3JetParFileName_{l3JetParFileName} {
 
  cout << "\nStarting class instantiation for L2L3!!" << endl;
  
  // Create the JetCorrectorParameter objects, the order does not matter.
  cout << "Creating the JetCorrectorParameter objects..." << endl; 
  l2JetPar_  = new JetCorrectorParameters(l2JetParFileName_.data());
  l3JetPar_  = new JetCorrectorParameters(l3JetParFileName_.data());

  // Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
  cout << "Load the JetCorrectorParameter objects into a vector..." << endl;
  vPar_.push_back(*l2JetPar_);
  vPar_.push_back(*l3JetPar_);

  // Construct a FactorizedJetCorrector object with the vector
  cout << "Construct a FactorizedJetCorrector object with the vector..." << endl;
  jetCorrector_ = new FactorizedJetCorrector(vPar_);

  cout << "Finished class instantiation for L2L3!!\n" << endl;
}

CrystalBallPar::~CrystalBallPar(){}

//double CrystalBallPar::operator()(double genPt, double genEta=1.) const{
double CrystalBallPar::operator()(double genPt, double genEta=1.) const{

  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCorFWLite

  // Set the necessary values, for every jet, inside the jet loop
  jetCorrector_->setJetPt(genPt);
  jetCorrector_->setJetEta(genEta);

  // Get a vector of the individual correction factors
  double correction = jetCorrector_->getCorrection();

  return correction;
}
