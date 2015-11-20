#ifndef TauAnalysis_SVfitTF_HadTauTFfromTGraph_h
#define TauAnalysis_SVfitTF_HadTauTFfromTGraph_h

/** \class HadTauTFfromTGraph
 *
 * Read transfer functions for hadronic tau pT as TGraph objects from ROOT file.
 *
 * \author Christian Veelken, Tallinn
 *
 */

#include "TauAnalysis/SVfitTF/interface/HadTauTFBase.h"

#include <TGraph.h>

#include <map>

class HadTauTFfromTGraph : public HadTauTFBase
{
 public:
  HadTauTFfromTGraph();
  HadTauTFfromTGraph(const TGraph* resolution);
  HadTauTFfromTGraph(const std::map<int, const TGraph*>& resolutionMap);
  ~HadTauTFfromTGraph();
  double operator()(double recPt, double genPt, double genEta) const;

  void setDecayMode(int decayMode) const;

  HadTauTFfromTGraph* Clone(const std::string& label) const;

 private:
  std::map<int, const TGraph*> resolutionMap_; // key = decayMode
  bool hasResolutionMap_;
  mutable int theDecayMode_;
  mutable const TGraph* theResolution_;
  mutable double xMin_; // range of validity of TGraph object (x-coordinate of first and last point)
  mutable double xMax_;
};

#endif
