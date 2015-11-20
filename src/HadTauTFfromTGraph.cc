#include "TauAnalysis/SVfitTF/interface/HadTauTFfromTGraph.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include <iostream>
#include <string>
#include <assert.h> 
#include <math.h>

namespace
{
  void initializeRange(const TGraph* graph, double& xMin, double& xMax)
  {
    int numPoints = graph->GetN();
    assert(numPoints >= 1);
    for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
      double x, y;
      graph->GetPoint(idxPoint, x, y);
      if ( idxPoint == 0 || x < xMin ) xMin = x;
      if ( idxPoint == 0 || x > xMax ) xMax = x;
    }
  }
}

HadTauTFfromTGraph::HadTauTFfromTGraph()
  : hasResolutionMap_(false),
    theDecayMode_(-1),
    theResolution_(0),
    xMin_(-1.),
    xMax_(-1.)
{}

HadTauTFfromTGraph::HadTauTFfromTGraph(const TGraph* resolution)
  : hasResolutionMap_(false),
    theDecayMode_(-1),
    theResolution_(resolution)
{ 
  resolutionMap_[-1] = theResolution_;
  initializeRange(theResolution_, xMin_, xMax_);
}

HadTauTFfromTGraph::HadTauTFfromTGraph(const std::map<int, const TGraph*>& resolutionMap)
  : resolutionMap_(resolutionMap),
    hasResolutionMap_(true),
    theDecayMode_(-1),
    theResolution_(0)
{}

HadTauTFfromTGraph::~HadTauTFfromTGraph()
{
  // CV: assume that TGraph objects are owned by calling code
  //     and do not delete them
}

void HadTauTFfromTGraph::setDecayMode(int decayMode) const
{
  if ( hasResolutionMap_ ) {
    if ( resolutionMap_.find(decayMode) == resolutionMap_.end() ) {
      std::cerr << "No tau pT transfer functions defined for decayMode = " << decayMode << " !!" << std::endl;
      assert(0);
    }
    theResolution_ = resolutionMap_.find(decayMode)->second;
    initializeRange(theResolution_, xMin_, xMax_);
  }
}

double HadTauTFfromTGraph::operator()(double recPt, double genPt, double genEta) const
{
  if ( !theResolution_ ) {
    std::cerr << "No tau pT transfer functions defined, call 'setDecayMode' function first !!" << std::endl;
    assert(0);
  }
  if ( genPt > 0. ) {
    double x = recPt/genPt;
    if ( x < 0. ) return 0.;
    if ( x < xMin_ ) x = xMin_;
    if ( x > xMax_ ) x = xMax_;
    return theResolution_->Eval(x);
  } else {
    return 0.;
  }
}

HadTauTFfromTGraph* HadTauTFfromTGraph::Clone(const std::string& label) const
{
  HadTauTFfromTGraph* clone = new HadTauTFfromTGraph();
  for ( std::map<int, const TGraph*>::const_iterator resolutionEntry = resolutionMap_.begin();
	resolutionEntry != resolutionMap_.end(); ++resolutionEntry ) {
    int decayMode = resolutionEntry->first;
    const TGraph* resolution = resolutionEntry->second;
    const TGraph* resolution_cloned = (TGraph*)resolution->Clone(Form("%s_%s", resolution->GetName(), label.data()));
    clone->resolutionMap_[decayMode] = resolution_cloned;
  }
  clone->hasResolutionMap_ = hasResolutionMap_;
  clone->theDecayMode_ = theDecayMode_;
  if ( clone->resolutionMap_.find(theDecayMode_) != clone->resolutionMap_.end() ) {
    clone->theResolution_ = clone->resolutionMap_.find(theDecayMode_)->second;
    initializeRange(clone->theResolution_, clone->xMin_, clone->xMax_);
  } else {
    clone->theResolution_ = 0;
  }  
  return clone;
}
