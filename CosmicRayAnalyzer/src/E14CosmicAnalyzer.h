#ifdef E14COSMICANALYZER_H_
#define E14COSMICANALYZER_H_

#include "TGraph.h"
#include "E14CosmicHough.h"
#include "E14CosmicChisq.h"
#include "TTree.h"

class E14CosmicAnalyzer{
 private:
  bool b_hough;
  bool b_chisq;
  double d_chisqRoh;
  double d_houghRoh;
  double d_chisqTheta;
  double d_houghTheta;

 public:
  E14CosmicHough* hough;
  E14CosmicChisq* chisq;
  TGraph* grCosmic; 

  E14CosmicAnalyzer();
  virtual ~E14CosmicAnalyzer();
  bool Init();
  bool Reset();
  bool GetResult(TGraph* gr, double& roh, double& theta );
  bool GetResult(TGraph* gr, double& roh, double& theta , double& roh1, double& theta1 );

};



