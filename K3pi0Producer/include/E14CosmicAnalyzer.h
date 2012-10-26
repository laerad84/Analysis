#ifndef E14COSMICANALYZER__H__
#define E14COSMICANALYZER__H__

#include "TGraph.h"

#include "HoughCsI.h"
#include "Chisq_cosmic.h"

class E14CosmicAnalyzer {
 private:
  bool mb_hough;
  bool mb_chisq;
  double md_chi2Roh;
  double md_chi2Theta;
  double md_houghRoh;
  double md_houghTheta;
 public:

  HoughCsI*     mc_hough;
  Chisq_cosmic* mc_chi2Cosmic;
  TGraph* m_gr;
  
  E14CosmicAnalyzer();
  virtual ~E14CosmicAnalyzer();
  
  bool Init();
  bool Reset();
  bool GetResult( TGraph*, double& , double& );
  bool GetResult( TGraph*, double& , double& , double&, double& );

};


#endif

