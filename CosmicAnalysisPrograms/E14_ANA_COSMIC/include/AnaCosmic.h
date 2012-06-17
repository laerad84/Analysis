#ifndef __ANACOSMIC__H__
#define __ANACOSMIC__H__

#include "CsIImage.h"
#include "HoughCsI.h"
#include "IDHandler.h"
#include "CsIImage.h"
#include "TH1.h"

class AnaCosmic{

 public:
  static const double sHisXminimum;
  static const double sHisXmaximum;
  static const double sHisNbin;

  
  E14ReadSumFile* m_Reader;
  IDHandler*      m_handler;
  CsIImage*       m_GainImage;
  CsIImage*       m_RmsImage;
  CsIImage*       m_nHitImage;
  TGraph*         m_GainGraph;
  TGraph*         m_RMSGraph;
  TGraph*         m_nHitGraph;

  TH1D*           m_MIPpeakHIS[N_TOTAL_CSI];
  
  double m_Trigger_threshold[N_TOTAL_Cosmic];
  int    m_Trigger_Hit[N_TOTAL_Cosmic];
  bool   m_Trigger_Up;
  bool   m_Trigger_Dn;

  TChain* m_chIn;
  TChain* m_chOut;


  



#endif

