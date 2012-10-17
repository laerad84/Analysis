//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 17 22:32:47 2012 by ROOT version 5.30/04
// from TTree WFTree/Waveform Analyzed Tree
// found on file: /Volume0/ExpData/2012_Feb_Beam/RootFile_wav/TEMPLATE_FIT_RESULT_1_4503.root
//////////////////////////////////////////////////////////

#ifndef E14WavReader_h
#define E14WavReader_h
#include <iostream>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class E14WavReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           RunNo;
   Int_t           EventNo;
   Int_t           TrigFlag;
   Int_t           CosmicTrig;
   Int_t           LaserTrig;
   Int_t           CVTrig;
   Int_t           CosmicTrigFlagUp;
   Int_t           CosmicTrigFlagDn;
   Int_t           CVTrigFlag;
   Int_t           LasertrigFlag;
   Double_t        TimePeak;
   Double_t        TimeSigma;
   Int_t           CsiNumber;
   Short_t         CsiID[2716];   //[CsiNumber]
   Short_t         CsiFitHeight[2716];   //[CsiNumber]
   Short_t         CsiFitTime[2716];   //[CsiNumber]
   Double_t        CsiChisq[2716];   //[CsiNumber]
   Short_t         CsiNDF[2716];   //[CsiNumber]
   Double_t        CsiFitShape[2716];   //[CsiNumber]
   Double_t        CsiDeltaDiff[2716];   //[CsiNumber]
   Double_t        CsiPedestal[2716];   //[CsiNumber]
   Double_t        CsiSignal[2716];   //[CsiNumber]
   Double_t        CsiTime[2716];   //[CsiNumber]
   Double_t        CsiHHTime[2716];   //[CsiNumber]
   Double_t        CsiParA[2716];   //[CsiNumber]
   Double_t        CsiParB[2716];   //[CsiNumber]
   Double_t        CsiADC[2716];   //[CsiNumber]
   Double_t        CsiFitADC[2716];   //[CsiNumber]
   Int_t           CC03Number;
   Short_t         CC03ID[32];   //[CC03Number]
   Short_t         CC03FitHeight[32];   //[CC03Number]
   Short_t         CC03FitTime[32];   //[CC03Number]
   Double_t        CC03Chisq[32];   //[CC03Number]
   Short_t         CC03NDF[32];   //[CC03Number]
   Double_t        CC03FitShape[32];   //[CC03Number]
   Double_t        CC03DeltaDiff[32];   //[CC03Number]
   Double_t        CC03Pedestal[32];   //[CC03Number]
   Double_t        CC03Signal[32];   //[CC03Number]
   Double_t        CC03Time[32];   //[CC03Number]
   Double_t        CC03HHTime[32];   //[CC03Number]
   Double_t        CC03ParA[32];   //[CC03Number]
   Double_t        CC03ParB[32];   //[CC03Number]
   Double_t        CC03ADC[32];   //[CC03Number]
   Double_t        CC03FitADC[32];   //[CC03Number]
   Int_t           OEVNumber;
   Short_t         OEVID[44];   //[OEVNumber]
   Short_t         OEVFitHeight[44];   //[OEVNumber]
   Short_t         OEVFitTime[44];   //[OEVNumber]
   Double_t        OEVChisq[44];   //[OEVNumber]
   Short_t         OEVNDF[44];   //[OEVNumber]
   Double_t        OEVFitShape[44];   //[OEVNumber]
   Double_t        OEVDeltaDiff[44];   //[OEVNumber]
   Double_t        OEVPedestal[44];   //[OEVNumber]
   Double_t        OEVSignal[44];   //[OEVNumber]
   Double_t        OEVTime[44];   //[OEVNumber]
   Double_t        OEVHHTime[44];   //[OEVNumber]
   Double_t        OEVParA[44];   //[OEVNumber]
   Double_t        OEVParB[44];   //[OEVNumber]
   Double_t        OEVADC[44];   //[OEVNumber]
   Double_t        OEVFitADC[44];   //[OEVNumber]
   Int_t           CVNumber;
   Short_t         CVID[10];   //[CVNumber]
   Short_t         CVFitHeight[10];   //[CVNumber]
   Short_t         CVFitTime[10];   //[CVNumber]
   Double_t        CVChisq[10];   //[CVNumber]
   Short_t         CVNDF[10];   //[CVNumber]
   Double_t        CVFitShape[10];   //[CVNumber]
   Double_t        CVDeltaDiff[10];   //[CVNumber]
   Double_t        CVPedestal[10];   //[CVNumber]
   Double_t        CVSignal[10];   //[CVNumber]
   Double_t        CVTime[10];   //[CVNumber]
   Double_t        CVHHTime[10];   //[CVNumber]
   Double_t        CVParA[10];   //[CVNumber]
   Double_t        CVParB[10];   //[CVNumber]
   Double_t        CVADC[10];   //[CVNumber]
   Double_t        CVFitADC[10];   //[CVNumber]
   Int_t           CosmicNumber;
   Short_t         CosmicID[20];   //[CosmicNumber]
   Short_t         CosmicFitHeight[20];   //[CosmicNumber]
   Short_t         CosmicFitTime[20];   //[CosmicNumber]
   Double_t        CosmicChisq[20];   //[CosmicNumber]
   Short_t         CosmicNDF[20];   //[CosmicNumber]
   Double_t        CosmicFitShape[20];   //[CosmicNumber]
   Double_t        CosmicDeltaDiff[20];   //[CosmicNumber]
   Double_t        CosmicPedestal[20];   //[CosmicNumber]
   Double_t        CosmicSignal[20];   //[CosmicNumber]
   Double_t        CosmicTime[20];   //[CosmicNumber]
   Double_t        CosmicHHTime[20];   //[CosmicNumber]
   Double_t        CosmicParA[20];   //[CosmicNumber]
   Double_t        CosmicParB[20];   //[CosmicNumber]
   Double_t        CosmicADC[20];   //[CosmicNumber]
   Double_t        CosmicFitADC[20];   //[CosmicNumber]
   Int_t           LaserNumber;
   Short_t         LaserID[5];   //[LaserNumber]
   Short_t         LaserFitHeight[5];   //[LaserNumber]
   Short_t         LaserFitTime[5];   //[LaserNumber]
   Double_t        LaserChisq[5];   //[LaserNumber]
   Short_t         LaserNDF[5];   //[LaserNumber]
   Double_t        LaserFitShape[5];   //[LaserNumber]
   Double_t        LaserDeltaDiff[5];   //[LaserNumber]
   Double_t        LaserPedestal[5];   //[LaserNumber]
   Double_t        LaserSignal[5];   //[LaserNumber]
   Double_t        LaserTime[5];   //[LaserNumber]
   Double_t        LaserHHTime[5];   //[LaserNumber]
   Double_t        LaserParA[5];   //[LaserNumber]
   Double_t        LaserParB[5];   //[LaserNumber]
   Double_t        LaserADC[5];   //[LaserNumber]
   Double_t        LaserFitADC[5];   //[LaserNumber]
   Int_t           EtcNumber;
   Short_t         EtcID[10];   //[EtcNumber]
   Short_t         EtcFitHeight[10];   //[EtcNumber]
   Short_t         EtcFitTime[10];   //[EtcNumber]
   Double_t        EtcChisq[10];   //[EtcNumber]
   Short_t         EtcNDF[10];   //[EtcNumber]
   Double_t        EtcFitShape[10];   //[EtcNumber]
   Double_t        EtcDeltaDiff[10];   //[EtcNumber]
   Double_t        EtcPedestal[10];   //[EtcNumber]
   Double_t        EtcSignal[10];   //[EtcNumber]
   Double_t        EtcTime[10];   //[EtcNumber]
   Double_t        EtcHHTime[10];   //[EtcNumber]
   Double_t        EtcParA[10];   //[EtcNumber]
   Double_t        EtcParB[10];   //[EtcNumber]
   Double_t        EtcADC[10];   //[EtcNumber]
   Double_t        EtcFitADC[10];   //[EtcNumber]

   // List of branches
   TBranch        *b_RunNo;   //!
   TBranch        *b_EventNo;   //!
   TBranch        *b_TrigFlag;   //!
   TBranch        *b_CosmicTrig;   //!
   TBranch        *b_LaserTrig;   //!
   TBranch        *b_CVTrig;   //!
   TBranch        *b_CosmicTrigFlagUp;   //!
   TBranch        *b_CosmicTrigFlagDn;   //!
   TBranch        *b_CVTrigFlag;   //!
   TBranch        *b_LaserTrigFlag;   //!
   TBranch        *b_TimePeak;   //!
   TBranch        *b_TimeSigma;   //!
   TBranch        *b_CsiNumber;   //!
   TBranch        *b_CsiID;   //!
   TBranch        *b_CsiFitHeight;   //!
   TBranch        *b_CsiFitTime;   //!
   TBranch        *b_CsiChisq;   //!
   TBranch        *b_CsiNDF;   //!
   TBranch        *b_CsiFitShape;   //!
   TBranch        *b_CsiDeltaDiff;   //!
   TBranch        *b_CsiPedestal;   //!
   TBranch        *b_CsiSignal;   //!
   TBranch        *b_CsiTime;   //!
   TBranch        *b_CsiHHTime;   //!
   TBranch        *b_CsiParA;   //!
   TBranch        *b_CsiParB;   //!
   TBranch        *b_CsiADC;   //!
   TBranch        *b_CsiFitADC;   //!
   TBranch        *b_CC03Number;   //!
   TBranch        *b_CC03ID;   //!
   TBranch        *b_CC03FitHeight;   //!
   TBranch        *b_CC03FitTime;   //!
   TBranch        *b_CC03Chisq;   //!
   TBranch        *b_CC03NDF;   //!
   TBranch        *b_CC03FitShape;   //!
   TBranch        *b_CC03DeltaDiff;   //!
   TBranch        *b_CC03Pedestal;   //!
   TBranch        *b_CC03Signal;   //!
   TBranch        *b_CC03Time;   //!
   TBranch        *b_CC03HHTime;   //!
   TBranch        *b_CC03ParA;   //!
   TBranch        *b_CC03ParB;   //!
   TBranch        *b_CC03ADC;   //!
   TBranch        *b_CC03FitADC;   //!
   TBranch        *b_OEVNumber;   //!
   TBranch        *b_OEVID;   //!
   TBranch        *b_OEVFitHeight;   //!
   TBranch        *b_OEVFitTime;   //!
   TBranch        *b_OEVChisq;   //!
   TBranch        *b_OEVNDF;   //!
   TBranch        *b_OEVFitShape;   //!
   TBranch        *b_OEVDeltaDiff;   //!
   TBranch        *b_OEVPedestal;   //!
   TBranch        *b_OEVSignal;   //!
   TBranch        *b_OEVTime;   //!
   TBranch        *b_OEVHHTime;   //!
   TBranch        *b_OEVParA;   //!
   TBranch        *b_OEVParB;   //!
   TBranch        *b_OEVADC;   //!
   TBranch        *b_OEVFitADC;   //!
   TBranch        *b_CVNumber;   //!
   TBranch        *b_CVID;   //!
   TBranch        *b_CVFitHeight;   //!
   TBranch        *b_CVFitTime;   //!
   TBranch        *b_CVChisq;   //!
   TBranch        *b_CVNDF;   //!
   TBranch        *b_CVFitShape;   //!
   TBranch        *b_CVDeltaDiff;   //!
   TBranch        *b_CVPedestal;   //!
   TBranch        *b_CVSignal;   //!
   TBranch        *b_CVTime;   //!
   TBranch        *b_CVHHTime;   //!
   TBranch        *b_CVParA;   //!
   TBranch        *b_CVParB;   //!
   TBranch        *b_CVADC;   //!
   TBranch        *b_CVFitADC;   //!
   TBranch        *b_CosmicNumber;   //!
   TBranch        *b_CosmicID;   //!
   TBranch        *b_CosmicFitHeight;   //!
   TBranch        *b_CosmicFitTime;   //!
   TBranch        *b_CosmicChisq;   //!
   TBranch        *b_CosmicNDF;   //!
   TBranch        *b_CosmicFitShape;   //!
   TBranch        *b_CosmicDeltaDiff;   //!
   TBranch        *b_CosmicPedestal;   //!
   TBranch        *b_CosmicSignal;   //!
   TBranch        *b_CosmicTime;   //!
   TBranch        *b_CosmicHHTime;   //!
   TBranch        *b_CosmicParA;   //!
   TBranch        *b_CosmicParB;   //!
   TBranch        *b_CosmicADC;   //!
   TBranch        *b_CosmicFitADC;   //!
   TBranch        *b_LaserNumber;   //!
   TBranch        *b_LaserID;   //!
   TBranch        *b_LaserFitHeight;   //!
   TBranch        *b_LaserFitTime;   //!
   TBranch        *b_LaserChisq;   //!
   TBranch        *b_LaserNDF;   //!
   TBranch        *b_LaserFitShape;   //!
   TBranch        *b_LaserDeltaDiff;   //!
   TBranch        *b_LaserPedestal;   //!
   TBranch        *b_LaserSignal;   //!
   TBranch        *b_LaserTime;   //!
   TBranch        *b_LaserHHTime;   //!
   TBranch        *b_LaserParA;   //!
   TBranch        *b_LaserParB;   //!
   TBranch        *b_LaserADC;   //!
   TBranch        *b_LaserFitADC;   //!
   TBranch        *b_EtcNumber;   //!
   TBranch        *b_EtcID;   //!
   TBranch        *b_EtcFitHeight;   //!
   TBranch        *b_EtcFitTime;   //!
   TBranch        *b_EtcChisq;   //!
   TBranch        *b_EtcNDF;   //!
   TBranch        *b_EtcFitShape;   //!
   TBranch        *b_EtcDeltaDiff;   //!
   TBranch        *b_EtcPedestal;   //!
   TBranch        *b_EtcSignal;   //!
   TBranch        *b_EtcTime;   //!
   TBranch        *b_EtcHHTime;   //!
   TBranch        *b_EtcParA;   //!
   TBranch        *b_EtcParB;   //!
   TBranch        *b_EtcADC;   //!
   TBranch        *b_EtcFitADC;   //!

   E14WavReader(TTree *tree=0);
   virtual ~E14WavReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef E14WavReader_cxx
E14WavReader::E14WavReader(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     std::cerr << " Tree is NULL " << std::endl;
   }
   Init(tree);
}

E14WavReader::~E14WavReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t E14WavReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t E14WavReader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void E14WavReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNo", &RunNo, &b_RunNo);
   fChain->SetBranchAddress("EventNo", &EventNo, &b_EventNo);
   fChain->SetBranchAddress("TrigFlag", &TrigFlag, &b_TrigFlag);
   fChain->SetBranchAddress("CosmicTrig", &CosmicTrig, &b_CosmicTrig);
   fChain->SetBranchAddress("LaserTrig", &LaserTrig, &b_LaserTrig);
   fChain->SetBranchAddress("CVTrig", &CVTrig, &b_CVTrig);
   fChain->SetBranchAddress("CosmicTrigFlagUp", &CosmicTrigFlagUp, &b_CosmicTrigFlagUp);
   fChain->SetBranchAddress("CosmicTrigFlagDn", &CosmicTrigFlagDn, &b_CosmicTrigFlagDn);
   fChain->SetBranchAddress("CVTrigFlag", &CVTrigFlag, &b_CVTrigFlag);
   fChain->SetBranchAddress("LasertrigFlag", &LasertrigFlag, &b_LaserTrigFlag);
   fChain->SetBranchAddress("TimePeak", &TimePeak, &b_TimePeak);
   fChain->SetBranchAddress("TimeSigma", &TimeSigma, &b_TimeSigma);
   fChain->SetBranchAddress("CsiNumber", &CsiNumber, &b_CsiNumber);
   fChain->SetBranchAddress("CsiID", CsiID, &b_CsiID);
   fChain->SetBranchAddress("CsiFitHeight", CsiFitHeight, &b_CsiFitHeight);
   fChain->SetBranchAddress("CsiFitTime", CsiFitTime, &b_CsiFitTime);
   fChain->SetBranchAddress("CsiChisq", CsiChisq, &b_CsiChisq);
   fChain->SetBranchAddress("CsiNDF", CsiNDF, &b_CsiNDF);
   fChain->SetBranchAddress("CsiFitShape", CsiFitShape, &b_CsiFitShape);
   fChain->SetBranchAddress("CsiDeltaDiff", CsiDeltaDiff, &b_CsiDeltaDiff);
   fChain->SetBranchAddress("CsiPedestal", CsiPedestal, &b_CsiPedestal);
   fChain->SetBranchAddress("CsiSignal", CsiSignal, &b_CsiSignal);
   fChain->SetBranchAddress("CsiTime", CsiTime, &b_CsiTime);
   fChain->SetBranchAddress("CsiHHTime", CsiHHTime, &b_CsiHHTime);
   fChain->SetBranchAddress("CsiParA", CsiParA, &b_CsiParA);
   fChain->SetBranchAddress("CsiParB", CsiParB, &b_CsiParB);
   fChain->SetBranchAddress("CsiADC", CsiADC, &b_CsiADC);
   fChain->SetBranchAddress("CsiFitADC", CsiFitADC, &b_CsiFitADC);
   fChain->SetBranchAddress("CC03Number", &CC03Number, &b_CC03Number);
   fChain->SetBranchAddress("CC03ID", CC03ID, &b_CC03ID);
   fChain->SetBranchAddress("CC03FitHeight", CC03FitHeight, &b_CC03FitHeight);
   fChain->SetBranchAddress("CC03FitTime", CC03FitTime, &b_CC03FitTime);
   fChain->SetBranchAddress("CC03Chisq", CC03Chisq, &b_CC03Chisq);
   fChain->SetBranchAddress("CC03NDF", CC03NDF, &b_CC03NDF);
   fChain->SetBranchAddress("CC03FitShape", CC03FitShape, &b_CC03FitShape);
   fChain->SetBranchAddress("CC03DeltaDiff", CC03DeltaDiff, &b_CC03DeltaDiff);
   fChain->SetBranchAddress("CC03Pedestal", CC03Pedestal, &b_CC03Pedestal);
   fChain->SetBranchAddress("CC03Signal", CC03Signal, &b_CC03Signal);
   fChain->SetBranchAddress("CC03Time", CC03Time, &b_CC03Time);
   fChain->SetBranchAddress("CC03HHTime", CC03HHTime, &b_CC03HHTime);
   fChain->SetBranchAddress("CC03ParA", CC03ParA, &b_CC03ParA);
   fChain->SetBranchAddress("CC03ParB", CC03ParB, &b_CC03ParB);
   fChain->SetBranchAddress("CC03ADC", CC03ADC, &b_CC03ADC);
   fChain->SetBranchAddress("CC03FitADC", CC03FitADC, &b_CC03FitADC);
   fChain->SetBranchAddress("OEVNumber", &OEVNumber, &b_OEVNumber);
   fChain->SetBranchAddress("OEVID", OEVID, &b_OEVID);
   fChain->SetBranchAddress("OEVFitHeight", OEVFitHeight, &b_OEVFitHeight);
   fChain->SetBranchAddress("OEVFitTime", OEVFitTime, &b_OEVFitTime);
   fChain->SetBranchAddress("OEVChisq", OEVChisq, &b_OEVChisq);
   fChain->SetBranchAddress("OEVNDF", OEVNDF, &b_OEVNDF);
   fChain->SetBranchAddress("OEVFitShape", OEVFitShape, &b_OEVFitShape);
   fChain->SetBranchAddress("OEVDeltaDiff", OEVDeltaDiff, &b_OEVDeltaDiff);
   fChain->SetBranchAddress("OEVPedestal", OEVPedestal, &b_OEVPedestal);
   fChain->SetBranchAddress("OEVSignal", OEVSignal, &b_OEVSignal);
   fChain->SetBranchAddress("OEVTime", OEVTime, &b_OEVTime);
   fChain->SetBranchAddress("OEVHHTime", OEVHHTime, &b_OEVHHTime);
   fChain->SetBranchAddress("OEVParA", OEVParA, &b_OEVParA);
   fChain->SetBranchAddress("OEVParB", OEVParB, &b_OEVParB);
   fChain->SetBranchAddress("OEVADC", OEVADC, &b_OEVADC);
   fChain->SetBranchAddress("OEVFitADC", OEVFitADC, &b_OEVFitADC);
   fChain->SetBranchAddress("CVNumber", &CVNumber, &b_CVNumber);
   fChain->SetBranchAddress("CVID", CVID, &b_CVID);
   fChain->SetBranchAddress("CVFitHeight", CVFitHeight, &b_CVFitHeight);
   fChain->SetBranchAddress("CVFitTime", CVFitTime, &b_CVFitTime);
   fChain->SetBranchAddress("CVChisq", CVChisq, &b_CVChisq);
   fChain->SetBranchAddress("CVNDF", CVNDF, &b_CVNDF);
   fChain->SetBranchAddress("CVFitShape", CVFitShape, &b_CVFitShape);
   fChain->SetBranchAddress("CVDeltaDiff", CVDeltaDiff, &b_CVDeltaDiff);
   fChain->SetBranchAddress("CVPedestal", CVPedestal, &b_CVPedestal);
   fChain->SetBranchAddress("CVSignal", CVSignal, &b_CVSignal);
   fChain->SetBranchAddress("CVTime", CVTime, &b_CVTime);
   fChain->SetBranchAddress("CVHHTime", CVHHTime, &b_CVHHTime);
   fChain->SetBranchAddress("CVParA", CVParA, &b_CVParA);
   fChain->SetBranchAddress("CVParB", CVParB, &b_CVParB);
   fChain->SetBranchAddress("CVADC", CVADC, &b_CVADC);
   fChain->SetBranchAddress("CVFitADC", CVFitADC, &b_CVFitADC);
   fChain->SetBranchAddress("CosmicNumber", &CosmicNumber, &b_CosmicNumber);
   fChain->SetBranchAddress("CosmicID", CosmicID, &b_CosmicID);
   fChain->SetBranchAddress("CosmicFitHeight", CosmicFitHeight, &b_CosmicFitHeight);
   fChain->SetBranchAddress("CosmicFitTime", CosmicFitTime, &b_CosmicFitTime);
   fChain->SetBranchAddress("CosmicChisq", CosmicChisq, &b_CosmicChisq);
   fChain->SetBranchAddress("CosmicNDF", CosmicNDF, &b_CosmicNDF);
   fChain->SetBranchAddress("CosmicFitShape", CosmicFitShape, &b_CosmicFitShape);
   fChain->SetBranchAddress("CosmicDeltaDiff", CosmicDeltaDiff, &b_CosmicDeltaDiff);
   fChain->SetBranchAddress("CosmicPedestal", CosmicPedestal, &b_CosmicPedestal);
   fChain->SetBranchAddress("CosmicSignal", CosmicSignal, &b_CosmicSignal);
   fChain->SetBranchAddress("CosmicTime", CosmicTime, &b_CosmicTime);
   fChain->SetBranchAddress("CosmicHHTime", CosmicHHTime, &b_CosmicHHTime);
   fChain->SetBranchAddress("CosmicParA", CosmicParA, &b_CosmicParA);
   fChain->SetBranchAddress("CosmicParB", CosmicParB, &b_CosmicParB);
   fChain->SetBranchAddress("CosmicADC", CosmicADC, &b_CosmicADC);
   fChain->SetBranchAddress("CosmicFitADC", CosmicFitADC, &b_CosmicFitADC);
   fChain->SetBranchAddress("LaserNumber", &LaserNumber, &b_LaserNumber);
   fChain->SetBranchAddress("LaserID", LaserID, &b_LaserID);
   fChain->SetBranchAddress("LaserFitHeight", LaserFitHeight, &b_LaserFitHeight);
   fChain->SetBranchAddress("LaserFitTime", LaserFitTime, &b_LaserFitTime);
   fChain->SetBranchAddress("LaserChisq", LaserChisq, &b_LaserChisq);
   fChain->SetBranchAddress("LaserNDF", LaserNDF, &b_LaserNDF);
   fChain->SetBranchAddress("LaserFitShape", LaserFitShape, &b_LaserFitShape);
   fChain->SetBranchAddress("LaserDeltaDiff", LaserDeltaDiff, &b_LaserDeltaDiff);
   fChain->SetBranchAddress("LaserPedestal", LaserPedestal, &b_LaserPedestal);
   fChain->SetBranchAddress("LaserSignal", LaserSignal, &b_LaserSignal);
   fChain->SetBranchAddress("LaserTime", LaserTime, &b_LaserTime);
   fChain->SetBranchAddress("LaserHHTime", LaserHHTime, &b_LaserHHTime);
   fChain->SetBranchAddress("LaserParA", LaserParA, &b_LaserParA);
   fChain->SetBranchAddress("LaserParB", LaserParB, &b_LaserParB);
   fChain->SetBranchAddress("LaserADC", LaserADC, &b_LaserADC);
   fChain->SetBranchAddress("LaserFitADC", LaserFitADC, &b_LaserFitADC);
   fChain->SetBranchAddress("EtcNumber", &EtcNumber, &b_EtcNumber);
   fChain->SetBranchAddress("EtcID", EtcID, &b_EtcID);
   fChain->SetBranchAddress("EtcFitHeight", EtcFitHeight, &b_EtcFitHeight);
   fChain->SetBranchAddress("EtcFitTime", EtcFitTime, &b_EtcFitTime);
   fChain->SetBranchAddress("EtcChisq", EtcChisq, &b_EtcChisq);
   fChain->SetBranchAddress("EtcNDF", EtcNDF, &b_EtcNDF);
   fChain->SetBranchAddress("EtcFitShape", EtcFitShape, &b_EtcFitShape);
   fChain->SetBranchAddress("EtcDeltaDiff", EtcDeltaDiff, &b_EtcDeltaDiff);
   fChain->SetBranchAddress("EtcPedestal", EtcPedestal, &b_EtcPedestal);
   fChain->SetBranchAddress("EtcSignal", EtcSignal, &b_EtcSignal);
   fChain->SetBranchAddress("EtcTime", EtcTime, &b_EtcTime);
   fChain->SetBranchAddress("EtcHHTime", EtcHHTime, &b_EtcHHTime);
   fChain->SetBranchAddress("EtcParA", EtcParA, &b_EtcParA);
   fChain->SetBranchAddress("EtcParB", EtcParB, &b_EtcParB);
   fChain->SetBranchAddress("EtcADC", EtcADC, &b_EtcADC);
   fChain->SetBranchAddress("EtcFitADC", EtcFitADC, &b_EtcFitADC);
   Notify();
}

Bool_t E14WavReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void E14WavReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t E14WavReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef E14WavReader_cxx
