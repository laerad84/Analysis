//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Oct 26 16:06:07 2012 by ROOT version 5.30/04
// from TTree Tree/
// found on file: /Volume0/ExpData/2012_Feb_Beam/RootFile_wav/run_wav_4543.root
//////////////////////////////////////////////////////////

#ifndef E14WavReader_V1_h
#define E14WavReader_V1_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class E14WavReader_V1 {
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
   Int_t           CsinTimeCluster;
   Double_t        CsiTimeClusterHead[24];   //[CsinTimeCluster]
   Double_t        CsiTimeClusterTail[24];   //[CsinTimeCluster]
   Int_t           CsiTotalEnergyInTimeCluster[24];   //[CsinTimeCluster]
   Double_t        CsiTotalEnergy;
   Int_t           CsinOverFlow;
   Int_t           CsiIDOverFlow[2716];   //[CsinOverFlow]
   Int_t           CsinUnderFlow;
   Int_t           CsiIDUnderFlow[2716];   //[CsinUnderFlow]
   Int_t           CsiNumber;
   Short_t         CsiID[2716];   //[CsiNumber]
   Int_t           CsiTimeClusterID[2716];   //[CsiNumber]
   Short_t         CsiFitHeight[2716];   //[CsiNumber]
   Short_t         CsiFitTime[2716];   //[CsiNumber]
   Double_t        CsiChisq[2716];   //[CsiNumber]
   Short_t         CsiNDF[2716];   //[CsiNumber]
   Double_t        CsiFitShape[2716];   //[CsiNumber]
   Double_t        CsiDeltaDiff[2716];   //[CsiNumber]
   Double_t        CsiPedestal[2716];   //[CsiNumber]
   Double_t        CsiSignal[2716];   //[CsiNumber]
   Double_t        CsiEne[2716];   //[CsiNumber]
   Double_t        CsiTime[2716];   //[CsiNumber]
   Double_t        CsiHHTime[2716];   //[CsiNumber]
   Double_t        CsiParA[2716];   //[CsiNumber]
   Double_t        CsiParB[2716];   //[CsiNumber]
   Double_t        CsiADC[2716];   //[CsiNumber]
   Double_t        CsiFitADC[2716];   //[CsiNumber]
   Double_t        Csiwav_SlopeDelta[2716];   //[CsiNumber]
   Double_t        Csiwav_Height[2716];   //[CsiNumber]
   Double_t        Csiwav_PeakTime[2716];   //[CsiNumber]
   Double_t        Csiwav_Pedestal[2716];   //[CsiNumber]
   Double_t        Csiwav_Width[2716];   //[CsiNumber]
   Double_t        Csiwav_FrontHalfTime[2716];   //[CsiNumber]
   Double_t        Csiwav_RearHalfTime[2716];   //[CsiNumber]
   Double_t        CsiFit_Pedestal[2716];   //[CsiNumber]
   Double_t        CsiFit_Time[2716];   //[CsiNumber]
   Double_t        CsiFit_Height[2716];   //[CsiNumber]
   Double_t        CsiFit_HHTime[2716];   //[CsiNumber]
   Double_t        CsiFit_ChisqNDF[2716];   //[CsiNumber]
   Double_t        CsiFit_ChisqPed[2716];   //[CsiNumber]
   Double_t        CsiFit_ChisqFront[2716];   //[CsiNumber]
   Double_t        CsiFit_ChisqRear[2716];   //[CsiNumber]
   Double_t        CsiFit_ChisqTail[2716];   //[CsiNumber]
   Double_t        CsiConv_Energy[2716];   //[CsiNumber]
   Double_t        CsiConv_Time[2716];   //[CsiNumber]
   Int_t           CC03nTimeCluster;
   Double_t        CC03TimeClusterHead[24];   //[CC03nTimeCluster]
   Double_t        CC03TimeClusterTail[24];   //[CC03nTimeCluster]
   Int_t           CC03TotalEnergyInTimeCluster[24];   //[CC03nTimeCluster]
   Double_t        CC03TotalEnergy;
   Int_t           CC03nOverFlow;
   Int_t           CC03IDOverFlow[32];   //[CC03nOverFlow]
   Int_t           CC03nUnderFlow;
   Int_t           CC03IDUnderFlow[32];   //[CC03nUnderFlow]
   Int_t           CC03Number;
   Short_t         CC03ID[32];   //[CC03Number]
   Int_t           CC03TimeClusterID[32];   //[CC03Number]
   Short_t         CC03FitHeight[32];   //[CC03Number]
   Short_t         CC03FitTime[32];   //[CC03Number]
   Double_t        CC03Chisq[32];   //[CC03Number]
   Short_t         CC03NDF[32];   //[CC03Number]
   Double_t        CC03FitShape[32];   //[CC03Number]
   Double_t        CC03DeltaDiff[32];   //[CC03Number]
   Double_t        CC03Pedestal[32];   //[CC03Number]
   Double_t        CC03Signal[32];   //[CC03Number]
   Double_t        CC03Ene[32];   //[CC03Number]
   Double_t        CC03Time[32];   //[CC03Number]
   Double_t        CC03HHTime[32];   //[CC03Number]
   Double_t        CC03ParA[32];   //[CC03Number]
   Double_t        CC03ParB[32];   //[CC03Number]
   Double_t        CC03ADC[32];   //[CC03Number]
   Double_t        CC03FitADC[32];   //[CC03Number]
   Int_t           OEVnTimeCluster;
   Double_t        OEVTimeClusterHead[24];   //[OEVnTimeCluster]
   Double_t        OEVTimeClusterTail[24];   //[OEVnTimeCluster]
   Int_t           OEVTotalEnergyInTimeCluster[24];   //[OEVnTimeCluster]
   Double_t        OEVTotalEnergy;
   Int_t           OEVnOverFlow;
   Int_t           OEVIDOverFlow[44];   //[OEVnOverFlow]
   Int_t           OEVnUnderFlow;
   Int_t           OEVIDUnderFlow[44];   //[OEVnUnderFlow]
   Int_t           OEVNumber;
   Short_t         OEVID[44];   //[OEVNumber]
   Int_t           OEVTimeClusterID[44];   //[OEVNumber]
   Short_t         OEVFitHeight[44];   //[OEVNumber]
   Short_t         OEVFitTime[44];   //[OEVNumber]
   Double_t        OEVChisq[44];   //[OEVNumber]
   Short_t         OEVNDF[44];   //[OEVNumber]
   Double_t        OEVFitShape[44];   //[OEVNumber]
   Double_t        OEVDeltaDiff[44];   //[OEVNumber]
   Double_t        OEVPedestal[44];   //[OEVNumber]
   Double_t        OEVSignal[44];   //[OEVNumber]
   Double_t        OEVEne[44];   //[OEVNumber]
   Double_t        OEVTime[44];   //[OEVNumber]
   Double_t        OEVHHTime[44];   //[OEVNumber]
   Double_t        OEVParA[44];   //[OEVNumber]
   Double_t        OEVParB[44];   //[OEVNumber]
   Double_t        OEVADC[44];   //[OEVNumber]
   Double_t        OEVFitADC[44];   //[OEVNumber]
   Int_t           CVnTimeCluster;
   Double_t        CVTimeClusterHead[24];   //[CVnTimeCluster]
   Double_t        CVTimeClusterTail[24];   //[CVnTimeCluster]
   Int_t           CVTotalEnergyInTimeCluster[24];   //[CVnTimeCluster]
   Double_t        CVTotalEnergy;
   Int_t           CVnOverFlow;
   Int_t           CVIDOverFlow[10];   //[CVnOverFlow]
   Int_t           CVnUnderFlow;
   Int_t           CVIDUnderFlow[10];   //[CVnUnderFlow]
   Int_t           CVNumber;
   Short_t         CVID[10];   //[CVNumber]
   Int_t           CVTimeClusterID[10];   //[CVNumber]
   Short_t         CVFitHeight[10];   //[CVNumber]
   Short_t         CVFitTime[10];   //[CVNumber]
   Double_t        CVChisq[10];   //[CVNumber]
   Short_t         CVNDF[10];   //[CVNumber]
   Double_t        CVFitShape[10];   //[CVNumber]
   Double_t        CVDeltaDiff[10];   //[CVNumber]
   Double_t        CVPedestal[10];   //[CVNumber]
   Double_t        CVSignal[10];   //[CVNumber]
   Double_t        CVEne[10];   //[CVNumber]
   Double_t        CVTime[10];   //[CVNumber]
   Double_t        CVHHTime[10];   //[CVNumber]
   Double_t        CVParA[10];   //[CVNumber]
   Double_t        CVParB[10];   //[CVNumber]
   Double_t        CVADC[10];   //[CVNumber]
   Double_t        CVFitADC[10];   //[CVNumber]
   Int_t           CosmicnTimeCluster;
   Double_t        CosmicTimeClusterHead[24];   //[CosmicnTimeCluster]
   Double_t        CosmicTimeClusterTail[24];   //[CosmicnTimeCluster]
   Int_t           CosmicTotalEnergyInTimeCluster[24];   //[CosmicnTimeCluster]
   Double_t        CosmicTotalEnergy;
   Int_t           CosmicnOverFlow;
   Int_t           CosmicIDOverFlow[20];   //[CosmicnOverFlow]
   Int_t           CosmicnUnderFlow;
   Int_t           CosmicIDUnderFlow[20];   //[CosmicnUnderFlow]
   Int_t           CosmicNumber;
   Short_t         CosmicID[20];   //[CosmicNumber]
   Int_t           CosmicTimeClusterID[20];   //[CosmicNumber]
   Short_t         CosmicFitHeight[20];   //[CosmicNumber]
   Short_t         CosmicFitTime[20];   //[CosmicNumber]
   Double_t        CosmicChisq[20];   //[CosmicNumber]
   Short_t         CosmicNDF[20];   //[CosmicNumber]
   Double_t        CosmicFitShape[20];   //[CosmicNumber]
   Double_t        CosmicDeltaDiff[20];   //[CosmicNumber]
   Double_t        CosmicPedestal[20];   //[CosmicNumber]
   Double_t        CosmicSignal[20];   //[CosmicNumber]
   Double_t        CosmicEne[20];   //[CosmicNumber]
   Double_t        CosmicTime[20];   //[CosmicNumber]
   Double_t        CosmicHHTime[20];   //[CosmicNumber]
   Double_t        CosmicParA[20];   //[CosmicNumber]
   Double_t        CosmicParB[20];   //[CosmicNumber]
   Double_t        CosmicADC[20];   //[CosmicNumber]
   Double_t        CosmicFitADC[20];   //[CosmicNumber]
   Int_t           LasernTimeCluster;
   Double_t        LaserTimeClusterHead[24];   //[LasernTimeCluster]
   Double_t        LaserTimeClusterTail[24];   //[LasernTimeCluster]
   Int_t           LaserTotalEnergyInTimeCluster[24];   //[LasernTimeCluster]
   Double_t        LaserTotalEnergy;
   Int_t           LasernOverFlow;
   Int_t           LaserIDOverFlow[5];   //[LasernOverFlow]
   Int_t           LasernUnderFlow;
   Int_t           LaserIDUnderFlow[5];   //[LasernUnderFlow]
   Int_t           LaserNumber;
   Short_t         LaserID[5];   //[LaserNumber]
   Int_t           LaserTimeClusterID[5];   //[LaserNumber]
   Short_t         LaserFitHeight[5];   //[LaserNumber]
   Short_t         LaserFitTime[5];   //[LaserNumber]
   Double_t        LaserChisq[5];   //[LaserNumber]
   Short_t         LaserNDF[5];   //[LaserNumber]
   Double_t        LaserFitShape[5];   //[LaserNumber]
   Double_t        LaserDeltaDiff[5];   //[LaserNumber]
   Double_t        LaserPedestal[5];   //[LaserNumber]
   Double_t        LaserSignal[5];   //[LaserNumber]
   Double_t        LaserEne[5];   //[LaserNumber]
   Double_t        LaserTime[5];   //[LaserNumber]
   Double_t        LaserHHTime[5];   //[LaserNumber]
   Double_t        LaserParA[5];   //[LaserNumber]
   Double_t        LaserParB[5];   //[LaserNumber]
   Double_t        LaserADC[5];   //[LaserNumber]
   Double_t        LaserFitADC[5];   //[LaserNumber]
   Int_t           EtcnTimeCluster;
   Double_t        EtcTimeClusterHead[24];   //[EtcnTimeCluster]
   Double_t        EtcTimeClusterTail[24];   //[EtcnTimeCluster]
   Int_t           EtcTotalEnergyInTimeCluster[24];   //[EtcnTimeCluster]
   Double_t        EtcTotalEnergy;
   Int_t           EtcnOverFlow;
   Int_t           EtcIDOverFlow[5];   //[EtcnOverFlow]
   Int_t           EtcnUnderFlow;
   Int_t           EtcIDUnderFlow[5];   //[EtcnUnderFlow]
   Int_t           EtcNumber;
   Short_t         EtcID[5];   //[EtcNumber]
   Int_t           EtcTimeClusterID[5];   //[EtcNumber]
   Short_t         EtcFitHeight[5];   //[EtcNumber]
   Short_t         EtcFitTime[5];   //[EtcNumber]
   Double_t        EtcChisq[5];   //[EtcNumber]
   Short_t         EtcNDF[5];   //[EtcNumber]
   Double_t        EtcFitShape[5];   //[EtcNumber]
   Double_t        EtcDeltaDiff[5];   //[EtcNumber]
   Double_t        EtcPedestal[5];   //[EtcNumber]
   Double_t        EtcSignal[5];   //[EtcNumber]
   Double_t        EtcEne[5];   //[EtcNumber]
   Double_t        EtcTime[5];   //[EtcNumber]
   Double_t        EtcHHTime[5];   //[EtcNumber]
   Double_t        EtcParA[5];   //[EtcNumber]
   Double_t        EtcParB[5];   //[EtcNumber]
   Double_t        EtcADC[5];   //[EtcNumber]
   Double_t        EtcFitADC[5];   //[EtcNumber]

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
   TBranch        *b_CsinTimeCluster;   //!
   TBranch        *b_CsiTimeClusterHead;   //!
   TBranch        *b_CsiTimeClusterTail;   //!
   TBranch        *b_CsiTotalEnergyInTimeCluster;   //!
   TBranch        *b_CsiTotalEnergy;   //!
   TBranch        *b_CsinOverFlow;   //!
   TBranch        *b_CsiIDOverFlow;   //!
   TBranch        *b_CsinUnderFlow;   //!
   TBranch        *b_CsiIDUnderFlow;   //!
   TBranch        *b_CsiNumber;   //!
   TBranch        *b_CsiID;   //!
   TBranch        *b_CsiTimeClusterID;   //!
   TBranch        *b_CsiFitHeight;   //!
   TBranch        *b_CsiFitTime;   //!
   TBranch        *b_CsiChisq;   //!
   TBranch        *b_CsiNDF;   //!
   TBranch        *b_CsiFitShape;   //!
   TBranch        *b_CsiDeltaDiff;   //!
   TBranch        *b_CsiPedestal;   //!
   TBranch        *b_CsiSignal;   //!
   TBranch        *b_CsiEne;   //!
   TBranch        *b_CsiTime;   //!
   TBranch        *b_CsiHHTime;   //!
   TBranch        *b_CsiParA;   //!
   TBranch        *b_CsiParB;   //!
   TBranch        *b_CsiADC;   //!
   TBranch        *b_CsiFitADC;   //!
   TBranch        *b_Csiwav_SlopeDelta;   //!
   TBranch        *b_Csiwav_Height;   //!
   TBranch        *b_Csiwav_PeakTime;   //!
   TBranch        *b_Csiwav_Pedestal;   //!
   TBranch        *b_Csiwav_Width;   //!
   TBranch        *b_Csiwav_FrontHalfTime;   //!
   TBranch        *b_Csiwav_RearHalfTime;   //!
   TBranch        *b_CsiFit_Pedestal;   //!
   TBranch        *b_CsiFit_Time;   //!
   TBranch        *b_CsiFit_Height;   //!
   TBranch        *b_CsiFit_HHTime;   //!
   TBranch        *b_CsiFit_ChisqNDF;   //!
   TBranch        *b_CsiFit_ChisqPed;   //!
   TBranch        *b_CsiFit_ChisqFront;   //!
   TBranch        *b_CsiFit_ChisqRear;   //!
   TBranch        *b_CsiFit_ChisqTail;   //!
   TBranch        *b_CsiConv_Energy;   //!
   TBranch        *b_CsiConv_Time;   //!
   TBranch        *b_CC03nTimeCluster;   //!
   TBranch        *b_CC03TimeClusterHead;   //!
   TBranch        *b_CC03TimeClusterTail;   //!
   TBranch        *b_CC03TotalEnergyInTimeCluster;   //!
   TBranch        *b_CC03TotalEnergy;   //!
   TBranch        *b_CC03nOverFlow;   //!
   TBranch        *b_CC03IDOverFlow;   //!
   TBranch        *b_CC03nUnderFlow;   //!
   TBranch        *b_CC03IDUnderFlow;   //!
   TBranch        *b_CC03Number;   //!
   TBranch        *b_CC03ID;   //!
   TBranch        *b_CC03TimeClusterID;   //!
   TBranch        *b_CC03FitHeight;   //!
   TBranch        *b_CC03FitTime;   //!
   TBranch        *b_CC03Chisq;   //!
   TBranch        *b_CC03NDF;   //!
   TBranch        *b_CC03FitShape;   //!
   TBranch        *b_CC03DeltaDiff;   //!
   TBranch        *b_CC03Pedestal;   //!
   TBranch        *b_CC03Signal;   //!
   TBranch        *b_CC03Ene;   //!
   TBranch        *b_CC03Time;   //!
   TBranch        *b_CC03HHTime;   //!
   TBranch        *b_CC03ParA;   //!
   TBranch        *b_CC03ParB;   //!
   TBranch        *b_CC03ADC;   //!
   TBranch        *b_CC03FitADC;   //!
   TBranch        *b_OEVnTimeCluster;   //!
   TBranch        *b_OEVTimeClusterHead;   //!
   TBranch        *b_OEVTimeClusterTail;   //!
   TBranch        *b_OEVTotalEnergyInTimeCluster;   //!
   TBranch        *b_OEVTotalEnergy;   //!
   TBranch        *b_OEVnOverFlow;   //!
   TBranch        *b_OEVIDOverFlow;   //!
   TBranch        *b_OEVnUnderFlow;   //!
   TBranch        *b_OEVIDUnderFlow;   //!
   TBranch        *b_OEVNumber;   //!
   TBranch        *b_OEVID;   //!
   TBranch        *b_OEVTimeClusterID;   //!
   TBranch        *b_OEVFitHeight;   //!
   TBranch        *b_OEVFitTime;   //!
   TBranch        *b_OEVChisq;   //!
   TBranch        *b_OEVNDF;   //!
   TBranch        *b_OEVFitShape;   //!
   TBranch        *b_OEVDeltaDiff;   //!
   TBranch        *b_OEVPedestal;   //!
   TBranch        *b_OEVSignal;   //!
   TBranch        *b_OEVEne;   //!
   TBranch        *b_OEVTime;   //!
   TBranch        *b_OEVHHTime;   //!
   TBranch        *b_OEVParA;   //!
   TBranch        *b_OEVParB;   //!
   TBranch        *b_OEVADC;   //!
   TBranch        *b_OEVFitADC;   //!
   TBranch        *b_CVnTimeCluster;   //!
   TBranch        *b_CVTimeClusterHead;   //!
   TBranch        *b_CVTimeClusterTail;   //!
   TBranch        *b_CVTotalEnergyInTimeCluster;   //!
   TBranch        *b_CVTotalEnergy;   //!
   TBranch        *b_CVnOverFlow;   //!
   TBranch        *b_CVIDOverFlow;   //!
   TBranch        *b_CVnUnderFlow;   //!
   TBranch        *b_CVIDUnderFlow;   //!
   TBranch        *b_CVNumber;   //!
   TBranch        *b_CVID;   //!
   TBranch        *b_CVTimeClusterID;   //!
   TBranch        *b_CVFitHeight;   //!
   TBranch        *b_CVFitTime;   //!
   TBranch        *b_CVChisq;   //!
   TBranch        *b_CVNDF;   //!
   TBranch        *b_CVFitShape;   //!
   TBranch        *b_CVDeltaDiff;   //!
   TBranch        *b_CVPedestal;   //!
   TBranch        *b_CVSignal;   //!
   TBranch        *b_CVEne;   //!
   TBranch        *b_CVTime;   //!
   TBranch        *b_CVHHTime;   //!
   TBranch        *b_CVParA;   //!
   TBranch        *b_CVParB;   //!
   TBranch        *b_CVADC;   //!
   TBranch        *b_CVFitADC;   //!
   TBranch        *b_CosmicnTimeCluster;   //!
   TBranch        *b_CosmicTimeClusterHead;   //!
   TBranch        *b_CosmicTimeClusterTail;   //!
   TBranch        *b_CosmicTotalEnergyInTimeCluster;   //!
   TBranch        *b_CosmicTotalEnergy;   //!
   TBranch        *b_CosmicnOverFlow;   //!
   TBranch        *b_CosmicIDOverFlow;   //!
   TBranch        *b_CosmicnUnderFlow;   //!
   TBranch        *b_CosmicIDUnderFlow;   //!
   TBranch        *b_CosmicNumber;   //!
   TBranch        *b_CosmicID;   //!
   TBranch        *b_CosmicTimeClusterID;   //!
   TBranch        *b_CosmicFitHeight;   //!
   TBranch        *b_CosmicFitTime;   //!
   TBranch        *b_CosmicChisq;   //!
   TBranch        *b_CosmicNDF;   //!
   TBranch        *b_CosmicFitShape;   //!
   TBranch        *b_CosmicDeltaDiff;   //!
   TBranch        *b_CosmicPedestal;   //!
   TBranch        *b_CosmicSignal;   //!
   TBranch        *b_CosmicEne;   //!
   TBranch        *b_CosmicTime;   //!
   TBranch        *b_CosmicHHTime;   //!
   TBranch        *b_CosmicParA;   //!
   TBranch        *b_CosmicParB;   //!
   TBranch        *b_CosmicADC;   //!
   TBranch        *b_CosmicFitADC;   //!
   TBranch        *b_LasernTimeCluster;   //!
   TBranch        *b_LaserTimeClusterHead;   //!
   TBranch        *b_LaserTimeClusterTail;   //!
   TBranch        *b_LaserTotalEnergyInTimeCluster;   //!
   TBranch        *b_LaserTotalEnergy;   //!
   TBranch        *b_LasernOverFlow;   //!
   TBranch        *b_LaserIDOverFlow;   //!
   TBranch        *b_LasernUnderFlow;   //!
   TBranch        *b_LaserIDUnderFlow;   //!
   TBranch        *b_LaserNumber;   //!
   TBranch        *b_LaserID;   //!
   TBranch        *b_LaserTimeClusterID;   //!
   TBranch        *b_LaserFitHeight;   //!
   TBranch        *b_LaserFitTime;   //!
   TBranch        *b_LaserChisq;   //!
   TBranch        *b_LaserNDF;   //!
   TBranch        *b_LaserFitShape;   //!
   TBranch        *b_LaserDeltaDiff;   //!
   TBranch        *b_LaserPedestal;   //!
   TBranch        *b_LaserSignal;   //!
   TBranch        *b_LaserEne;   //!
   TBranch        *b_LaserTime;   //!
   TBranch        *b_LaserHHTime;   //!
   TBranch        *b_LaserParA;   //!
   TBranch        *b_LaserParB;   //!
   TBranch        *b_LaserADC;   //!
   TBranch        *b_LaserFitADC;   //!
   TBranch        *b_EtcnTimeCluster;   //!
   TBranch        *b_EtcTimeClusterHead;   //!
   TBranch        *b_EtcTimeClusterTail;   //!
   TBranch        *b_EtcTotalEnergyInTimeCluster;   //!
   TBranch        *b_EtcTotalEnergy;   //!
   TBranch        *b_EtcnOverFlow;   //!
   TBranch        *b_EtcIDOverFlow;   //!
   TBranch        *b_EtcnUnderFlow;   //!
   TBranch        *b_EtcIDUnderFlow;   //!
   TBranch        *b_EtcNumber;   //!
   TBranch        *b_EtcID;   //!
   TBranch        *b_EtcTimeClusterID;   //!
   TBranch        *b_EtcFitHeight;   //!
   TBranch        *b_EtcFitTime;   //!
   TBranch        *b_EtcChisq;   //!
   TBranch        *b_EtcNDF;   //!
   TBranch        *b_EtcFitShape;   //!
   TBranch        *b_EtcDeltaDiff;   //!
   TBranch        *b_EtcPedestal;   //!
   TBranch        *b_EtcSignal;   //!
   TBranch        *b_EtcEne;   //!
   TBranch        *b_EtcTime;   //!
   TBranch        *b_EtcHHTime;   //!
   TBranch        *b_EtcParA;   //!
   TBranch        *b_EtcParB;   //!
   TBranch        *b_EtcADC;   //!
   TBranch        *b_EtcFitADC;   //!

   E14WavReader_V1(TTree *tree=0);
   virtual ~E14WavReader_V1();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef E14WavReader_V1_cxx
E14WavReader_V1::E14WavReader_V1(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   Init(tree);
}

E14WavReader_V1::~E14WavReader_V1()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t E14WavReader_V1::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t E14WavReader_V1::LoadTree(Long64_t entry)
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

void E14WavReader_V1::Init(TTree *tree)
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
   fChain->SetBranchAddress("CsinTimeCluster", &CsinTimeCluster, &b_CsinTimeCluster);
   fChain->SetBranchAddress("CsiTimeClusterHead", CsiTimeClusterHead, &b_CsiTimeClusterHead);
   fChain->SetBranchAddress("CsiTimeClusterTail", CsiTimeClusterTail, &b_CsiTimeClusterTail);
   fChain->SetBranchAddress("CsiTotalEnergyInTimeCluster", CsiTotalEnergyInTimeCluster, &b_CsiTotalEnergyInTimeCluster);
   fChain->SetBranchAddress("CsiTotalEnergy", &CsiTotalEnergy, &b_CsiTotalEnergy);
   fChain->SetBranchAddress("CsinOverFlow", &CsinOverFlow, &b_CsinOverFlow);
   fChain->SetBranchAddress("CsiIDOverFlow", &CsiIDOverFlow, &b_CsiIDOverFlow);
   fChain->SetBranchAddress("CsinUnderFlow", &CsinUnderFlow, &b_CsinUnderFlow);
   fChain->SetBranchAddress("CsiIDUnderFlow", &CsiIDUnderFlow, &b_CsiIDUnderFlow);
   fChain->SetBranchAddress("CsiNumber", &CsiNumber, &b_CsiNumber);
   fChain->SetBranchAddress("CsiID", CsiID, &b_CsiID);
   fChain->SetBranchAddress("CsiTimeClusterID", CsiTimeClusterID, &b_CsiTimeClusterID);
   fChain->SetBranchAddress("CsiFitHeight", CsiFitHeight, &b_CsiFitHeight);
   fChain->SetBranchAddress("CsiFitTime", CsiFitTime, &b_CsiFitTime);
   fChain->SetBranchAddress("CsiChisq", CsiChisq, &b_CsiChisq);
   fChain->SetBranchAddress("CsiNDF", CsiNDF, &b_CsiNDF);
   fChain->SetBranchAddress("CsiFitShape", CsiFitShape, &b_CsiFitShape);
   fChain->SetBranchAddress("CsiDeltaDiff", CsiDeltaDiff, &b_CsiDeltaDiff);
   fChain->SetBranchAddress("CsiPedestal", CsiPedestal, &b_CsiPedestal);
   fChain->SetBranchAddress("CsiSignal", CsiSignal, &b_CsiSignal);
   fChain->SetBranchAddress("CsiEne", CsiEne, &b_CsiEne);
   fChain->SetBranchAddress("CsiTime", CsiTime, &b_CsiTime);
   fChain->SetBranchAddress("CsiHHTime", CsiHHTime, &b_CsiHHTime);
   fChain->SetBranchAddress("CsiParA", CsiParA, &b_CsiParA);
   fChain->SetBranchAddress("CsiParB", CsiParB, &b_CsiParB);
   fChain->SetBranchAddress("CsiADC", CsiADC, &b_CsiADC);
   fChain->SetBranchAddress("CsiFitADC", CsiFitADC, &b_CsiFitADC);
   fChain->SetBranchAddress("Csiwav_SlopeDelta", Csiwav_SlopeDelta, &b_Csiwav_SlopeDelta);
   fChain->SetBranchAddress("Csiwav_Height", Csiwav_Height, &b_Csiwav_Height);
   fChain->SetBranchAddress("Csiwav_PeakTime", Csiwav_PeakTime, &b_Csiwav_PeakTime);
   fChain->SetBranchAddress("Csiwav_Pedestal", Csiwav_Pedestal, &b_Csiwav_Pedestal);
   fChain->SetBranchAddress("Csiwav_Width", Csiwav_Width, &b_Csiwav_Width);
   fChain->SetBranchAddress("Csiwav_FrontHalfTime", Csiwav_FrontHalfTime, &b_Csiwav_FrontHalfTime);
   fChain->SetBranchAddress("Csiwav_RearHalfTime", Csiwav_RearHalfTime, &b_Csiwav_RearHalfTime);
   fChain->SetBranchAddress("CsiFit_Pedestal", CsiFit_Pedestal, &b_CsiFit_Pedestal);
   fChain->SetBranchAddress("CsiFit_Time", CsiFit_Time, &b_CsiFit_Time);
   fChain->SetBranchAddress("CsiFit_Height", CsiFit_Height, &b_CsiFit_Height);
   fChain->SetBranchAddress("CsiFit_HHTime", CsiFit_HHTime, &b_CsiFit_HHTime);
   fChain->SetBranchAddress("CsiFit_ChisqNDF", CsiFit_ChisqNDF, &b_CsiFit_ChisqNDF);
   fChain->SetBranchAddress("CsiFit_ChisqPed", CsiFit_ChisqPed, &b_CsiFit_ChisqPed);
   fChain->SetBranchAddress("CsiFit_ChisqFront", CsiFit_ChisqFront, &b_CsiFit_ChisqFront);
   fChain->SetBranchAddress("CsiFit_ChisqRear", CsiFit_ChisqRear, &b_CsiFit_ChisqRear);
   fChain->SetBranchAddress("CsiFit_ChisqTail", CsiFit_ChisqTail, &b_CsiFit_ChisqTail);
   fChain->SetBranchAddress("CsiConv_Energy", CsiConv_Energy, &b_CsiConv_Energy);
   fChain->SetBranchAddress("CsiConv_Time", CsiConv_Time, &b_CsiConv_Time);
   fChain->SetBranchAddress("CC03nTimeCluster", &CC03nTimeCluster, &b_CC03nTimeCluster);
   fChain->SetBranchAddress("CC03TimeClusterHead", &CC03TimeClusterHead, &b_CC03TimeClusterHead);
   fChain->SetBranchAddress("CC03TimeClusterTail", &CC03TimeClusterTail, &b_CC03TimeClusterTail);
   fChain->SetBranchAddress("CC03TotalEnergyInTimeCluster", &CC03TotalEnergyInTimeCluster, &b_CC03TotalEnergyInTimeCluster);
   fChain->SetBranchAddress("CC03TotalEnergy", &CC03TotalEnergy, &b_CC03TotalEnergy);
   fChain->SetBranchAddress("CC03nOverFlow", &CC03nOverFlow, &b_CC03nOverFlow);
   fChain->SetBranchAddress("CC03IDOverFlow", &CC03IDOverFlow, &b_CC03IDOverFlow);
   fChain->SetBranchAddress("CC03nUnderFlow", &CC03nUnderFlow, &b_CC03nUnderFlow);
   fChain->SetBranchAddress("CC03IDUnderFlow", &CC03IDUnderFlow, &b_CC03IDUnderFlow);
   fChain->SetBranchAddress("CC03Number", &CC03Number, &b_CC03Number);
   fChain->SetBranchAddress("CC03ID", CC03ID, &b_CC03ID);
   fChain->SetBranchAddress("CC03TimeClusterID", CC03TimeClusterID, &b_CC03TimeClusterID);
   fChain->SetBranchAddress("CC03FitHeight", CC03FitHeight, &b_CC03FitHeight);
   fChain->SetBranchAddress("CC03FitTime", CC03FitTime, &b_CC03FitTime);
   fChain->SetBranchAddress("CC03Chisq", CC03Chisq, &b_CC03Chisq);
   fChain->SetBranchAddress("CC03NDF", CC03NDF, &b_CC03NDF);
   fChain->SetBranchAddress("CC03FitShape", CC03FitShape, &b_CC03FitShape);
   fChain->SetBranchAddress("CC03DeltaDiff", CC03DeltaDiff, &b_CC03DeltaDiff);
   fChain->SetBranchAddress("CC03Pedestal", CC03Pedestal, &b_CC03Pedestal);
   fChain->SetBranchAddress("CC03Signal", CC03Signal, &b_CC03Signal);
   fChain->SetBranchAddress("CC03Ene", CC03Ene, &b_CC03Ene);
   fChain->SetBranchAddress("CC03Time", CC03Time, &b_CC03Time);
   fChain->SetBranchAddress("CC03HHTime", CC03HHTime, &b_CC03HHTime);
   fChain->SetBranchAddress("CC03ParA", CC03ParA, &b_CC03ParA);
   fChain->SetBranchAddress("CC03ParB", CC03ParB, &b_CC03ParB);
   fChain->SetBranchAddress("CC03ADC", CC03ADC, &b_CC03ADC);
   fChain->SetBranchAddress("CC03FitADC", CC03FitADC, &b_CC03FitADC);
   fChain->SetBranchAddress("OEVnTimeCluster", &OEVnTimeCluster, &b_OEVnTimeCluster);
   fChain->SetBranchAddress("OEVTimeClusterHead", &OEVTimeClusterHead, &b_OEVTimeClusterHead);
   fChain->SetBranchAddress("OEVTimeClusterTail", &OEVTimeClusterTail, &b_OEVTimeClusterTail);
   fChain->SetBranchAddress("OEVTotalEnergyInTimeCluster", &OEVTotalEnergyInTimeCluster, &b_OEVTotalEnergyInTimeCluster);
   fChain->SetBranchAddress("OEVTotalEnergy", &OEVTotalEnergy, &b_OEVTotalEnergy);
   fChain->SetBranchAddress("OEVnOverFlow", &OEVnOverFlow, &b_OEVnOverFlow);
   fChain->SetBranchAddress("OEVIDOverFlow", &OEVIDOverFlow, &b_OEVIDOverFlow);
   fChain->SetBranchAddress("OEVnUnderFlow", &OEVnUnderFlow, &b_OEVnUnderFlow);
   fChain->SetBranchAddress("OEVIDUnderFlow", &OEVIDUnderFlow, &b_OEVIDUnderFlow);
   fChain->SetBranchAddress("OEVNumber", &OEVNumber, &b_OEVNumber);
   fChain->SetBranchAddress("OEVID", OEVID, &b_OEVID);
   fChain->SetBranchAddress("OEVTimeClusterID", OEVTimeClusterID, &b_OEVTimeClusterID);
   fChain->SetBranchAddress("OEVFitHeight", OEVFitHeight, &b_OEVFitHeight);
   fChain->SetBranchAddress("OEVFitTime", OEVFitTime, &b_OEVFitTime);
   fChain->SetBranchAddress("OEVChisq", OEVChisq, &b_OEVChisq);
   fChain->SetBranchAddress("OEVNDF", OEVNDF, &b_OEVNDF);
   fChain->SetBranchAddress("OEVFitShape", OEVFitShape, &b_OEVFitShape);
   fChain->SetBranchAddress("OEVDeltaDiff", OEVDeltaDiff, &b_OEVDeltaDiff);
   fChain->SetBranchAddress("OEVPedestal", OEVPedestal, &b_OEVPedestal);
   fChain->SetBranchAddress("OEVSignal", OEVSignal, &b_OEVSignal);
   fChain->SetBranchAddress("OEVEne", OEVEne, &b_OEVEne);
   fChain->SetBranchAddress("OEVTime", OEVTime, &b_OEVTime);
   fChain->SetBranchAddress("OEVHHTime", OEVHHTime, &b_OEVHHTime);
   fChain->SetBranchAddress("OEVParA", OEVParA, &b_OEVParA);
   fChain->SetBranchAddress("OEVParB", OEVParB, &b_OEVParB);
   fChain->SetBranchAddress("OEVADC", OEVADC, &b_OEVADC);
   fChain->SetBranchAddress("OEVFitADC", OEVFitADC, &b_OEVFitADC);
   fChain->SetBranchAddress("CVnTimeCluster", &CVnTimeCluster, &b_CVnTimeCluster);
   fChain->SetBranchAddress("CVTimeClusterHead", &CVTimeClusterHead, &b_CVTimeClusterHead);
   fChain->SetBranchAddress("CVTimeClusterTail", &CVTimeClusterTail, &b_CVTimeClusterTail);
   fChain->SetBranchAddress("CVTotalEnergyInTimeCluster", &CVTotalEnergyInTimeCluster, &b_CVTotalEnergyInTimeCluster);
   fChain->SetBranchAddress("CVTotalEnergy", &CVTotalEnergy, &b_CVTotalEnergy);
   fChain->SetBranchAddress("CVnOverFlow", &CVnOverFlow, &b_CVnOverFlow);
   fChain->SetBranchAddress("CVIDOverFlow", &CVIDOverFlow, &b_CVIDOverFlow);
   fChain->SetBranchAddress("CVnUnderFlow", &CVnUnderFlow, &b_CVnUnderFlow);
   fChain->SetBranchAddress("CVIDUnderFlow", &CVIDUnderFlow, &b_CVIDUnderFlow);
   fChain->SetBranchAddress("CVNumber", &CVNumber, &b_CVNumber);
   fChain->SetBranchAddress("CVID", CVID, &b_CVID);
   fChain->SetBranchAddress("CVTimeClusterID", CVTimeClusterID, &b_CVTimeClusterID);
   fChain->SetBranchAddress("CVFitHeight", CVFitHeight, &b_CVFitHeight);
   fChain->SetBranchAddress("CVFitTime", CVFitTime, &b_CVFitTime);
   fChain->SetBranchAddress("CVChisq", CVChisq, &b_CVChisq);
   fChain->SetBranchAddress("CVNDF", CVNDF, &b_CVNDF);
   fChain->SetBranchAddress("CVFitShape", CVFitShape, &b_CVFitShape);
   fChain->SetBranchAddress("CVDeltaDiff", CVDeltaDiff, &b_CVDeltaDiff);
   fChain->SetBranchAddress("CVPedestal", CVPedestal, &b_CVPedestal);
   fChain->SetBranchAddress("CVSignal", CVSignal, &b_CVSignal);
   fChain->SetBranchAddress("CVEne", CVEne, &b_CVEne);
   fChain->SetBranchAddress("CVTime", CVTime, &b_CVTime);
   fChain->SetBranchAddress("CVHHTime", CVHHTime, &b_CVHHTime);
   fChain->SetBranchAddress("CVParA", CVParA, &b_CVParA);
   fChain->SetBranchAddress("CVParB", CVParB, &b_CVParB);
   fChain->SetBranchAddress("CVADC", CVADC, &b_CVADC);
   fChain->SetBranchAddress("CVFitADC", CVFitADC, &b_CVFitADC);
   fChain->SetBranchAddress("CosmicnTimeCluster", &CosmicnTimeCluster, &b_CosmicnTimeCluster);
   fChain->SetBranchAddress("CosmicTimeClusterHead", &CosmicTimeClusterHead, &b_CosmicTimeClusterHead);
   fChain->SetBranchAddress("CosmicTimeClusterTail", &CosmicTimeClusterTail, &b_CosmicTimeClusterTail);
   fChain->SetBranchAddress("CosmicTotalEnergyInTimeCluster", &CosmicTotalEnergyInTimeCluster, &b_CosmicTotalEnergyInTimeCluster);
   fChain->SetBranchAddress("CosmicTotalEnergy", &CosmicTotalEnergy, &b_CosmicTotalEnergy);
   fChain->SetBranchAddress("CosmicnOverFlow", &CosmicnOverFlow, &b_CosmicnOverFlow);
   fChain->SetBranchAddress("CosmicIDOverFlow", &CosmicIDOverFlow, &b_CosmicIDOverFlow);
   fChain->SetBranchAddress("CosmicnUnderFlow", &CosmicnUnderFlow, &b_CosmicnUnderFlow);
   fChain->SetBranchAddress("CosmicIDUnderFlow", &CosmicIDUnderFlow, &b_CosmicIDUnderFlow);
   fChain->SetBranchAddress("CosmicNumber", &CosmicNumber, &b_CosmicNumber);
   fChain->SetBranchAddress("CosmicID", CosmicID, &b_CosmicID);
   fChain->SetBranchAddress("CosmicTimeClusterID", CosmicTimeClusterID, &b_CosmicTimeClusterID);
   fChain->SetBranchAddress("CosmicFitHeight", CosmicFitHeight, &b_CosmicFitHeight);
   fChain->SetBranchAddress("CosmicFitTime", CosmicFitTime, &b_CosmicFitTime);
   fChain->SetBranchAddress("CosmicChisq", CosmicChisq, &b_CosmicChisq);
   fChain->SetBranchAddress("CosmicNDF", CosmicNDF, &b_CosmicNDF);
   fChain->SetBranchAddress("CosmicFitShape", CosmicFitShape, &b_CosmicFitShape);
   fChain->SetBranchAddress("CosmicDeltaDiff", CosmicDeltaDiff, &b_CosmicDeltaDiff);
   fChain->SetBranchAddress("CosmicPedestal", CosmicPedestal, &b_CosmicPedestal);
   fChain->SetBranchAddress("CosmicSignal", CosmicSignal, &b_CosmicSignal);
   fChain->SetBranchAddress("CosmicEne", CosmicEne, &b_CosmicEne);
   fChain->SetBranchAddress("CosmicTime", CosmicTime, &b_CosmicTime);
   fChain->SetBranchAddress("CosmicHHTime", CosmicHHTime, &b_CosmicHHTime);
   fChain->SetBranchAddress("CosmicParA", CosmicParA, &b_CosmicParA);
   fChain->SetBranchAddress("CosmicParB", CosmicParB, &b_CosmicParB);
   fChain->SetBranchAddress("CosmicADC", CosmicADC, &b_CosmicADC);
   fChain->SetBranchAddress("CosmicFitADC", CosmicFitADC, &b_CosmicFitADC);
   fChain->SetBranchAddress("LasernTimeCluster", &LasernTimeCluster, &b_LasernTimeCluster);
   fChain->SetBranchAddress("LaserTimeClusterHead", &LaserTimeClusterHead, &b_LaserTimeClusterHead);
   fChain->SetBranchAddress("LaserTimeClusterTail", &LaserTimeClusterTail, &b_LaserTimeClusterTail);
   fChain->SetBranchAddress("LaserTotalEnergyInTimeCluster", &LaserTotalEnergyInTimeCluster, &b_LaserTotalEnergyInTimeCluster);
   fChain->SetBranchAddress("LaserTotalEnergy", &LaserTotalEnergy, &b_LaserTotalEnergy);
   fChain->SetBranchAddress("LasernOverFlow", &LasernOverFlow, &b_LasernOverFlow);
   fChain->SetBranchAddress("LaserIDOverFlow", &LaserIDOverFlow, &b_LaserIDOverFlow);
   fChain->SetBranchAddress("LasernUnderFlow", &LasernUnderFlow, &b_LasernUnderFlow);
   fChain->SetBranchAddress("LaserIDUnderFlow", &LaserIDUnderFlow, &b_LaserIDUnderFlow);
   fChain->SetBranchAddress("LaserNumber", &LaserNumber, &b_LaserNumber);
   fChain->SetBranchAddress("LaserID", LaserID, &b_LaserID);
   fChain->SetBranchAddress("LaserTimeClusterID", LaserTimeClusterID, &b_LaserTimeClusterID);
   fChain->SetBranchAddress("LaserFitHeight", LaserFitHeight, &b_LaserFitHeight);
   fChain->SetBranchAddress("LaserFitTime", LaserFitTime, &b_LaserFitTime);
   fChain->SetBranchAddress("LaserChisq", LaserChisq, &b_LaserChisq);
   fChain->SetBranchAddress("LaserNDF", LaserNDF, &b_LaserNDF);
   fChain->SetBranchAddress("LaserFitShape", LaserFitShape, &b_LaserFitShape);
   fChain->SetBranchAddress("LaserDeltaDiff", LaserDeltaDiff, &b_LaserDeltaDiff);
   fChain->SetBranchAddress("LaserPedestal", LaserPedestal, &b_LaserPedestal);
   fChain->SetBranchAddress("LaserSignal", LaserSignal, &b_LaserSignal);
   fChain->SetBranchAddress("LaserEne", LaserEne, &b_LaserEne);
   fChain->SetBranchAddress("LaserTime", LaserTime, &b_LaserTime);
   fChain->SetBranchAddress("LaserHHTime", LaserHHTime, &b_LaserHHTime);
   fChain->SetBranchAddress("LaserParA", LaserParA, &b_LaserParA);
   fChain->SetBranchAddress("LaserParB", LaserParB, &b_LaserParB);
   fChain->SetBranchAddress("LaserADC", LaserADC, &b_LaserADC);
   fChain->SetBranchAddress("LaserFitADC", LaserFitADC, &b_LaserFitADC);
   fChain->SetBranchAddress("EtcnTimeCluster", &EtcnTimeCluster, &b_EtcnTimeCluster);
   fChain->SetBranchAddress("EtcTimeClusterHead", &EtcTimeClusterHead, &b_EtcTimeClusterHead);
   fChain->SetBranchAddress("EtcTimeClusterTail", &EtcTimeClusterTail, &b_EtcTimeClusterTail);
   fChain->SetBranchAddress("EtcTotalEnergyInTimeCluster", &EtcTotalEnergyInTimeCluster, &b_EtcTotalEnergyInTimeCluster);
   fChain->SetBranchAddress("EtcTotalEnergy", &EtcTotalEnergy, &b_EtcTotalEnergy);
   fChain->SetBranchAddress("EtcnOverFlow", &EtcnOverFlow, &b_EtcnOverFlow);
   fChain->SetBranchAddress("EtcIDOverFlow", &EtcIDOverFlow, &b_EtcIDOverFlow);
   fChain->SetBranchAddress("EtcnUnderFlow", &EtcnUnderFlow, &b_EtcnUnderFlow);
   fChain->SetBranchAddress("EtcIDUnderFlow", &EtcIDUnderFlow, &b_EtcIDUnderFlow);
   fChain->SetBranchAddress("EtcNumber", &EtcNumber, &b_EtcNumber);
   fChain->SetBranchAddress("EtcID", EtcID, &b_EtcID);
   fChain->SetBranchAddress("EtcTimeClusterID", EtcTimeClusterID, &b_EtcTimeClusterID);
   fChain->SetBranchAddress("EtcFitHeight", EtcFitHeight, &b_EtcFitHeight);
   fChain->SetBranchAddress("EtcFitTime", EtcFitTime, &b_EtcFitTime);
   fChain->SetBranchAddress("EtcChisq", EtcChisq, &b_EtcChisq);
   fChain->SetBranchAddress("EtcNDF", EtcNDF, &b_EtcNDF);
   fChain->SetBranchAddress("EtcFitShape", EtcFitShape, &b_EtcFitShape);
   fChain->SetBranchAddress("EtcDeltaDiff", EtcDeltaDiff, &b_EtcDeltaDiff);
   fChain->SetBranchAddress("EtcPedestal", EtcPedestal, &b_EtcPedestal);
   fChain->SetBranchAddress("EtcSignal", EtcSignal, &b_EtcSignal);
   fChain->SetBranchAddress("EtcEne", EtcEne, &b_EtcEne);
   fChain->SetBranchAddress("EtcTime", EtcTime, &b_EtcTime);
   fChain->SetBranchAddress("EtcHHTime", EtcHHTime, &b_EtcHHTime);
   fChain->SetBranchAddress("EtcParA", EtcParA, &b_EtcParA);
   fChain->SetBranchAddress("EtcParB", EtcParB, &b_EtcParB);
   fChain->SetBranchAddress("EtcADC", EtcADC, &b_EtcADC);
   fChain->SetBranchAddress("EtcFitADC", EtcFitADC, &b_EtcFitADC);
   Notify();
}

Bool_t E14WavReader_V1::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void E14WavReader_V1::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t E14WavReader_V1::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef E14WavReader_V1_cxx
