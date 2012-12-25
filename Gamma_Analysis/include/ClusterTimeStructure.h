#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include <vector>
#include <list>
#include "cluster/Cluster.h"
#include "gamma/Gamma.h"
#include "klong/Klong.h"


#include "gnana/E14GNAnaDataContainer.h"
#include "gnana/E14GNAnaFunction.h"
#include "rec2g/Rec2g.h"

#include "User_Function.h"

#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TF1.h"

#include "csimap/CsiMap.h"
#include "E14Fsim/E14FsimFunction.h"

class ClusterTime{
 public:
  ClusterTime();
  ~ClusterTime();
  bool operator  < (const ClusterTime& ) const;
  bool operator == (const ClusterTime& ) const;
  int Compare( const ClusterTime& clus) const;

  int id()               const { return m_id;}
  int status()           const { return m_status; }
  double RThreshold()    const { return m_Rthreshold; }

  void  SetClusterMeanTime( double MeanTime ){ m_MeanTime = MeanTime ; }
  void  SetClusterR( double R ){ m_RCluster = R; }
  void  SetClusterPhi( double Phi ){ m_PhiCluster = Phi; }
  void  SetCrystalData( double cID, double cTime, double cR, double cPhi){
    m_cidVec.push_back( cID );
    m_cTimeDeltaVec.push_back( cTime );
    m_cRVec.push_back( cR );
    m_cPhiVec.push_back( cPhi );
  }
  void  Reset(){
    m_id = -1;
    m_status = -1; 
    m_MeanTime = 0.;
    m_RCluster = 0.;
    m_PhiCluster = 0.;
    
    m_cidVec.clear();
    m_cTimeDeltaVec.clear();
    m_cRVec.clear();
    m_cPhiVec.clear();
  }

  double GetClusterTime()   const { return m_MeanTime; }
  double GetClusterR()   const { return m_RCluster; }
  double GetClusterPhi() const { return m_PhiCluster; }
  const std::vector<double>& clusterIDVec()        const { return m_cidVec; }
  const std::vector<double>& clusterTimeDeltaVec() const { return m_cTimeDeltaVec; }
  const std::vector<double>& clusterPhiVec()        const { return m_cPhiVec; }
  const std::vector<double>& clusterRVec()         const { return m_cRVec; }
  

 private:
  int    m_id;
  int    m_status;
  double m_Rthreshold;
  double m_MeanTime;
  double m_RCluster;
  double m_PhiCluster;

  std::vector<double> m_cidVec;
  std::vector<double> m_cTimeDeltaVec;
  std::vector<double> m_cRVec;
  std::vector<double> m_cPhiVec;

};

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

class ClusterTimeAnalyzer {
 public:
  ClusterTimeAnalyzer();
  virtual ~ClusterTimeAnalyzer();  
  const Double_t RadiusThreshold;
  int Convert( std::list<Cluster>& clist, std::list<ClusterTime>& clusterTimeList );
  int ConvertData(Cluster& clus, ClusterTime& clusterTime );
  Double_t GetMeanTime( Cluster clus );
 private:
  Double_t Posx[3000];
  Double_t Posy[3000];
  Double_t Width[3000];
};

