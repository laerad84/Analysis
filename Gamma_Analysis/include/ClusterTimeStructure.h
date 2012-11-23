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

#include "csimap/CsiMap.h"
#include "E14Fsim/E14FsimFunction.h"



class ClusterTime{
 public:
  clusterTime();
  ~ClusterTime();
  bool operator  < (const ClusterTime& ) const;
  bool operator == (const ClusterTime& ) const;
  std::ostream& operator << ( std::ostream& out , const Cluster& clus );
  int Compare( const ClusterTime& clus) const;

  int id()               const { return m_id;}
  int status()           const { return m_status; }
  double RThreshold()    const { return m_Rthreshold; }

  void  SetClusterMeanTime( double MeanTime ){ m_MeanTime = MenTime ; }
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
  const std::vector<int>&    clusterIDVec()        const { return m_cidVec; }
  const std::vector<double>& clusterTimeDeltaVec() const { return m_cTimeDeltaVec; }
  const std::vector<double>& clusterPiVec()        const { return m_cPhiVec; }
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

ClusterTime::ClusterTime() : m_id(-1),
  m_status = -1,
  m_Rthreshold = TMath::Sqrt(2)*5./6.*25.,
  m_MeanTime = 0,
  m_nMeanCrystal = 0 
{
  m_cidVec.reserve( 300 );
  m_cTimeDeltaVec.reserve( 300 );
  m_cRvec.reserve( 300 );
  m_cPhiVec.reserve( 300 );

}  
  ClusterTime::~ClusterTime()
{  
  ;
}

bool ClusterTime::operator <  (const ClusterTime& clus ) const
{
  if( compare( clus ) < 0){
    return true;
  }else{
    return false;
  }
}
bool ClusterTime::operator == (const ClusterTime& clus ) const 
{
  if( compare( clus ) == 0){
    return true; 
  }else{ 
    return false;
  }
}
std::ostream& operator     << ( std::ostream& out, const Cluster& clus )
{
  out << "ClusterTime"  << "\n"
      << "id : "        << clus.m_id         << "\n"
      << "status : "    << clus.m_status     << "\n"
      << "Rthreshold :" << clus.m_Rthreshold << "\n"
      << "MeanTime: "   << clus.m_MeanTime   << "\n"
      << std::flush;
  out << "idVec : " ;
  for( int i = 0; i < clus.m_cidVec.size() ;i++){
    out << clus.m_cidVec[i] << " " ;
  }
  out << "\n";
  out << "RVec : " ;
  for( int i = 0; i < clus.m_cRVec.size(); i++){
    out << clus.m_cRVec[i] << " " ;
  }
  out << "\n";
  out << "PhiVec : " ; 
  for( int i = 0; i < clus.m_cPhiVec.size(); i++){
    out << clus.m_cPhiVec[i] << " ";
  }
  out << "\n" << std::flush;
  return out; 
}
int ClusterTime::Compare( const ClusterTime& clus ) const 
{
  if (m_MeanTime < clus.m_MeanTime ){
    return 1; 
  }else if( m_energy > clus.m_energy ){
    return -1;
  }
  return 0; 
}

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

class ClusterTimeAnalyzer {

  ClusterTimeAnalyzer();
  virtual ~ClusterTimeAnalyzer();  
  const Double_t RadiusThreshold;
  std::list<ClusterTime>& ClusterTimeAnalyzer::Convert( std::list<Cluster> clist );
  Double_t GetMeanTime( Cluster clus );
 private:
  Double_t Posx[3000];
  Double_t Posy[3000];
  Double_t Width[3000];
  
};

ClusterTimeAnalyzer::ClusterTimeAnalyzer() : RadiusThreshold( 25.*5./6.*TMath::Sqrt(2)){
  // Copy X,Y positionArray //
  Int_t IDArr[2716];
  for( int IDIndex = 0; IDIndex < 2716; IDIndex++){
    IDArr[IDIndex] = IDIndex; 
  }
  CsiMap::getCsiMap()->getXYWarray( 2716, IDArr, Posx, Posy, Width ); 
}
ClusterTimeAnalyzer::~ClusterTimeAnalyzer(){
  ;
}
std::list<ClusterTime>& ClusterTimeAnalyzer::Convert(std::list<Cluster> clist)
{
  m_clusterTimeList.clear();
  for( std::list<Cluster>::iterator itCluster;
       itCluster != clist.end();
       itCluster++){
    ClusterTime clusTime;        
    ConvertData( clus, clusTime );
    m_clusterTimeList.push_back(clusTime);    
  }
  return m_clusterTimeList;
}

int  ClusterTimeAnalyzer::ConvertData( Cluster clus, ClusterTime& clustTime){
  Double_t ClusterR   = 0.;
  Double_t ClusterPhi = 0.;
  Double_t coex = clus.x();
  Double_t coey = clus.y();
  std::vector<int>    IDVec = clus.clusterIdVec();
  std::vector<double> EVec  = clus.clusterEVec();
  std::vector<double> TVec  = clus.clusterTimeVec();
  std::vector<double> RVec;
  std::vector<double> PhiVec;

  Int_t    IDInCluster    = 0;
  Double_t TInCluster  = 0.;
  Double_t EInCluster  = 0.;
  Double_t MeanTime    = 0.;
  Int_t    nCenterCrystal = 0;

  ClusterR = TMath::Sqrt( coex*coex + coey*coey );
  ClusterPhi = TMath::ATan2( coey , coex );

  for( int vecIndex = 0; vecIndex < IDVec.size(); vecIndex++){
    Double_t CrystalR   =  TMath::Sqrt(( Posx[ IDVec[vecIndex] ] - coex )*( Posx[ IDVec[vecIndex] ] - coex ) + ( Posy[ IDVec[vecIndex] ] - coey )*( Posy[ IDVec[vecIndex] ] ));
    Double_t CrystalPhi =  TMath::ATan2( (Posy[ IDVec[vecIndex] ] - coey ) , (Posx[ IDVec[vecIndex] ] - coex ) );
    Double_t CrystalDeltaPhi = 0.;
    if( TMath::Abs( CrystalPhi - ClusterPhi ) > TMath::Pi() ){
      if( CrystalPhi - ClusterPhi < 0 ){ 
	CrystalDeltaPhi = 2*TMath::Pi() + ( CrystalPhi - ClusterPhi );
      }else{ //CrystalPhi -Clusterphi  > 0
	CrystalDeltaPhi = -2*TMath::Pi() + ( CrystalPhi - ClusterPhi );
      }
    }
    RVec.push_back( CrystalR );
    PhiVec.push_back( CrystalDeltaPhi);    
    if( DistanceFromCOE < RadiusThrehold ){
      MeanTime += TVec[ vecIndex ] ;
      nCenterCrystal++; 
    } 
  }
  MeanTime = MeanTime/nCenterCrystal;

  clusterTime.SetClusterMeanTime( MeanTime );
  clusterTime.SetClusterR( ClusterR );
  clusterTime.SetClusterPhi( ClusterPhi );
  for( int vecIndex  =0; vecIndex < IDVec.size(); vecIndex ++){
    clusterTime.SetCrystalData(  IDVec[vecIndex], TVec[vecIndex]-MeanTime, RVec[vecIndex], PhiVec[vecIndex]);
  }
  return IDVec.size();
}




