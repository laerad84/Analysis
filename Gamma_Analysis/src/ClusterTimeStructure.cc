#include "ClusterTimeStructure.h"

ClusterTime::ClusterTime(): m_id(-1),
			    m_status(-1),
			    m_Rthreshold(TMath::Sqrt(2)*5./6.*25.),
			    m_MeanTime(0)
{
  m_cidVec.reserve( 300 );
  m_cTimeDeltaVec.reserve( 300 );
  m_cRVec.reserve( 300 );
  m_cPhiVec.reserve( 300 );

}
ClusterTime::~ClusterTime(){
  ;
}
bool ClusterTime::operator <  (const ClusterTime& clus ) const
{
  if( Compare( clus ) < 0){
    return true;
  }else{
    return false;
  }
}
bool ClusterTime::operator == (const ClusterTime& clus ) const 
{
  if( Compare( clus ) == 0){
    return true; 
  }else{ 
    return false;
  }
}
int ClusterTime::Compare( const ClusterTime& clus ) const 
{
  if (m_MeanTime < clus.m_MeanTime ){
    return 1; 
  }else if( m_MeanTime > clus.m_MeanTime ){
    return -1;
  }
  return 0; 
}
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
int ClusterTimeAnalyzer::Convert(std::list<Cluster>& clist , std::list<ClusterTime>& ctlist)
{
  //m_clusterTimeList.clear();

  for( std::list<Cluster>::iterator itCluster= clist.begin();
       itCluster != clist.end();
       itCluster++){
    ClusterTime clusTime;  
    //ConvertData( (*itCluster), clusTime );
    
    Double_t ClusterR   = 0.;
    Double_t ClusterPhi = 0.;
    std::cout<< __LINE__ << std::endl; 
    Double_t coex = (*itCluster).x();
    Double_t coey = (*itCluster).y();
    std::cout<< __LINE__ << std::endl; 
    std::vector<int>    IDVec = (*itCluster).clusterIdVec();
    std::vector<double> EVec  = (*itCluster).clusterEVec();
    std::vector<double> TVec  = (*itCluster).clusterTimeVec();
    std::vector<double> RVec;
    std::vector<double> PhiVec;
    std::cout<< __LINE__ << std::endl; 
    Int_t    IDInCluster = 0;
    Double_t TInCluster  = 0.;
    Double_t EInCluster  = 0.;
    Double_t MeanTime    = 0.;
    Int_t    nCenterCrystal = 0;
    std::cout<< __LINE__ << std::endl; 
    ClusterR = TMath::Sqrt( coex*coex + coey*coey );
    ClusterPhi = TMath::ATan2( coey , coex );
    std::cout<< __LINE__ << std::endl; 
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
      if( CrystalR  < RadiusThreshold ){
	MeanTime += TVec[ vecIndex ] ;
	nCenterCrystal++; 
      } 
    }
    MeanTime = MeanTime/nCenterCrystal;
    
    clusTime.SetClusterMeanTime( MeanTime );
    clusTime.SetClusterR( ClusterR );
    clusTime.SetClusterPhi( ClusterPhi );
    for( int vecIndex = 0; vecIndex < IDVec.size(); vecIndex++ ){
      clusTime.SetCrystalData(  IDVec[vecIndex], TVec[vecIndex]-MeanTime, RVec[vecIndex], PhiVec[vecIndex]);
    }
    std::cout<< __LINE__ << std::endl;
    ctlist.push_back(clusTime);    
  }
  std::cout<< __LINE__ << std::endl;
  return ctlist.size();
}
int  ClusterTimeAnalyzer::ConvertData( Cluster& clus, ClusterTime& clusterTime){

  Double_t ClusterR   = 0.;
  Double_t ClusterPhi = 0.;
  std::cout<< __LINE__ << std::endl; 
  Double_t coex = clus.x();
  Double_t coey = clus.y();
  std::cout<< __LINE__ << std::endl; 
  std::vector<int>    IDVec = clus.clusterIdVec();
  std::vector<double> EVec  = clus.clusterEVec();
  std::vector<double> TVec  = clus.clusterTimeVec();
  std::vector<double> RVec;
  std::vector<double> PhiVec;
  std::cout<< __LINE__ << std::endl; 
  Int_t    IDInCluster = 0;
  Double_t TInCluster  = 0.;
  Double_t EInCluster  = 0.;
  Double_t MeanTime    = 0.;
  Int_t    nCenterCrystal = 0;
  std::cout<< __LINE__ << std::endl; 
  ClusterR = TMath::Sqrt( coex*coex + coey*coey );
  ClusterPhi = TMath::ATan2( coey , coex );
  std::cout<< __LINE__ << std::endl; 
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
    if( CrystalR  < RadiusThreshold ){
      MeanTime += TVec[ vecIndex ] ;
      nCenterCrystal++; 
    } 
  }
  MeanTime = MeanTime/nCenterCrystal;

  clusterTime.SetClusterMeanTime( MeanTime );
  clusterTime.SetClusterR( ClusterR );
  clusterTime.SetClusterPhi( ClusterPhi );
  for( int vecIndex = 0; vecIndex < IDVec.size(); vecIndex++ ){
    clusterTime.SetCrystalData(  IDVec[vecIndex], TVec[vecIndex]-MeanTime, RVec[vecIndex], PhiVec[vecIndex]);
  }
  return IDVec.size();
}



