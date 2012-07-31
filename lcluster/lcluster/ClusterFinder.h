// -*- C++ -*-
//
// ClusterFinder.h
// csi clustering
//
// Author:  Kazufumi Sato
// Created: Fri Aug  26 10:48:00 JST 2011
//
// $Id: ClusterFinder.h,v 1.1 $
//
// $Log: ClusterFinder.h,v $
//
// Revision 1.1  2011/08/26 Sato
// New for cluster finding.
//


#ifndef CLUSTERFINDER_H_INCLUDED
#define CLUSTERFINDER_H_INCLUDED

#include <list>
#include <vector>
#include <map>
#include "CLHEP/Vector/ThreeVector.h"
#include "cluster/Cluster.h"
#include "TClonesArray.h"
#include "GsimData/GsimDigiData.h"


class ClusterFinder {
public:
  // constructor
  ClusterFinder(double threshold = 3, int debugLevel=0 );
  
  // destructor
  ~ClusterFinder();
  
  // extractor
  std::list<Cluster>&   clusterList() { return m_clusterList; }
  
  // method
  std::list<Cluster>&   findCluster( TClonesArray const *tcl );  // GsimDigiData
  std::list<Cluster>&   findCluster( int nDigi,int *ID,double *E,double *time);
  double getThreshold() const{ return m_threshold;}
  void setThreshold(double thre){ m_threshold = thre;}
  
  std::list<Cluster>&  findClusterKS(  TClonesArray const *tcl );// GsimDigiData
  std::list<Cluster>&  findClusterKS( int num,int *ID,double *E,double *time ); 

private:
  // data member
  int  m_debugLevel;
  double  m_threshold;
  double  m_conneDist;
  
  std::list<Cluster>   m_clusterList;
  


  ///////// followings are for Shiomi_version.//////////
private:
  void initKS(std::map<int,int>&XYtoIDmap,std::map<int,int>& IDtoXYmap);
  
  void setCsIDigiIDArray(int count,int id) { m_csi_id[count]=id;};
  void setCsIEneArray(int count,double energy) { m_csi_energy[count]=energy;};
  void setCsITimeArray(int count,double time) { m_csi_time[count]=time;};
  void setCsIPosArray(int count, double x,double y,double z) {m_csi_pos[count]=CLHEP::Hep3Vector(x,y,z); };
  void setCsIWidthArray(int count,double width) { m_csi_width[count]=width;};  
  int getClusterNumber() { return m_nCluster;};
    
private:
  // data member
  void constructKS();
  void deleteKS();
  int                 m_ncsi;
  int                 m_nCluster;
  
  int*                m_csi_id;
  CLHEP::Hep3Vector*         m_csi_pos;
  double*             m_csi_dxy;
  double*             m_csi_energy;
  double*             m_csi_time;
  double*             m_csi_width;
  std::list<int>*     m_csi_idlist;
  /////////////////////////////////////////////////////


};

//
#endif // CLUSTERFINDER_H_INCLUDED




