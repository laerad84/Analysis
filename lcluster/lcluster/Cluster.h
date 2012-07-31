// -*- C++ -*-
//
// Cluster.h
// cluster container class ( imported from E391a library)
//
// Author:  Kazufumi Ssato
// Created: Sun Nov 25 16:02 JST 2011
//
// $Id: Cluster.h,v 1.1 2011/11/25 16:02 sato Exp $
//
//
// $Log: Cluster.h,v $
//
//

#ifndef CLUSTER_H_INCLUDED
#define CLUSTER_H_INCLUDED

#include <iostream>
#include <list>
#include <vector>
#include <cstdlib>
#include "GsimPersistency/GsimMessage.h"
#include "CLHEP/Vector/ThreeVector.h"

class Cluster {
public:
  // constructor
  Cluster();

  // destructor
  ~Cluster();

  // operator
  bool    operator<( const Cluster& ) const;
  bool    operator==( const Cluster& ) const;

  // extractor
  int                      id()     const { return m_id; }
  int                      status() const { return m_status; }
  double                   threshold() const{ return m_threshold; }
  double                   e()      const { return m_energy; }
  const CLHEP::Hep3Vector& pos()    const { return m_pos; }
  double                   x()      const { return m_pos.x(); }
  double                   y()      const { return m_pos.y(); }
  double                   z()      const { return m_pos.z(); }
  double                   t()      const { return m_time; }
  double                   rms()    const { return m_rms; }

  //  const std::list<int>&    clusterIdList() const { return m_cidList; }
  const std::vector<int>&    clusterIdVec() const { return m_cidVec; }
  const std::vector<double>&    clusterEVec() const { return m_eVec; }
  const std::vector<double>&    clusterTimeVec() const { return m_timeVec; }

  // method
  void  setId( int id ) { m_id = id; }
  void  setStatus(int stat){ m_status  = stat; }
  void  setThreshold(double thre){ m_threshold = thre; }
  void  setEnergy( double energy ){m_energy = energy;}
  void  setPos( const CLHEP::Hep3Vector& pos ){m_pos = pos;}
  void  setPos( double x, double y, double z )
  { setPos(CLHEP::Hep3Vector( x,y,z )); }
  void  setTime( double time ){m_time = time;}
  void  setRms( double rms ) { m_rms = rms; }
  void  setClusterIdList( const std::list<int>& cidList ){m_cidList = cidList;}
  void  setClusterIdVec( const std::vector<int>& cidVec ){m_cidVec = cidVec;}
  void  setClusterEVec( const std::vector<double>& eVec ){m_eVec = eVec;}
  void  setClusterTimeVec( const std::vector<double>& timeVec ){m_timeVec = timeVec;}


  // friend functions
  friend std::ostream& operator<<( std::ostream& out, const Cluster& clus );

private:
  int            m_id;
  int            m_status;
  double         m_threshold;
  double         m_energy;
  CLHEP::Hep3Vector     m_pos;
  double         m_time;
  double         m_rms;
  std::vector<int> m_cidVec;
  std::vector<double> m_eVec;
  std::vector<double> m_timeVec;

  std::list<int> m_cidList; // <- for ClusterFinderKoShi

  int     compare( const Cluster& ) const;
};

//
#endif // CLUSTER_H_INCLUDED

