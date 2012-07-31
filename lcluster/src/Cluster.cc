
#include <iostream>
#include <list>
#include <math.h>
#include "cluster/Cluster.h"


/////
Cluster::Cluster()
  : m_id( 0 ),
    m_pos( CLHEP::Hep3Vector( 0,0,0 ) ),
    m_energy( 0.0 ),
    m_time( 0.0 ),
    m_rms(  0.0 ),
    m_status( 0 ),
    m_threshold( 0 )
{
  //
  //std::cout << "Cluster::Cluster()" << std::endl;

}


/////
Cluster::~Cluster()
{
  //
  //std::cout << "Cluster::~Cluster()" << std::endl;

}

/////
bool
Cluster::operator<( const Cluster& clus ) const 
{
  if( compare( clus ) < 0 ) {
    return( true );
  }
  return( false );

}

/////
bool
Cluster::operator==( const Cluster& clus ) const
{
  if( compare( clus ) == 0 ) {
    return( true );
  }
  return( false );

}

/////
std::ostream& operator<<( std::ostream& out, const Cluster& clus )
{
  //
  out << "Cluster::dump : " << std::endl
      << " id       = " << clus.m_id     << std::endl
      << " status   = " << clus.m_status << std::endl
      << " threshold   = " << clus.m_threshold << std::endl
      << " e      = " << clus.m_energy   << std::endl
      << " x      = " << clus.m_pos.x()  << std::endl
      << " y      = " << clus.m_pos.y()  << std::endl
      << " z      = " << clus.m_pos.z()  << std::endl
      << " t  = " << clus.m_time    << std::endl
      << " rms  = " << clus.m_rms    << std::endl;

  out << " idList   = ";
  for( std::list<int>::const_iterator p=clus.m_cidList.begin();
       p!=clus.m_cidList.end(); p++ ) {
    out << *p << " ";
  }
  out<<std::endl
     << " idVec   = ";
  for( int i=0;i<clus.m_cidVec.size();i++){
    out << clus.m_cidVec[i] << " ";
  }
  out<<std::endl
     << " eVec   = ";
  for( int i=0;i<clus.m_eVec.size();i++){
    out << clus.m_eVec[i] << " ";
  }
  out<<std::endl
     << " timeVec   = ";
  for( int i=0;i<clus.m_timeVec.size();i++){
    out << clus.m_timeVec[i] << " ";
  }

  out << std::endl;
  

  return( out );
}

/////
int 
Cluster::compare( const Cluster& clus ) const
{
  if( m_energy < clus.m_energy ) {
    return( 1 );
  }
  else if( m_energy > clus.m_energy ) {
    return( -1 );
  }

  return( 0 ); // same
}

