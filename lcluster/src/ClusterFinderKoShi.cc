// -*- C++ -*-
//
// ClusterFinderKoShi.cc
// Wrapper for CsI cluster finding routine by K.Shiomi
//
// Author:  Koji Shiomi
//
// $Id: ClusterFinderKoShi.cc,v 1.1 2011/12/25  Kazufumi Sato Exp $
//
// $Log: ClusterFinderKoShi.cc,v $
//
// Revision 1.1  2004/05/29 08:30:10  Sato
// clustering routine originally made by Shiomi Koji.
// In this method, csi threshold is fixed to 5 MeV.
//

#include <iostream>
#include <list>
#include <vector>
#include "cluster/Cluster.h"
#include "cluster/ClusterFinder.h"
#include "csimap/CsiMap.h"


void ClusterFinder::constructKS(){
  m_ncsi = 0;
  m_nCluster = 0;
  int const ncsi = CsiMap::getCsiMap()->getN();
  m_csi_id     = new int [ncsi];
  m_csi_pos    = new CLHEP::Hep3Vector[ncsi];
  m_csi_dxy    = new double[ncsi];
  m_csi_idlist = new std::list<int>[ncsi];
  m_csi_energy = new double[ncsi];
  m_csi_time   = new double[ncsi];
  m_csi_width   = new double[ncsi];

  for( int i=0; i<ncsi; i++ ) {
    m_csi_pos[i]    = CLHEP::Hep3Vector( 0,0,0 );
    m_csi_dxy[i]    = 0.;
    m_csi_idlist[i].clear();
    m_csi_energy[i] = 0.;
    m_csi_time[i]   = 0.;
  }

}


void ClusterFinder::deleteKS(){
  delete [] m_csi_id;
  delete [] m_csi_pos;
  delete [] m_csi_dxy;
  delete [] m_csi_idlist;
  delete [] m_csi_energy;
  delete [] m_csi_time;
  delete [] m_csi_width;
}

/////
std::list<Cluster>& 
ClusterFinder::findClusterKS(TClonesArray const *csiDigi )
{
  int nDigi = csiDigi->GetEntries();
  static int ID[3000];
  static double E[3000];
  static double time[3000];
  for(int i=0;i<nDigi;i++){
    GsimDigiData* digi = (GsimDigiData*)csiDigi->At(i);
    ID[i] = digi->modID;
    E[i] = digi->energy;
    time[i] = digi->time;
  }
  return findClusterKS(nDigi,ID,E,time );
}



std::list<Cluster>&
ClusterFinder::findClusterKS(int nDigi,int *ID,double *E,double *time )
{
  ////// add by sato //////
  m_ncsi = nDigi;
  m_nCluster = 0;
  for(int i=0; i<m_ncsi; i++){
    double x,y,width;
    CsiMap::getCsiMap()->getXYW(ID[i],x,y,width);
    //    double z = 6148; // mm
    double z = CsiMap::getCsiMap()->getZSurface(ID[i]);
    
    setCsIDigiIDArray(i,ID[i]);
    setCsIEneArray(i,E[i]);
    setCsITimeArray(i,time[i]);
    setCsIPosArray(i,x,y,z);   
    setCsIWidthArray(i,width);
  }
  /////////////////////////
  
  m_threshold=5;// E deposit threshold [MeV]
  m_conneDist=25.3*sqrt(8.5); // distance [mm]
  m_clusterList.clear();

  double  Csir_Len     = 180; // CsI region
  double  Csir_Zmin    = 609.8; // CsI region (updated : Mon Jun 21 22:45:34\ JST 2004)
  double  Csic_Zmin    = 614.8; // CsI face   (updated : Mon Jun 21 22:45:34\ JST 2004)
  double  Csic_XY      = 7.0;
  double  Csic_Len     = 30.0;
  double  Csic_Gap     = 0.03;
  double  KtevCsic_XY  = 5.0;
  double  KtevCsic_Len = 50.0;
  double  KtevCsic_Gap = 0.03;

  double  KtevCsicFine_XY  = 2.5;
  double  KtevCsicFine_Len = 50.0;
  double  KtevCsicFine_Gap = KtevCsic_Gap*2.5/5.;

  int nBlockFine=48;
  int nBeamHole=8;
  int nLayer=38;

  static std::map<int,int> XYtoIDmap;
  static std::map<int,int> IDtoXYmap;
  static int init=0;
  
  if(init==0){
    initKS(XYtoIDmap,IDtoXYmap);
    init=1;
  }

  std::list<int> templist;
  for(int i=0;i<m_ncsi;i++){
    if(m_csi_energy[i]>m_threshold){
      templist.push_back(i);
    }
  }
  std::list<int> indexList;
  std::list<int> clusterIDList;
  std::list<double> timeList;
  
  while(templist.size()>0){

    Cluster clstr;

    indexList.clear();
    clusterIDList.clear();
    timeList.clear();

    int firstIndex = templist.front();
    indexList.push_back(firstIndex);
    clusterIDList.push_back(firstIndex);
    timeList.push_back(m_csi_time[firstIndex]);
    
    double sumEdep = m_csi_energy[firstIndex];
    double coex = sumEdep*m_csi_pos[firstIndex].x();
    double coey = sumEdep*m_csi_pos[firstIndex].y();
    templist.pop_front();
    
    int size=1;
    for(std::list<int>::iterator i=indexList.begin();i!=indexList.end();i++){
      int stackID,stackXID,stackYID,startxid,startyid,endxid,endyid;
      int tempID,tempXID,tempYID;
      int regionflag;

      std::map<int,int>::iterator stackmap;
      stackmap = IDtoXYmap.find(m_csi_id[*i]);
      if( stackmap != IDtoXYmap.end()) stackID= stackmap->second;
      else std::cout<<"CsI IDtoXYmap mapping error"<<std::endl;
      
      for(std::list<int>::iterator j=templist.begin();j!=templist.end();j++){		

	std::map<int,int>::iterator tempmap;
	tempmap = IDtoXYmap.find(m_csi_id[*j]);
	if( tempmap != IDtoXYmap.end()) tempID= tempmap->second;
	else std::cout<<"CsI IDtoXYmap mapping error"<<std::endl;

        double tempX = m_csi_pos[*j].x();
        double tempY = m_csi_pos[*j].y();
	
	if(m_csi_id[*i]<2240){
	  if(m_csi_id[*j]<2240){
	    stackXID = stackID%nBlockFine;
	    stackYID = int(stackID/nBlockFine);		    	    
	    tempXID = tempID%nBlockFine;
	    tempYID = int(tempID/nBlockFine);		    	    
	    
	    regionflag=0;
	  } else {
	    stackXID = (stackID%nBlockFine)/2+7;
	    stackYID = (stackID/nBlockFine)/2+7;	    
	    tempXID = tempID%nLayer;
	    tempYID = tempID/nLayer;     
	    regionflag=1;
	  }
	} else{
	  if(m_csi_id[*j]<2240){
	    stackXID = stackID%nLayer;
	    stackYID = stackID/nLayer;     
	    tempXID =  (tempID%nBlockFine)/2+7;
	    //	    tempYID =  (tempID%nBlockFine)/2+7;
	    tempYID =  (tempID/nBlockFine)/2+7;
	  } else {
	    stackXID = stackID%nLayer;
	    stackYID = stackID/nLayer;     
	    tempXID =  tempID%nLayer;
	    tempYID =  tempID/nLayer;	    
	  }
	  regionflag=2;
	}
	
	if(regionflag==0){
	  startxid=stackXID-2;
	  startyid=stackYID-2;
	  endxid=stackXID+2;
	  endyid=stackYID+2;
	} else {
	  startxid=stackXID-1;
	  startyid=stackYID-1;
	  endxid=stackXID+1;
	  endyid=stackYID+1;
	}


	if(startxid<=tempXID&&tempXID<=endxid){
	  if(startyid<=tempYID&&tempYID<=endyid){	    	    
	    
	    double tempE = m_csi_energy[*j];
	    sumEdep+=tempE;
	    coex+=tempE*tempX;
	    coey+=tempE*tempY;
	    
	    size++;
	    
	    indexList.push_back(*j);
	    clusterIDList.push_back(*j);
	    timeList.push_back(m_csi_time[*j]);

	    std::list<int>::iterator dj = j--;
	    templist.erase(dj);
	  }
	}
      }
    }
    clstr.setEnergy(sumEdep);
    clstr.setClusterIdList(clusterIDList);
    
    timeList.sort();
    //std::list<double>::iterator p = timeList.end();
    std::list<double>::iterator p = timeList.begin();

    
    clstr.setTime(*p);
    
    
    //static double const param[2]={0.00361837,0.037974};
    //double F = param[0]+param[1]/sqrt(sumEdep/1000);
    //clstr.setEnergy((1+F)*sumEdep);
    //double f = 1+0.06334+0.01815/sqrt(sumEdep/1000.)-0.01634*log(sumEdep/1000.);
    //clstr.setEnergy(f*sumEdep);
    //clstr.setEnergy(sumEdep);

    //    clstr.setPos(coex/sumEdep,coey/sumEdep,6148);
    clstr.setPos(coex/sumEdep,coey/sumEdep,
		 m_csi_pos[clusterIDList.front()].z()); // change by sato

    double RMS=0;
    std::vector<int> idVec;
    std::vector<double> eVec;
    std::vector<double> tVec;
    for(std::list<int>::iterator i=clusterIDList.begin();i!=clusterIDList.end();i++){
      double R = sqrt(pow(clstr.x()-m_csi_pos[*i].x(),2)+	         
		      pow(clstr.y()-m_csi_pos[*i].y(),2));
      RMS += m_csi_energy[*i]*R/sumEdep;

      ////// add by sato/////////
      eVec.push_back(m_csi_energy[*i]);
      tVec.push_back(m_csi_time[*i]);
      idVec.push_back(m_csi_id[*i]);
      ///////////////////////////
    }
    clstr.setRms(RMS);
    ////// add by sato/////////
    clstr.setClusterIdVec(idVec);
    clstr.setClusterEVec(eVec);
    clstr.setClusterTimeVec(tVec);
    clstr.setThreshold(m_threshold);
    ////////////////////////////

    
    if(clusterIDList.size()>1){
      m_clusterList.push_back(clstr);
      m_nCluster++;
    }
  }
  
  m_clusterList.sort();// sorting : cedep[0]>cedep[1]>cede[2]>... 
  return m_clusterList;
}

    
//if( m_debugLevel ) {
//  std::cout << "clist.size() = " << clist.size() << std::endl;
//  for( std::list<CSI>::iterator i=clist.begin();
//       i!=clist.end(); i++ ) {
//    int    id = i->id;
//    double e  = i->energy;
//    printf("%d %1.4e\n",id,e);
//  }
//}



/////
void
ClusterFinder::initKS(std::map<int,int>&XYtoIDmap,std::map<int,int>& IDtoXYmap)
{
  //
  if( m_debugLevel )
    std::cout << "ClusterFinder::initKoShi()2" << std::endl;


  double  Csir_Len     = 180; // CsI region
  double  Csir_Zmin    = 609.8; // CsI region (updated : Mon Jun 21 22:45:34\ JST 2004)
  double  Csic_Zmin    = 614.8; // CsI face   (updated : Mon Jun 21 22:45:34\ JST 2004)
  double  Csic_XY      = 7.0;
  double  Csic_Len     = 30.0;
  double  Csic_Gap     = 0.03;
  double  KtevCsic_XY  = 5.0;
  double  KtevCsic_Len = 50.0;
  double  KtevCsic_Gap = 0.03;
  
  double  KtevCsicFine_XY  = 2.5;
  double  KtevCsicFine_Len = 50.0;
  double  KtevCsicFine_Gap = KtevCsic_Gap*2.5/5.;
  
  
  int nBlockFine=48;
  int nBeamHole=8;
  int id=0;
  for(int i=0;i<nBlockFine;i++) {
    
    for(int j=0;j<nBlockFine;j++) {
      {
	//beam hole
	if( i>=(nBlockFine-nBeamHole)/2. &&
	    i<(nBlockFine-nBeamHole)/2.+nBeamHole &&
	    j>=(nBlockFine-nBeamHole)/2. &&
	    j<(nBlockFine-nBeamHole)/2.+nBeamHole )
	  continue;
	
	
	IDtoXYmap.insert(std::make_pair(id,j+i*nBlockFine));
	XYtoIDmap.insert(std::make_pair(j+i*nBlockFine,id));
	id++;
      }
    }
  }

  int csiStack[38]
    = { 0,12,16,20,22,24,
        26,28,30,32,32,34,
        34,36,36,36,36,36,36,
        36,36,36,36,36,36,34,
        34,32,32,30,28,26,
        24,22,20,16,12,0 };
  
  int nLayer=38;
  std::cout<<id<<std::endl;
  
  for(int i=0;i<nLayer;i++) {
    int nBlock = csiStack[i];
    
    for(int j=0;j<nBlock;j++) {
      
      //KtevFine
      if( i>=(nLayer-nBlockFine/2.)/2. &&
          i<(nLayer-nBlockFine/2.)/2.+nBlockFine/2. &&
          j>=(nBlock-nBlockFine/2.)/2. &&
          j<(nBlock-nBlockFine/2.)/2.+nBlockFine/2. )
        continue;
      
      
      IDtoXYmap.insert(std::make_pair(id,(nLayer-nBlock)/2+j+i*nLayer));
      XYtoIDmap.insert(std::make_pair((nLayer-nBlock)/2+j+i*nLayer,id));
      
      
      id++;
    }
  }
  
}
  

