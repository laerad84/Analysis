// -*- C++ -*-

#include <iostream>
#include <list>
#include "E14Fsim/E14FsimFunction.h"
#include "cluster/ClusterFinder.h"
#include "csimap/CsiMap.h"

struct CSI{

  int ID;
  double Edep;
  double time;   

};             


bool 
csi_compare_energy( CSI& a, CSI& b )
{
  if( a.Edep > b.Edep )
    return( true );
  else 
    return( false );
}

/////
ClusterFinder::ClusterFinder( double threshold,int debugLevel )
  : m_debugLevel( debugLevel )
{
  m_threshold=threshold;// E deposit threshold [MeV] def=3
  m_conneDist=25.*2.8; // distance [mm] def=;
  constructKS();
}

/////
ClusterFinder::~ClusterFinder()
{
  deleteKS();
}


/*
double ClusterFinder::showerLeakage(double thre,double E,double ws,double wl){
  static int isWarning = true;
  static double const par[5][6]
    ={{-0.0124228, 0.0535973, -0.00570087, -0.0272626, 0.0492056, -0.0143836},
      {0.0135769, 0.0516771, 0.00738265,-0.0200259, 0.0532812, -0.00951445},
      {0.0451242, 0.0446529, 0.0220706,-0.00882033, 0.0543621, -0.00272351},
      {0.0632153, 0.050921, 0.0275368,0.0051316, 0.0534536, 0.00466786},
      {0.0705192, 0.063754, 0.0273419,0.0162835, 0.0534304, 0.010156} };
  static double const parthre[5] = {1.5,3.0,5.0,7.5,10};
  
  int ID = 0;
  double weight = 1;
  
  if(thre<parthre[0]||thre>=parthre[4]){
    if(isWarning){
      std::cout<<"warning : CsI energy threshold=="<<thre<<". The function for cluster energy correction is not prepared in this range"<<std::endl;
      isWarning = false;
    }
    if(thre>=parthre[4]){
      ID = 3;
      weight = 0;
    }
  }else{
    bool isFound = false;
    for(int i=1;i<5;i++){
      if(thre<parthre[i]){
	ID = i-1;
	weight = (parthre[i]-thre)/(parthre[i]-parthre[i-1]);
	isFound = true;
	//	std::cout<<"i:"<<i<<" ID:"<<ID<<std::endl;
	break;
      }
    }
    if( !isFound ){
      std::cout<<"Cluserer::showerLeakage() : Error"<<std::endl;
      std::cout<<thre<<" "<<E<<std::endl;
      exit(1);
    }
  }
  //  std::cout<<"check ID:"<<ID<<std::endl;
  double f0 = ws*(par[ID][0]+par[ID][1]/sqrt(E/1000.)-par[ID][2]*log(E/1000.))
    +wl*(par[ID][3]+par[ID][4]/sqrt(E/1000.)-par[ID][5]*log(E/1000.));
  ID++;
  double f1 = ws*(par[ID][0]+par[ID][1]/sqrt(E/1000.)-par[ID][2]*log(E/1000.))
    +wl*(par[ID][3]+par[ID][4]/sqrt(E/1000.)-par[ID][5]*log(E/1000.));

  //  std::cout<<"thre:"<<thre<<" weight:"<<weight<<" ID:"<<ID<<std::endl;
  double f = weight*f0+(1-weight)*f1;
  
  return f;
  }
*/
std::list<Cluster>& 
ClusterFinder::findCluster(TClonesArray const *csiDigi )
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
  return findCluster(nDigi,ID,E,time );
}

std::list<Cluster>& 
ClusterFinder::findCluster(int nDigi,int *ID,double *E,double *time )
{
  m_clusterList.clear();
  
  std::list<int> templist;
  for(int i=0;i<nDigi;i++){
    if(E[i]>m_threshold){
      templist.push_back(i);
    }
  }
  double posx[3000]={0};
  double posy[3000]={0};
  double width[3000]={0};
  CsiMap::getCsiMap()->getXYWarray(nDigi,ID,posx,posy,width);
  
  std::list<int> indexList;
  std::list<CSI> csiList;
  CSI csiData;

  while(templist.size()>0){
    Cluster clstr;
    
    indexList.clear();
    csiList.clear();
    
    int firstIndex = templist.front();
    indexList.push_back(firstIndex);
    csiData.ID=ID[firstIndex];
    csiData.Edep=E[firstIndex];
    csiData.time=time[firstIndex];
    csiList.push_back(csiData);
    double minTime = csiData.time;
    
    double sumEdep = E[firstIndex];
    double coex = sumEdep*posx[firstIndex];
    double coey = sumEdep*posy[firstIndex];
    templist.pop_front();

    for(std::list<int>::iterator i=indexList.begin();i!=indexList.end();i++){
      double posX = posx[*i];
      double posY = posy[*i];
      for(std::list<int>::iterator j=templist.begin();j!=templist.end();j++){
        double tempX = posx[*j];
        double tempY = posy[*j];

	double distX = fabs(tempX-posX);
	double distY = fabs(tempY-posY);
	if(distX>m_conneDist||distY>m_conneDist) continue;
	csiData.ID=ID[*j];
	csiData.Edep=E[*j];
	csiData.time=time[*j];

        double tempE = csiData.Edep = E[*j];
        sumEdep+=tempE;
        coex+=tempE*tempX;
        coey+=tempE*tempY;

	indexList.push_back(*j);
	csiList.push_back(csiData);
	if(minTime>csiData.time) minTime = csiData.time;
	std::list<int>::iterator dj = j--;
        templist.erase(dj);
      }
    }
    clstr.setThreshold(m_threshold);
    clstr.setEnergy(sumEdep);
    //    clstr.setClusterIdList(indexList);
    csiList.sort( csi_compare_energy );
    std::vector<int> idVec;
    std::vector<double> eVec;
    std::vector<double> tVec;
    for(std::list<CSI>::iterator i=csiList.begin();i!=csiList.end();i++){
      idVec.push_back(i->ID);
      eVec.push_back(i->Edep);
      tVec.push_back(i->time);
    }
    clstr.setClusterIdVec(idVec);
    clstr.setClusterEVec(eVec);
    clstr.setClusterTimeVec(tVec);
    
    clstr.setTime(minTime);

    double zsurface = CsiMap::getCsiMap()->getZSurface(clstr.clusterIdVec().at(0));
    clstr.setPos(coex/sumEdep,coey/sumEdep,zsurface);

    double RMS=0;
    for(std::list<int>::iterator i=indexList.begin();i!=indexList.end();i++){
      double R = sqrt(pow(clstr.x()-posx[*i],2)+ pow(clstr.y()-posy[*i],2));
      RMS += E[*i]*R/sumEdep;
    }
    
    clstr.setRms(RMS);

    if(csiList.size()>1 ){
      m_clusterList.push_back(clstr);
    }
  }

  m_clusterList.sort();// sorting : cedep[0]>cedep[1]>cede[2]>... 
  int id=0;
  for( std::list<Cluster>::iterator clus=m_clusterList.begin();
       clus!=m_clusterList.end(); clus++ ){
    clus->setId(id++);
  }

  return m_clusterList;
}





  



