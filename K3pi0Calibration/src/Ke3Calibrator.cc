#include "Ke3Calibrator.h"

Ke3Calibrator::Ke3Calibrator(int numRequest)
  : m_numRequest(numRequest)
{
  for(int i=0;i<s_maxArrSize;i++){
    m_calibFactor[i] = 0;
    m_counter[i] = 0;
    m_matrix[i] = new double[s_maxArrSize+1];
    for(int j=0;j<s_maxArrSize+1;j++){
      m_matrix[i][j] = 0;
    }
  }
}

Ke3Calibrator::~Ke3Calibrator(){
  delete m_calibFactor;
  for(int i=0;i<s_maxArrSize;i++)  delete m_matrix[i];
}


void Ke3Calibrator::fillGamma(Gamma const &gam,double refEnergy,double sigmaUser){
  std::vector<int> const &idvec = gam.clusterIdVec();
  std::vector<double> const &evec = gam.clusterEVec();
  int size = idvec.size();
  double EGeV = gam.e()/1000.;
  double correctFactor = gam.e()/gam.edep();
  double sigma = (sigmaUser>0)?sigmaUser:sqrt(0.01*0.01*EGeV*EGeV+0.02*0.02*EGeV)*1000;
  double sigma2 = sigma*sigma;
  for(int icsi0=0;icsi0<size;icsi0++){
    int ID0 = idvec.at(icsi0);
    double e0 = evec.at(icsi0);
    m_counter[ID0]++;
    for(int icsi=0;icsi<size;icsi++){
      int IDi = idvec.at(icsi); 
      double e = evec.at(icsi); 
      m_matrix[ID0][IDi]+=e*e0*correctFactor/sigma2;//*correctFactor;
    }
    m_matrix[ID0][3000]+=refEnergy*e0/sigma2;//*correctFactor;
  }
}


bool Ke3Calibrator::gaussianElimination(){
  std::cout<<"Ke3Calibrator : start gaussian elimination"<<std::endl;
  
  int OrigToNewID[s_maxArrSize]={0},NewToOrigID[s_maxArrSize]={0};
  int nparam=0;
  for(int i=0;i<s_maxArrSize;i++){
    if(m_counter[i]>m_numRequest){
      OrigToNewID[i] = nparam;
      NewToOrigID[nparam] = i;
      nparam++;
    }
  }

  std::cout<<"# of CsI:"<<nparam<<std::endl;
  double* matrix[s_maxArrSize];
  
  for(int i=0;i<nparam;i++){
    matrix[i] = new double[nparam+1];
    for(int j=0;j<nparam;j++)
      matrix[i][j]=m_matrix[NewToOrigID[i]][NewToOrigID[j]];
    matrix[i][nparam]=m_matrix[NewToOrigID[i]][s_maxArrSize];
  }
  
  for(int i=0;i<nparam;i++){
    if(i%100==0) std::cout<<"proceed "<<i<<" / "<<nparam<<std::endl;
    if(matrix[i][i]==0){
      std::cout<<"gaussian elimination failed : this is not a regular matrix."<<std::endl;
      return false;
    }
    for(int icol=nparam;icol>=i;icol--){
      matrix[i][icol]/=matrix[i][i];
      if(matrix[i][icol]==0) continue;
      for(int irow=0;irow<nparam;irow++){
	if(irow==i) continue;
	if(fabs(matrix[irow][icol]-=matrix[i][icol]*matrix[irow][i])<1e-10)
	  matrix[irow][icol]=0;
      }
    }
  }
  for(int i=0;i<nparam;i++){
    int ID = NewToOrigID[i];
    if(!isnan(matrix[i][nparam]))
      m_calibFactor[ID]=matrix[i][nparam];
    else{
      std::cout<<"nan:"<<" i="<<i<<std::endl;
    }
  }

  
  for(int i=0;i<nparam;i++) delete [] matrix[i];
  for(int i=0;i<s_maxArrSize;i++){
    m_counter[i] = 0;
    for(int j=0;j<s_maxArrSize+1;j++){
      m_matrix[i][j] = 0;
    }
  }

  std::cout<<"finish gaussian elimination"<<std::endl;
  return true;
}
