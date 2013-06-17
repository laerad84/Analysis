#include "TChain.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TSpline.h"
#include "TCut.h"
#include "TEventList.h"
#include "TH1F.h"
#include "math.h"
#include "feb2012/global.h"
#include "csimap/CsiMap.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "gnana/E14GNAnaFunction.h"
#include "TDirectory.h"
#include "TRandom.h"
#include <iostream>
#include <string>

std::string func_chamberReso="sqrt(pow(0.0159,2)+pow(0.0273*x/1000,2))*x";
double const edepOnMaterial = 4.62;

double ChamberResolution(double const &PM){
  static TF1 func("localFunc_chamber",func_chamberReso.c_str(),0,5000);
  return func.Eval(PM); //MeV
}
double CsiResolution(double const &E){
  double EGeV = E/1000;
  return sqrt(0.01*0.01*EGeV*EGeV+0.02*0.02*EGeV)*1000;
}
int fidutialMode = 0;
bool csiFidutial(CLHEP::Hep3Vector const &pos){
  static bool isFirstTime = true;
  static bool juneConfig = false;
  if(isFirstTime){
    isFirstTime = false;
    std::string csiJuneConfig = chamber::getConfigDir()+"/csiJuneConfig.dat";
    std::ifstream ifs(csiJuneConfig.c_str());
    if(ifs){
      juneConfig = true;
      std::cout<<"csi june configuration"<<std::endl;
    }else std::cout<<"!csi june configuration"<<std::endl;
  }
  //  std::cout<<pos.x()<<" "<<pos.y()<<" "<<pos.r()<<std::endl;
  if(fabs(pos.y())<550&&pos.perp()<850&&(fabs(pos.x())>150||fabs(pos.y())>150)){
    if(!juneConfig) return true;
    if(pos.y()<150&&pos.x()>150) return true;
    if(pos.y()>150&&pos.x()>-50) return true;
  }
  return false;
}

int main(int argc,char** argv){
  std::string prefix="";
  if(argc==7){
    prefix = argv[6];
  }else if(argc!=6){
    std::cout<<"arg err"<<std::endl
	     <<"inputDir output version seed nuse[1-10] "<<std::endl;
    return 1;
  }

  double const calibDiviation = 0.05;
  int const nItr = 3;
  int const seed = atoi(argv[4]);
  int const nuse = atoi(argv[5]);
  std::cout<<"randomSeed:"<<seed<<" nuse:"<<nuse<<std::endl;
  std::string version = argv[3];
  bool isSimulation = false;
  if(version.find("ver")==std::string::npos){
    isSimulation = true;
  }
  std::cout<<"isSimulation:"<<isSimulation<<std::endl;
  chamber::setVersion(version);
  TChain *trin = new TChain("goodtr");
  TChain *trin_gam = new TChain("gammatr");
  TChain *trin_com = new TChain("comptr");
  std::string inputDir = argv[1];
  std::string fname = inputDir+"/run*.root";
  trin->Add(fname.c_str());
  std::cout<<"read "<<fname<<" # of entry:"<<trin->GetEntries()<<std::endl;
  //  fname = inputDir+"/"+prefix+"gamrun*.root";
  fname = prefix+"*.root";
  trin_gam->Add(fname.c_str());
  std::cout<<"read "<<fname<<" # of entry:"<<trin_gam->GetEntries()<<std::endl;
  fname = inputDir+"/comp*.root";
  trin_com->Add(fname.c_str());
  std::cout<<"read "<<fname<<" # of entry:"<<trin_com->GetEntries()<<std::endl;

  TFile f(argv[2],"recreate");
  std::cout<<"create "<<argv[2]<<std::endl;
  bool active[3000];
  double calibConstant[3000],x[3000],y[3000],w[3000];  
  double MeVtoCnt[3000],cosmiCalib[3000];
  for(int i=0;i<3000;i++){
    calibConstant[i] = MeVtoCnt[i] =  cosmiCalib[i] = -1;
    x[i]=y[i]=w[i] = -1;
    active[i] = false;
  }
  {
    int ncsi = CsiMap::getCsiMap()->getN();
    for(int i=0;i<ncsi;i++){
      CsiMap::getCsiMap()->getXYW(i,x[i],y[i],w[i]);
      CLHEP::Hep3Vector vec(x[i],y[i],0);
      if(csiFidutial(vec)) active[i] = true;

    }
  }
  
  {
    if(isSimulation){
      gRandom->SetSeed(seed);
      for(int i=0;i<CsiMap::getCsiMap()->getN();i++){
	calibConstant[i] = gRandom->Gaus(1,calibDiviation);
      }
    }else {
      for(int i=0;i<CsiMap::getCsiMap()->getN();i++){
	calibConstant[i] = 1.;
      }
    }
  }

  int const nloop = trin->GetEntries();
  std::cout<<"nloop:"<<nloop<<std::endl;

  E14GNAnaDataContainer clusData;
  clusData.setBranchAddress(trin_gam);
  TwoChargedEvent trackData;
  trackData.setBranchAddress(trin);
  trackData.setBranchAddress(trin_gam);
  int hitStat[2],EventNo=-1,runid=-1;
  trin_com->SetBranchAddress("hitStat",hitStat);
  trin->SetBranchAddress("specEvt",&EventNo);
  trin->SetBranchAddress("runid",&runid);
      
  TTree tr("calib","");
  int iteration = 0;
  double correction[3000];
  for(int i=0;i<3000;i++) correction[i] = 1;
  tr.Branch("nitr",&iteration,"nitr/I");
  tr.Branch("correction",correction,"correction[3000]/D");
  tr.Branch("calib",calibConstant,"calib[3000]/D");
  tr.Branch("cosmicalib",cosmiCalib,"cosmicalib[3000]/D");
  tr.Branch("x",x,"x[3000]/D");
  tr.Branch("y",y,"y[3000]/D");
  tr.Branch("w",w,"w[3000]/D");
  tr.Branch("active",active,"active[3000]/O");
  tr.Fill();
  

  TTree *evttr[100],*evttr2[100];
  for(int iitr=0;iitr<nItr;iitr++){
    evttr[iitr] = new TTree(Form("tr%d",iitr),"events used for calibration");
    //    evttr2[iitr] = new TTree(Form("tr%dB",iitr),"events not used for calibration");
    bool bad;
    int csize;
    double edep,ene,PM,X,Y,hitX,hitY,maxCnt;
    evttr[iitr]->Branch("bad",&bad,"bad/O");
    evttr[iitr]->Branch("csize",&csize,"csize/I");
    evttr[iitr]->Branch("edep",&edep,"edep/D");
    evttr[iitr]->Branch("ene",&ene,"ene/D");
    evttr[iitr]->Branch("PM",&PM,"PM/D");
    evttr[iitr]->Branch("X",&X,"X/D");
    evttr[iitr]->Branch("Y",&Y,"Y/D");
    evttr[iitr]->Branch("hitX",&hitX,"hitX/D");
    evttr[iitr]->Branch("hitY",&hitY,"hitY/D");
    evttr[iitr]->Branch("maxCnt",&maxCnt,"maxCnt/D");
    /*
    evttr2[iitr]->Branch("bad",&bad,"bad/O");
    evttr2[iitr]->Branch("csize",&csize,"csize/I");
    evttr2[iitr]->Branch("edep",&edep,"edep/D");
    evttr2[iitr]->Branch("ene",&ene,"ene/D");
    evttr2[iitr]->Branch("PM",&PM,"PM/D");
    evttr2[iitr]->Branch("X",&X,"X/D");
    evttr2[iitr]->Branch("Y",&Y,"Y/D");
    evttr2[iitr]->Branch("hitX",&hitX,"hitX/D");
    evttr2[iitr]->Branch("hitY",&hitY,"hitY/D");
    evttr2[iitr]->Branch("maxCnt",&maxCnt,"maxCnt/D");
    */
    double* tmpmatrix[3000];
    for(int i=0;i<3000;i++){
      tmpmatrix[i] = new double[3000+1];
      for(int j=0;j<3000+1;j++) tmpmatrix[i][j]=0;
    }
    int counter[3000]={0};
    int nelectron = 0;
    for(int ievt=0;ievt<nloop;ievt++){
      if(nloop>100&&ievt%(nloop/100)==0)std::cout<<ievt/(nloop/100)<<"%"<<std::endl;
      //      std::cout<<"a"<<std::endl;
      bool useCalib = false;
      //      bool fillTree = false;
      //      if(ievt%100>=70){
      //	fillTree = true;
      //      }
      if(nuse>10&&nelectron>nuse) break;
      
      trin->GetEntry(ievt);
      if(EventNo%10<=nuse){
	useCalib = true;
      }


      //      if(!fillTree&&!useCalib) continue;
      if(!useCalib) continue;

      trin_gam->GetEntry(ievt);
      trin_com->GetEntry(ievt);
      std::list<Gamma> glist;
      clusData.getData(glist);
      //      std::cout<<"z"<<std::endl;
      if(hitStat[0]!=0||hitStat[1]!=0) continue;
      //      std::cout<<"a"<<std::endl;
      if(trackData.nTrack!=2) continue;
      //      std::cout<<"b"<<std::endl;
      if(trackData.kpm0>-9000 || trackData.L<0 || trackData.L>50 ||
	 std::max(trackData.chi2[0],trackData.chi2[1])>30 || 
	 trackData.clusSize[0]>=80 || trackData.clusSize[1]>=80 ||
	 trackData.vtx[2]<-500 || trackData.vtx[2]>1400 || 
	 trackData.trigStat < 2) continue;
      //      std::cout<<"c"<<std::endl;
      //	 trackData.status!=0 ) continue;
      for(int itra=0;itra<2;itra++){
	//	std::cout<<"b"<<std::endl;
	int ID = trackData.clusId[itra];
	if(ID<0) continue;
	//	std::cout<<"ID:"<<ID<<std::endl;
	//	std::cout<<"diffPos:"<<trackData.diffPos[itra]<<std::endl;
	if(trackData.diffPos[itra]>100) continue;
	//	std::cout<<"runid,EvetNo:"<<runid<<" "<<EventNo<<std::endl;
	//	std::cout<<"glist.size():"<<glist.size()<<std::endl;
	std::list<Gamma>::iterator git = glist.begin();
	for(int i=0;i<ID;i++) git++;
	Gamma gam = *(git);
	bad = true;
	//	std::cout<<"bb"<<std::endl;
	if(gam.clusterIdVec().size()>4 && gam.chisq()<2.5 &&
	   csiFidutial(gam.pos())==true ) bad = false;
	//	std::cout<<gam.clusterIdVec().size()<<" "<<gam.chisq()<<" "<<csiFidutial(gam.pos())<<std::endl;
	//	std::cout<<"bbb"<<std::endl;
	csize = gam.clusterIdVec().size();
	if(trackData.clusSize[itra]!=csize)
	  std::cout<<"ID:"<<ID<<" clusSize:"<<trackData.clusSize[itra]<<" GamClusSize:"<<gam.clusterIdVec().size()<<std::endl;

	maxCnt = 0;
	X=gam.x();
	Y=gam.y();
	hitX=trackData.hitX[itra][0];
	//	if(seed==0) hitX+=40.*trackData.hitP[itra][0]/trackData.hitP[itra][2];
	hitY=trackData.hitX[itra][1];
	PM=trackData.pMag[itra];
	double sigChamber = ChamberResolution(PM);
	//	std::cout<<"c"<<std::endl;
	{
	  Cluster clus = gam.cluster();
	  int ndigi = clus.clusterEVec().size();
	  std::vector<double> evec;
	  double esum = 0;
	  for(int i=0;i<ndigi;i++){
	    int id = clus.clusterIdVec().at(i);
	    double e = clus.clusterEVec().at(i);
	    e = e * calibConstant[id]*correction[id];
	    esum+=e;
	    evec.push_back(e);
	  }
	  clus.setClusterEVec(evec);
	  clus.setEnergy(esum);
	  gam.setCluster(clus);
	  //	  std::cout<<"gam.e():"<<gam.e()<<" "<<gam.edep()<<std::endl;
	  E14GNAnaFunction::getFunction()->correctEnergy(gam);
	  //	  std::cout<<"2 gam.e():"<<gam.e()<<std::endl;
	  E14GNAnaFunction::getFunction()->correctEnergyWithAngle(gam);
	  //	  std::cout<<"3 gam.e():"<<gam.e()<<std::endl;
	}	 
	//	std::cout<<"d"<<std::endl;
	double sigCsi = gam.sigmaE();
	double correctFactor = gam.e()/gam.edep();
	double sigma2 = sigCsi*sigCsi+sigChamber*sigChamber;
	ene = gam.e();
	//	double scale = 1;
	//	if(seed==0) scale = 0.99; 
	//	if( bad  || fabs(ene/(PM*scale-3)-1)>0.1 ){
	double acceptRange =(iitr==0)?0.2:0.1;
	//	if( fillTree ) evttr2[iitr]->Fill();
	//	if(!useCalib) continue;
	if( bad  || fabs(ene/(PM-edepOnMaterial)-1)>acceptRange ){
	  evttr[iitr]->Fill();
	  //	  std::cout<<"e"<<std::endl;
	  continue;
	}
	//	std::cout<<"f"<<std::endl;
	nelectron++;
	{
	  std::vector<int> const &idvec = gam.clusterIdVec();
	  std::vector<double> const &evec = gam.clusterEVec();
	  int size = idvec.size();
	  for(int icsi0=0;icsi0<size;icsi0++){
	    int ID0 = idvec.at(icsi0);
	    double e0 = evec.at(icsi0);
	    counter[ID0]++;
	    for(int icsi=0;icsi<size;icsi++){
	      int IDi = idvec.at(icsi); 
	      double e = evec.at(icsi); 
	      tmpmatrix[ID0][IDi]+=e*e0*correctFactor/sigma2;//*correctFactor;
	    }
	    tmpmatrix[ID0][3000]+=(PM-edepOnMaterial)*e0/sigma2;//*correctFactor;
	  }
	}
	//	std::cout<<"z"<<std::endl;
	evttr[iitr]->Fill();
      }
    }
    std::cout<<"# of electron:"<<nelectron<<std::endl;
    // if(iitr>=nItr-2) continue;
    int OrigToThisID[3000]={0},ThisToOrigID[3000]={0};
    int nCsIWithHitOver10=0;
    for(int i=0;i<3000;i++){
      if(counter[i]>10){
	OrigToThisID[i] = nCsIWithHitOver10;
	ThisToOrigID[nCsIWithHitOver10] = i;
	nCsIWithHitOver10++;
      }
    }
    int const nparam = nCsIWithHitOver10;
    std::cout<<"# of CsI:"<<nparam<<std::endl;
    double* matrix[3000];

    for(int i=0;i<nparam;i++){
      matrix[i] = new double[nparam+1];
      for(int j=0;j<nparam;j++)
	matrix[i][j]=tmpmatrix[ThisToOrigID[i]][ThisToOrigID[j]];
      matrix[i][nparam]=tmpmatrix[ThisToOrigID[i]][3000];
    }
    
    std::cout<<"start gauss's elimination"<<std::endl;
    for(int i=0;i<nparam;i++){
      if(i%100==0) std::cout<<"proceed "<<i<<" / "<<nparam<<std::endl;
      if(matrix[i][i]==0){
	std::cout<<"kono Gyouretsu ha SEISOKU deha ARIMASEN."<<std::endl;
	std::cout<<"i:"<<i<<std::endl;
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
      int ID = ThisToOrigID[i];
      if(!isnan(matrix[i][nparam]))
	correction[ID]*=matrix[i][nparam];
      else{
	std::cout<<"nan:"<<" i="<<i<<std::endl;
      }
    }
    iteration = iitr+1;
    tr.Fill();
    for(int i=0;i<nparam;i++){
      int ID = ThisToOrigID[i];
      if(correction[ID]>0&&csiFidutial(CLHEP::Hep3Vector(x[ID],y[ID],0))==true){
	;
      }else{
	correction[ID]=1;
      }
    }
    for(int i=0;i<nparam;i++) delete [] matrix[i];
    for(int i=0;i<3000;i++) delete [] tmpmatrix[i];
    
  }
  for(int iitr=0;iitr<nItr;iitr++)  evttr[iitr]->Write();
  //  for(int iitr=0;iitr<nItr;iitr++)  evttr2[iitr]->Write();
  tr.Write();
  f.Close();
}
