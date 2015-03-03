/*
 *  $Id: DmpAlgCalMipsBar.cc, 2015-03-02 14:44:50 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 19/07/2014
*/

//#include <stdio.h>
//#include <boost/lexical_cast.hpp>
#include <algorithm>

#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"

#include "DmpEvtHeader.h"
#include "DmpEvtBgoRaw.h"
#include "DmpAlgCalMipsBar.h"
#include "DmpDataBuffer.h"
#include "DmpParameterBgo.h"
#include "DmpParameterPsd.h"
#include "DmpBgoBase.h"
#include "DmpPsdBase.h"
#include "DmpCore.h"
#include "DmpTimeConvertor.h"
#include "MyFunctions.h"

//-------------------------------------------------------------------
DmpAlgCalMipsBar::DmpAlgCalMipsBar()
 :DmpVAlg("Cal/Bgo/Ped"),
  fEvtHeader(0),
  fEvtBgo(0),
  fEvtPsd(0),
  fFirstEvtTime(-1),
  fLastEvtTime(-1)
{
  gRootIOSvc->SetOutputKey("ped");
}

//-------------------------------------------------------------------
DmpAlgCalMipsBar::~DmpAlgCalMipsBar(){
}

//-------------------------------------------------------------------
bool DmpAlgCalMipsBar::Initialize(){
  // read input data
  std::string inpath = "Event/Rec0/";
  fEvtHeader = dynamic_cast<DmpEvtHeader*>(gDataBuffer->ReadObject(inpath+"EventHeader"));
  if(0 == fEvtHeader){
    fEvtHeader = new DmpEvtHeader();
    gDataBuffer->LinkRootFile(inpath+"EventHeader",fEvtHeader);
  }
  fEvtBgo = dynamic_cast<DmpEvtBgoRaw*>(gDataBuffer->ReadObject(inpath+"Bgo"));
  if(0 == fEvtBgo){
    fEvtBgo = new DmpEvtBgoRaw();
    gDataBuffer->LinkRootFile(inpath+"Bgo",fEvtBgo);
  }
  fEvtPsd = dynamic_cast<DmpEvtPsdRaw*>(gDataBuffer->ReadObject(inpath+"Psd"));
  if(0 == fEvtPsd){
    fEvtPsd = new DmpEvtPsdRaw();
    gDataBuffer->LinkRootFile(inpath+"Psd",fEvtPsd);
  }
  // create Hist map
  short layerNo = DmpParameterBgo::kPlaneNo*2;
  short barNo = DmpParameterBgo::kBarNo;
  for(short l=0;l<layerNo;++l){
    for(short b=0;b<barNo;++b){
          char name[50];
          short gid_bar = DmpBgoBase::ConstructGlobalBarID(l,b);
          snprintf(name,50,"BgoMip_%05d-L%02d_B%02d",gid_bar,l,b);
          fBgoMipHist.insert(std::make_pair(gid_bar,new TH1D(name,name,2000,0,800000)));
    }
  }

  layerNo = DmpParameterPsd::kPlaneNo*2;
  barNo = DmpParameterPsd::kStripNo;
  for(short l=0;l<layerNo;++l){
    for(short b=0;b<barNo;++b){
          char name[50];
          short gid_bar = DmpPsdBase::ConstructGlobalStripID(l,b);
          snprintf(name,50,"PsdPed_%05d-L%02d_S%02d",gid_bar,l,b);
          fPsdMipHist.insert(std::make_pair(gid_bar,new TH1D(name,name,1000,0,1000)));
  //std::cout<<"DEBUG: "<<__FILE__<<"("<<__LINE__<<")\t"<<name<<std::endl;
    }
  }

  return true;
}

//-------------------------------------------------------------------
bool DmpAlgCalMipsBar::ProcessThisEvent(){
  if(fEvtHeader->GetSecond() < gCore->GetStartTime() || fEvtHeader->GetSecond() > gCore->GetStopTime()){
    return false;
  }
  if(fEvtHeader->EnabledPeriodTrigger() && fEvtHeader->GeneratedPeriodTrigger()){
    return false;
  }
  // bgo Mip bar
  short layerNo = DmpParameterBgo::kPlaneNo*2;
  short barNo = DmpParameterBgo::kBarNo;
  for(short l=0;l<layerNo;++l){
    for(short b=0;b<barNo;++b){
      std::vector<short>::iterator s0dy8 = std::find(fEvtBgo->fGlobalDynodeID.begin(),fEvtBgo->fGlobalDynodeID.end(),DmpBgoBase::ConstructGlobalDynodeID(l,b,0,8));
      std::vector<short>::iterator s1dy8 = std::find(fEvtBgo->fGlobalDynodeID.begin(),fEvtBgo->fGlobalDynodeID.end(),DmpBgoBase::ConstructGlobalDynodeID(l,b,1,8));
      if(s0dy8 == fEvtBgo->fGlobalDynodeID.end() || s1dy8 == fEvtBgo->fGlobalDynodeID.end()){
        continue;
      }
      int i0 = s0dy8 - fEvtBgo->fGlobalDynodeID.begin();
      int i1 = s1dy8 - fEvtBgo->fGlobalDynodeID.begin();
      if(fEvtBgo->fADC.at(i0) > 1500 || fEvtBgo->fADC.at(i1) > 1500){
        continue;
      }
      fBgoMipHist[DmpBgoBase::ConstructGlobalBarID(l,b)]->Fill(fEvtBgo->fADC.at(i0) * fEvtBgo->fADC.at(i1));
    }
  }
  // Psd Mip bar
  layerNo = DmpParameterPsd::kPlaneNo*2;
  barNo = DmpParameterPsd::kStripNo;

  if(fFirstEvtTime == -1){
    fFirstEvtTime = fEvtHeader->GetSecond();
  }
  fLastEvtTime = fEvtHeader->GetSecond();
  return true;
}

//-------------------------------------------------------------------
bool DmpAlgCalMipsBar::Finalize(){
  TF1 *lxg_f = gMyFunctions->GetLanGau();
  std::string histFileName = gRootIOSvc->GetOutputPath()+gRootIOSvc->GetInputStem()+"_MipsBarHist.root";
  TFile *histFile = new TFile(histFileName.c_str(),"RECREATE");
  histFile->mkdir("Bgo");
  histFile->mkdir("Psd");

  // create output txtfile      BGO
  histFile->cd("Bgo");
  std::string name = "BgoBar_"+gRootIOSvc->GetInputStem()+".mip";
  OutBgoMipData.open(name.c_str(),std::ios::out);
  OutBgoMipData<<gRootIOSvc->GetInputFileName()<<std::endl;
  OutBgoMipData<<DmpTimeConvertor::Second2Date(fFirstEvtTime)<<std::endl;
  OutBgoMipData<<DmpTimeConvertor::Second2Date(fLastEvtTime)<<std::endl;
  OutBgoMipData<<"globalBarID\tlayer\tbar\t\tWidth\tMP\tArea\tGSigma"<<std::endl;
  for(std::map<short,TH1D*>::iterator aHist=fBgoMipHist.begin();aHist!=fBgoMipHist.end();++aHist){
    // Fit and save output data
      double mean = aHist->second->GetMean(), sigma = aHist->second->GetRMS();
      lxg_f->SetRange(mean-2*sigma,mean+5*sigma);
      //aHist->second->Fit(lxg_f,"RQ");
      OutBgoMipData<<aHist->first<<"\t"<<DmpBgoBase::GetLayerID(aHist->first)<<"\t"<<DmpBgoBase::GetBarID(aHist->first)<<"\t";
      for(int ip=0;ip<lxg_f->GetNumberFreeParameters();++ip){
        OutBgoMipData<<"\t"<<lxg_f->GetParameter(ip);
      }
      OutBgoMipData<<std::endl;
      aHist->second->Write();
      delete aHist->second;
  }
  OutBgoMipData.close();

  // create output txtfile      PSD

  delete histFile;
  return true;
}


