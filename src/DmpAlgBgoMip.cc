/*
 *  $Id: DmpAlgBgoMip.cc, 2014-09-11 15:46:28 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 11/09/2014
*/

#include <stdio.h>

#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"

#include "DmpEvtHeader.h"
#include "DmpEvtBgoRaw.h"
#include "DmpEvtBgoMip.h"
#include "DmpAlgBgoMip.h"
#include "DmpDataBuffer.h"
#include "DmpParameterBgo.h"
#include "DmpBgoBase.h"

//-------------------------------------------------------------------
DmpAlgBgoMip::DmpAlgBgoMip()
 :DmpVAlg("Cal/Bgo/Mip"),
  fEvtHeader(0),
  fBgoRaw(0),
  fBgoMip(0)
{
}

//-------------------------------------------------------------------
DmpAlgBgoMip::~DmpAlgBgoMip(){
}

//-------------------------------------------------------------------
bool DmpAlgBgoMip::Initialize(){
  gRootIOSvc->Set("OutData/FileName","./"+gRootIOSvc->GetInputFileName()+"_mip.root");
  // read input data
  fEvtHeader = new DmpEvtHeader();
  if(not gDataBuffer->ReadObject("Event/Rdc/EventHeader",fEvtHeader)){
    return false;
  }
  fBgoRaw = new DmpEvtBgoRaw();
  if(not gDataBuffer->ReadObject("Event/Rdc/Bgo",fBgoRaw)){
    return false;
  }
  // create output data holder
  fBgoMip = new DmpEvtBgoMip();
  if(not gDataBuffer->RegisterObject("Calibration/Bgo/Mip",fBgoMip,"DmpEvtBgoMip")){
    return false;
  }
  fBgoMip->UsedFileName = gRootIOSvc->GetInputFileName();
  gRootIOSvc->PrepareEvent(0);
  fBgoMip->StartTime = fEvtHeader->GetSecond();
  // create Hist map
  short layerNo = DmpParameterBgo::kPlaneNo*2;
  short barNo = DmpParameterBgo::kBarNo;
  for(short l=0;l<layerNo;++l){
    for(short b=0;b<barNo;++b){
      for(short s=0;s<DmpParameterBgo::kSideNo;++s){
        char name[50];
        short gid_dy = DmpBgoBase::ConstructGlobalDynodeID(l,b,s,8);
        snprintf(name,50,"BgoMip_%05d-L%02d_B%02d_Dy%02d",gid_dy,l,b,s*10+8);
        fMipHist.insert(std::make_pair(gid_dy,new TH1F(name,name,1000,-500,1500)));
      }
    }
  }
  return true;
}

//-------------------------------------------------------------------
bool DmpAlgBgoMip::ProcessThisEvent(){
  short nSignal = fBgoRaw->GetSignalSize();
  //std::map<short,>
  short gid = 0,adc = -999,l=-1,b=-1,s=-1;
  short adcTmp[DmpParameterBgo::kPlaneNo*2][DmpParameterBgo::kSideNo][DmpParameterBgo::kBarNo];
  for(short i=0;i<nSignal;++i){
    if(fBgoRaw->GetSignal(i,gid,adc)){
      if(DmpBgoBase::GetDynodeID(gid) == 8){
        adcTmp[DmpBgoBase::GetLayerID][DmpBgoBase::GetSideID][DmpBgoBase::GetBarID] = adc;
        fMipHist[gid]->Fill(adc);
      }
    }
  }
  return true;
}

//-------------------------------------------------------------------
bool DmpAlgBgoMip::Finalize(){
  TF1 *gausFit = new TF1("GausFit","gaus",-500,1500);
  std::string histFileName = gRootIOSvc->GetInputFileName()+"_mip_Hist.root";
  TFile *histFile = new TFile(histFileName.c_str(),"RECREATE");
  fBgoMip->StopTime = fEvtHeader->GetSecond();
  for(std::map<short,TH1F*>::iterator aHist=fMipHist.begin();aHist!=fMipHist.end();++aHist){
      fBgoMip->GlobalDynodeID.push_back(aHist->first);
// *
// *  TODO: fit and save output data 
// *
      float mean = aHist->second->GetMean(), sigma = aHist->second->GetRMS();
      for(short i = 0;i<3;++i){
        gausFit->SetRange(mean-2*sigma,mean+2*sigma);
        aHist->second->Fit(gausFit,"RQ");
        mean = gausFit->GetParameter(1);
        sigma = gausFit->GetParameter(2);
      }
      fBgoMip->Mean.push_back(mean);
      fBgoMip->Sigma.push_back(sigma);
      if((mean > 1600 || mean<-1600) && sigma >300){
         DmpLogWarning<<"GID = "<<aHist->first<<"\tmean = "<<mean<<"\tsigma = "<<sigma<<DmpLogEndl;
      }
      aHist->second->Write();
      delete aHist->second;
  }
  delete histFile;
  return true;
}


