/*
 *  $Id: DmpAlgBgoMip.cc, 2014-09-19 13:09:57 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 11/09/2014
*/

#include <stdio.h>

#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"

#include "DmpEvtHeader.h"
#include "DmpEvtBgoRaw.h"
#include "DmpEvtBgoMip.h"
#include "DmpAlgBgoMip.h"
#include "DmpDataBuffer.h"
#include "DmpParameterBgo.h"
#include "DmpBgoBase.h"
#include "DmpCore.h"

//-------------------------------------------------------------------
DmpAlgBgoMip::DmpAlgBgoMip()
 :DmpVAlg("Cal/Bgo/Mip"),
  fEvtHeader(0),
  fBgoRaw(0),
  fBgoMip(0)
{
  gRootIOSvc->Set("Output/Key","mip");
}

//-------------------------------------------------------------------
DmpAlgBgoMip::~DmpAlgBgoMip(){
}

//-------------------------------------------------------------------
bool DmpAlgBgoMip::Initialize(){
  // read input data
  fEvtHeader = new DmpEvtHeader();
  if(not gDataBuffer->ReadObject("Event/Rdc/EventHeader",fEvtHeader)){
    return false;
  }
  fBgoRaw = new DmpEvtBgoRaw();
  if(not gDataBuffer->ReadObject("Event/Raw/Bgo",fBgoRaw)){
    return false;
  }
  // create output data holder
  fBgoMip = new DmpEvtBgoMip();
  if(not gDataBuffer->RegisterObject("Calibration/Bgo/Mip",fBgoMip,"DmpEvtBgoMip")){
    return false;
  }
  fBgoMip->UsedFileName = gRootIOSvc->GetInputFileName();
  gRootIOSvc->PrepareEvent(gCore->GetCurrentEventID());
  fBgoMip->StartTime = fEvtHeader->fSecond;
  // create Hist map
  short layerNo = DmpParameterBgo::kPlaneNo*2;
  short barNo = DmpParameterBgo::kBarNo;
  for(short l=0;l<layerNo;++l){
    for(short b=0;b<barNo;++b){
      for(short s=0;s<DmpParameterBgo::kSideNo;++s){
        char name[50];
        short gid_dy = DmpBgoBase::ConstructGlobalDynodeID(l,b,s,8);
        snprintf(name,50,"BgoMip_%05d-L%02d_B%02d_Dy%02d",gid_dy,l,b,s*10+8);
        fMipHist.insert(std::make_pair(gid_dy,new TH1D(name,name,150,-500,2500)));
      }
    }
  }
  return true;
}

//-------------------------------------------------------------------
bool DmpAlgBgoMip::ProcessThisEvent(){
  short max_adc[DmpParameterBgo::kPlaneNo*2][DmpParameterBgo::kSideNo];
  short barID_max_adc[DmpParameterBgo::kPlaneNo*2][DmpParameterBgo::kSideNo];
  for(short layer=0;layer<DmpParameterBgo::kPlaneNo*2;++layer){
    for(short side =0;side<DmpParameterBgo::kSideNo;++side){
      max_adc[layer][side] = -9999;
      barID_max_adc[layer][side] = DmpParameterBgo::kBarNo; // id 0 ~ (kBarNo-1)
    }
  }
//-------------------------------------------------------------------
  short gid = 0,adc = -999,l=-1,b=-1,s=-1,d=-1;
  short nSignal = fBgoRaw->fADC.size();
  for(short i=0;i<nSignal;++i){
    gid = fBgoRaw->fGlobalDynodeID[i];
    adc = fBgoRaw->fADC[i];
    DmpBgoBase::LoadLBSDID(gid,l,b,s,d);
    if(b == DmpParameterBgo::kBarNo || b == DmpParameterBgo::kBarNo+1){
      continue;
    }
    if(d == 8){
      if(adc > max_adc[l][s]){
        barID_max_adc[l][s] = b;
        max_adc[l][s] = adc;
      }
    }
  }
//-------------------------------------------------------------------
/*
  short side_0_barID = barID_max_adc[0][0];
  short side_1_barID = barID_max_adc[0][1];
  for(short i=0;i<DmpParameterBgo::kPlaneNo;++i){
    if(side_0_barID != barID_max_adc[i][0]){
      return false;
    }
    if(side_0_barID != barID_max_adc[i][1]){
      return false;
    }
  }
  */
//-------------------------------------------------------------------
  for(short layer=0;layer<DmpParameterBgo::kPlaneNo*2;++layer){
    for(short side =0;side<DmpParameterBgo::kSideNo;++side){
      fMipHist[DmpBgoBase::ConstructGlobalDynodeID(layer,barID_max_adc[layer][side],side,8)]->Fill(max_adc[layer][side]);
    }
  }
  return true;
}

//-------------------------------------------------------------------
bool DmpAlgBgoMip::Finalize(){
  TF1 *gausFit = new TF1("GausFit","gaus",-100,2500);
  std::string histFileName = gRootIOSvc->GetOutputPath()+gRootIOSvc->GetOutputStem()+"_Hist.root";
  TFile *histFile = new TFile(histFileName.c_str(),"RECREATE");
  fBgoMip->StopTime = fEvtHeader->fSecond;
  for(std::map<short,TH1D*>::iterator aHist=fMipHist.begin();aHist!=fMipHist.end();++aHist){
      fBgoMip->GlobalDynodeID.push_back(aHist->first);
      double mean = aHist->second->GetMean(), sigma = aHist->second->GetRMS();
      for(short i = 0;i<3;++i){
        gausFit->SetRange(mean-1.5*sigma,mean+1.5*sigma);
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



