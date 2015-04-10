/*
 *  $Id: DmpAlgCalibrationMips.cc, 2015-03-03 23:17:07 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 19/07/2014
*/

//#include <stdio.h>
//#include <boost/lexical_cast.hpp>
#include <algorithm>

#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TFile.h"

#include "DmpEvtBgoRaw.h"
#include "DmpAlgCalibrationMips.h"
#include "DmpDataBuffer.h"
#include "DmpParameterBgo.h"
#include "DmpParameterPsd.h"
#include "DmpLoadParameters.h"
#include "DmpBgoBase.h"
#include "DmpPsdBase.h"
#include "DmpCore.h"
#include "DmpTimeConvertor.h"
#include "MyFunctions.h"

#define CertainMip 10

//-------------------------------------------------------------------
DmpAlgCalibrationMips::DmpAlgCalibrationMips()
 :DmpVAlg("DmpAlgCalibrationMips"),
  fEvtBgo(0),
  fEvtPsd(0)
{
  gRootIOSvc->SetOutFileKey("CalMip");
}

//-------------------------------------------------------------------
DmpAlgCalibrationMips::~DmpAlgCalibrationMips(){
}

//-------------------------------------------------------------------
bool DmpAlgCalibrationMips::Initialize(){
  // read input data
  std::string inpath = "Event/Rec0/";
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

  char name[50];
  short gid = -1;
  // Bgo create Hist map
  for(short l=0;l<DmpParameterBgo::kPlaneNo*2;++l){
    for(short b=0;b<DmpParameterBgo::kBarNo;++b){
      gid = DmpBgoBase::ConstructGlobalBarID(l,b);
      snprintf(name,50,"BgoMip_%05d-L%02d_B%02d",gid,l,b);
      fBgoMipsHist_Bar.insert(std::make_pair(gid,new TH1D(name,name,500,0,2000)));
      for(short s=0;s<DmpParameterBgo::kSideNo;++s){
        gid = DmpBgoBase::ConstructGlobalDynodeID(l,b,s,8);
        snprintf(name,50,"BgoMip_%05d-L%02d_B%02d_Dy%02d",gid,l,b,s*10+8);
        fBgoMipsHist_Dy.insert(std::make_pair(gid,new TH1D(name,name,500,0,2000)));
      }
    }
  }
  // Psd
  for(short l=0;l<DmpParameterPsd::kPlaneNo*2;++l){
    for(short b=0;b<DmpParameterPsd::kStripNo;++b){
      gid = DmpPsdBase::ConstructGlobalStripID(l,b);
      snprintf(name,50,"PsdMip_%05d-L%02d_S%02d",gid,l,b);
      fPsdMipsHist_Bar.insert(std::make_pair(gid,new TH1D(name,name,500,0,2000)));
      for(short s=0;s<DmpParameterPsd::kSideNo;++s){
        gid = DmpPsdBase::ConstructGlobalDynodeID(l,b,s,8);
        snprintf(name,50,"PsdMip_%05d-L%02d_S%02d_Dy%02d",gid,l,b,s*10+8);
        fPsdMipsHist_Dy.insert(std::make_pair(gid,new TH1D(name,name,500,0,2000)));
      }
    }
  }

  return true;
}

//-------------------------------------------------------------------
bool DmpAlgCalibrationMips::ProcessThisEvent(){
  if(gCore->GetEventHeader()->EnabledPeriodTrigger() && gCore->GetEventHeader()->GeneratedPeriodTrigger()){
    return false;
  }
  //-------------------------------------------------------------------
  // bgo Mip bar
  static std::vector<short>   barIDs;
  static short dyid0,dyid1;
  for(short l=0;l<DmpParameterBgo::kPlaneNo*2;++l){
    GetBarIDOfLayer_bgo(barIDs,l);
    if(barIDs.size() > 3) return false;
  }
  for(short l=0;l<DmpParameterBgo::kPlaneNo*2;++l){
    GetBarIDOfLayer_bgo(barIDs,l);
    if(barIDs.size() > 2)continue;
    for(size_t ib = 0;ib<barIDs.size();++ib){
      dyid0 = DmpBgoBase::ConstructGlobalDynodeID(l,barIDs[ib],0,8);
      dyid1 = DmpBgoBase::ConstructGlobalDynodeID(l,barIDs[ib],1,8);
      fBgoMipsHist_Dy[dyid0]->Fill(fEvtBgo->fADC[dyid0]);
      fBgoMipsHist_Dy[dyid1]->Fill(fEvtBgo->fADC[dyid1]);
      fBgoMipsHist_Bar[DmpBgoBase::ConstructGlobalBarID(l,barIDs[ib])]->Fill(TMath::Sqrt(fEvtBgo->fADC[dyid0] * fEvtBgo->fADC[dyid1]));
    }
  }

  // Psd Mip bar
  for(short l=0;l<DmpParameterPsd::kPlaneNo*2;++l){
    GetBarIDOfLayer_psd(barIDs,l);
    if(barIDs.size() > 2)continue;
    for(size_t ib = 0;ib<barIDs.size();++ib){
      dyid0 = DmpPsdBase::ConstructGlobalDynodeID(l,barIDs[ib],0,8);
      dyid1 = DmpPsdBase::ConstructGlobalDynodeID(l,barIDs[ib],1,8);
      fPsdMipsHist_Dy[dyid0]->Fill(fEvtPsd->fADC[dyid0]);
      fPsdMipsHist_Dy[dyid1]->Fill(fEvtPsd->fADC[dyid1]);
      fPsdMipsHist_Bar[DmpPsdBase::ConstructGlobalStripID(l,barIDs[ib])]->Fill(TMath::Sqrt(fEvtPsd->fADC[dyid0] * fEvtPsd->fADC[dyid1]));
    }
  }

  return true;
}

//-------------------------------------------------------------------
bool DmpAlgCalibrationMips::Finalize(){
  TF1 *lxg_f = gMyFunctions->GetLanGau();
  double sv[4]={20,500,5000,50};
  double pllo[4]={10,0,1000,10};
  double plhi[4]={100,2000,100000,200};
  lxg_f->SetParameters(sv);
  for (int i=0; i<4; ++i) {
    lxg_f->SetParLimits(i, pllo[i], plhi[i]);
  }
  std::string histFileName = gRootIOSvc->GetOutputPath()+gRootIOSvc->GetInputStem()+"_MipsBarHist.root";
  TFile *histFile = gRootIOSvc->GetOutputRootFile();//new TFile(histFileName.c_str(),"RECREATE");

  // create output txtfile      BGO
  histFile->mkdir("Bgo");
  histFile->cd("Bgo");
  std::string name = "Bgo_"+gRootIOSvc->GetInputStem()+".mip";
  o_MipData_Bgo.open(name.c_str(),std::ios::out);
  o_MipData_Bgo<<Mark_S<<"\nFileName="<<gRootIOSvc->GetInputFileName()<<std::endl;
  o_MipData_Bgo<<"StartTime="<<gCore->GetTimeFirstOutput()<<"\nStopTime="<<gCore->GetTimeLastOutput()<<std::endl;
  o_MipData_Bgo<<Mark_D<<std::endl;
  lxg_f->SetRange(100,1500);
  short layerNo = DmpParameterBgo::kPlaneNo*2;
  short gid_bar = -1;
  short gid_dy = -1;
  double mpv = 0;
  for(short l=0;l<layerNo;++l){
    for(short b=0;b<DmpParameterBgo::kBarNo;++b){
      gid_bar = DmpBgoBase::ConstructGlobalBarID(l,b);
      o_MipData_Bgo<<gid_bar<<"\t\t"<<l<<"\t\t"<<b<<"\t\t";
      fBgoMipsHist_Bar[gid_bar]->Fit(lxg_f,"RQB");
      for(int ip=0;ip<lxg_f->GetNumberFreeParameters();++ip){
        o_MipData_Bgo<<"\t\t"<<lxg_f->GetParameter(ip);
      }
      mpv = lxg_f->GetParameter(1);
      o_MipData_Bgo<<"\t\t"<<lxg_f->GetMaximumX(0.8*mpv,1.2*mpv)<<"\t\t"<<lxg_f->GetChisquare() / lxg_f->GetNDF()<<"\t\t"<<fBgoMipsHist_Bar[gid_bar]->GetEntries()<<std::endl;
      fBgoMipsHist_Bar[gid_bar]->Write();
      delete fBgoMipsHist_Bar[gid_bar];

      for(short s = 0;s<DmpParameterBgo::kSideNo;++s){
        gid_dy = DmpBgoBase::ConstructGlobalDynodeID(l,b,s,8);
        o_MipData_Bgo<<gid_dy<<"\t\t"<<l<<"\t\t"<<b<<"\t\t"<<s;
        fBgoMipsHist_Dy[gid_dy]->Fit(lxg_f,"RQB");
        for(int ip=0;ip<lxg_f->GetNumberFreeParameters();++ip){
          o_MipData_Bgo<<"\t\t"<<lxg_f->GetParameter(ip);
        }
        mpv = lxg_f->GetParameter(1);
        o_MipData_Bgo<<"\t\t"<<lxg_f->GetMaximumX(0.8*mpv,1.2*mpv)<<"\t\t"<<lxg_f->GetChisquare() / lxg_f->GetNDF()<<"\t\t"<<fBgoMipsHist_Dy[gid_dy]->GetEntries()<<std::endl;
        fBgoMipsHist_Dy[gid_dy]->Write();
        delete fBgoMipsHist_Dy[gid_dy];
      }
    }
  }
  o_MipData_Bgo<<Mark_N<<std::endl;
  o_MipData_Bgo.close();

  // create output txtfile      PSD
  histFile->mkdir("Psd");
  histFile->cd("Psd");
  name= "Psd_"+gRootIOSvc->GetInputStem()+".mip";
  o_MipData_Psd.open(name.c_str(),std::ios::out);
  o_MipData_Psd<<Mark_S<<"\nFileName="<<gRootIOSvc->GetInputFileName()<<std::endl;
  o_MipData_Psd<<"StartTime="<<gCore->GetTimeFirstOutput()<<"\nStopTime="<<gCore->GetTimeLastOutput()<<std::endl;
  o_MipData_Psd<<Mark_D<<std::endl;
  lxg_f->SetRange(100,1000);
  layerNo = DmpParameterPsd::kPlaneNo*2;
  for(short l=0;l<layerNo;++l){
    for(short b=0;b<DmpParameterPsd::kStripNo;++b){
      gid_bar = DmpPsdBase::ConstructGlobalStripID(l,b);
      o_MipData_Psd<<gid_bar<<"\t\t"<<l<<"\t\t"<<b<<"\t\t";
      fPsdMipsHist_Bar[gid_bar]->Fit(lxg_f,"RQB");
      for(int ip=0;ip<lxg_f->GetNumberFreeParameters();++ip){
        o_MipData_Psd<<"\t\t"<<lxg_f->GetParameter(ip);
      }
      mpv = lxg_f->GetParameter(1);
      o_MipData_Psd<<"\t\t"<<lxg_f->GetMaximumX(0.8*mpv,1.2*mpv)<<"\t\t"<<lxg_f->GetChisquare() / lxg_f->GetNDF()<<"\t\t"<<fPsdMipsHist_Bar[gid_bar]->GetEntries()<<std::endl;
      fPsdMipsHist_Bar[gid_bar]->Write();
      delete fPsdMipsHist_Bar[gid_bar];

      for(short s = 0;s<DmpParameterPsd::kSideNo;++s){
        gid_dy = DmpPsdBase::ConstructGlobalDynodeID(l,b,s,8);
        o_MipData_Psd<<gid_dy<<"\t\t"<<l<<"\t\t"<<b<<"\t\t"<<s;
        fPsdMipsHist_Dy[gid_dy]->Fit(lxg_f,"RQB");
        for(int ip=0;ip<lxg_f->GetNumberFreeParameters();++ip){
          o_MipData_Psd<<"\t\t"<<lxg_f->GetParameter(ip);
        }
        mpv = lxg_f->GetParameter(1);
        o_MipData_Psd<<"\t\t"<<lxg_f->GetMaximumX(0.8*mpv,1.2*mpv)<<"\t\t"<<lxg_f->GetChisquare() / lxg_f->GetNDF()<<"\t\t"<<fPsdMipsHist_Dy[gid_dy]->GetEntries()<<std::endl;
        fPsdMipsHist_Dy[gid_dy]->Write();
        delete fPsdMipsHist_Dy[gid_dy];
      }
    }
  }
  o_MipData_Psd<<Mark_N<<std::endl;
  o_MipData_Psd.close();

  return true;
}

void DmpAlgCalibrationMips::GetBarIDOfLayer_bgo(std::vector<short> &ret,short l,short side)const
{
  ret.clear();
  static short gid_dy =0;
  static short gid_dy5 =0;
  if(side == 0  || side == 1){
    for(short b=0;b<DmpParameterBgo::kBarNo;++b){
      gid_dy = DmpBgoBase::ConstructGlobalDynodeID(l,b,side,8);
      if(fEvtBgo->fADC.find(gid_dy) == fEvtBgo->fADC.end()){
        continue;
      }
      gid_dy5 = DmpBgoBase::ConstructGlobalDynodeID(l,b,side,5);
      if(fEvtBgo->fADC.find(gid_dy5) != fEvtBgo->fADC.end()){
        if(fEvtBgo->fADC[gid_dy5] > 30){
std::cout<<"DEBUG: dy8 error\t"<<__FILE__<<"("<<__LINE__<<")"<<std::endl;
          continue;
        }
      }
      if(fEvtBgo->fADC[gid_dy] > CertainMip){
        ret.push_back(b);
      }
    }
  }else{
    std::vector<short> bs0,bs1;
    this->GetBarIDOfLayer_bgo(bs0,l,0);
    this->GetBarIDOfLayer_bgo(bs1,l,1);
    for(short i=0;i<bs0.size();++i){
      if(std::find(bs1.begin(),bs1.end(),bs0[i]) != bs1.end()){
        ret.push_back(bs0[i]);
      }
    }
  }
}

void DmpAlgCalibrationMips::GetBarIDOfLayer_psd(std::vector<short> &ret,short l,short side)const
{
  ret.clear();
  static short gid_dy =0;
  if(side == 0  || side == 1){
    for(short b=0;b<DmpParameterPsd::kStripNo;++b){
      gid_dy = DmpPsdBase::ConstructGlobalDynodeID(l,b,side,8);
      if(fEvtPsd->fADC.find(gid_dy) == fEvtPsd->fADC.end()){
        continue;
      }
      if(fEvtPsd->fADC[gid_dy] > CertainMip){
        ret.push_back(b);
      }
    }
  }else{
    std::vector<short> bs0,bs1;
    this->GetBarIDOfLayer_psd(bs0,l,0);
    this->GetBarIDOfLayer_psd(bs1,l,1);
    for(short i=0;i<bs0.size();++i){
      if(std::find(bs1.begin(),bs1.end(),bs0[i]) != bs1.end()){
        ret.push_back(bs0[i]);
      }
    }
  }
}

