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

#include "DmpEvtHeader.h"
#include "DmpEvtBgoRaw.h"
#include "DmpAlgCalibrationMips.h"
#include "DmpDataBuffer.h"
#include "DmpParameterBgo.h"
#include "DmpParameterPsd.h"
#include "DmpBgoBase.h"
#include "DmpPsdBase.h"
#include "DmpCore.h"
#include "DmpTimeConvertor.h"
#include "MyFunctions.h"

#define  Psd_Bgo_Gap 300

//-------------------------------------------------------------------
DmpAlgCalibrationMips::DmpAlgCalibrationMips()
 :DmpVAlg("DmpAlgCalibrationMips"),
  fEvtHeader(0),
  fEvtBgo(0),
  fEvtPsd(0),
  fFirstEvtTime(-1),
  fLastEvtTime(-1),
  fRange_lo(100),
  fRange_hi(1600),
  fBinNo(200)
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

  char name[50];
  short gid = -1;
  // Bgo create Hist map
  for(short l=0;l<DmpParameterBgo::kPlaneNo*2;++l){
    for(short b=0;b<DmpParameterBgo::kBarNo;++b){
      gid = DmpBgoBase::ConstructGlobalBarID(l,b);
      snprintf(name,50,"BgoMip_%05d-L%02d_B%02d",gid,l,b);
      fBgoMipsHist_Bar.insert(std::make_pair(gid,new TH1D(name,name,fBinNo,fRange_lo,fRange_hi)));
      for(short s=0;s<DmpParameterBgo::kSideNo;++s){
        gid = DmpBgoBase::ConstructGlobalDynodeID(l,b,s,8);
        snprintf(name,50,"BgoMip_%05d-L%02d_B%02d_Dy%02d",gid,l,b,s*10+8);
        fBgoMipsHist_Dy.insert(std::make_pair(gid,new TH1D(name,name,fBinNo,fRange_lo,fRange_hi)));
      }
    }
  }
  // Psd
  for(short l=0;l<DmpParameterPsd::kPlaneNo*2;++l){
    for(short b=0;b<DmpParameterPsd::kStripNo;++b){
      gid = DmpPsdBase::ConstructGlobalStripID(l,b);
      snprintf(name,50,"PsdMip_%05d-L%02d_S%02d",gid,l,b);
      fPsdMipsHist_Bar.insert(std::make_pair(gid,new TH1D(name,name,fBinNo,fRange_lo,fRange_hi-Psd_Bgo_Gap)));
      for(short s=0;s<DmpParameterPsd::kSideNo;++s){
        gid = DmpPsdBase::ConstructGlobalDynodeID(l,b,s,8);
        snprintf(name,50,"PsdMip_%05d-L%02d_S%02d_Dy%02d",gid,l,b,s*10+8);
        fPsdMipsHist_Dy.insert(std::make_pair(gid,new TH1D(name,name,fBinNo,fRange_lo,fRange_hi-Psd_Bgo_Gap)));
      }
    }
  }

  return true;
}

//-------------------------------------------------------------------
bool DmpAlgCalibrationMips::ProcessThisEvent(){
  if(fEvtHeader->GetSecond() < gCore->GetStartTime() || fEvtHeader->GetSecond() > gCore->GetStopTime()){
    return false;
  }
  if(fEvtHeader->EnabledPeriodTrigger() && fEvtHeader->GeneratedPeriodTrigger()){
    return false;
  }
  if(fFirstEvtTime == -1){
    fFirstEvtTime = fEvtHeader->GetSecond();
  }
  fLastEvtTime = fEvtHeader->GetSecond();

  //-------------------------------------------------------------------
  // bgo Mip bar
  static short bgo_layerNo = DmpParameterBgo::kPlaneNo*2;
  for(short l=0;l<bgo_layerNo;++l){
    for(short b=0;b<DmpParameterBgo::kBarNo;++b){
      short gid_dy_s0 = DmpBgoBase::ConstructGlobalDynodeID(l,b,0,8);
      short gid_dy_s1 = DmpBgoBase::ConstructGlobalDynodeID(l,b,1,8);
      std::vector<short>::iterator s0dy8 = std::find(fEvtBgo->fGlobalDynodeID.begin(),fEvtBgo->fGlobalDynodeID.end(),gid_dy_s0);
      std::vector<short>::iterator s1dy8 = std::find(fEvtBgo->fGlobalDynodeID.begin(),fEvtBgo->fGlobalDynodeID.end(),gid_dy_s1);
      if(s0dy8 == fEvtBgo->fGlobalDynodeID.end() || s1dy8 == fEvtBgo->fGlobalDynodeID.end()){
        continue;
      }
      int i0 = s0dy8 - fEvtBgo->fGlobalDynodeID.begin();
      int i1 = s1dy8 - fEvtBgo->fGlobalDynodeID.begin();
      if(fEvtBgo->fADC.at(i0) < fRange_lo || fEvtBgo->fADC.at(i0) > fRange_hi || fEvtBgo->fADC.at(i1) < fRange_lo || fEvtBgo->fADC.at(i1) > fRange_hi){
        continue;
      }
      fBgoMipsHist_Dy[gid_dy_s0]->Fill(fEvtBgo->fADC.at(i0));
      fBgoMipsHist_Dy[gid_dy_s1]->Fill(fEvtBgo->fADC.at(i1));
      fBgoMipsHist_Bar[DmpBgoBase::ConstructGlobalBarID(l,b)]->Fill(TMath::Sqrt(fEvtBgo->fADC.at(i0) * fEvtBgo->fADC.at(i1)));
    }
  }

  // Psd Mip bar
  static short psd_layerNo = DmpParameterPsd::kPlaneNo*2;
  for(short l=0;l<psd_layerNo;++l){
    for(short b=0;b<DmpParameterPsd::kStripNo;++b){
      short gid_dy_s0 = DmpPsdBase::ConstructGlobalDynodeID(l,b,0,8);
      short gid_dy_s1 = DmpPsdBase::ConstructGlobalDynodeID(l,b,1,8);
      std::vector<short>::iterator s0dy8 = std::find(fEvtPsd->fGlobalDynodeID.begin(),fEvtPsd->fGlobalDynodeID.end(),gid_dy_s0);
      std::vector<short>::iterator s1dy8 = std::find(fEvtPsd->fGlobalDynodeID.begin(),fEvtPsd->fGlobalDynodeID.end(),gid_dy_s1);
      if(s0dy8 == fEvtPsd->fGlobalDynodeID.end() || s1dy8 == fEvtPsd->fGlobalDynodeID.end()){
        continue;
      }
      int i0 = s0dy8 - fEvtPsd->fGlobalDynodeID.begin();
      int i1 = s1dy8 - fEvtPsd->fGlobalDynodeID.begin();
      if(fEvtPsd->fADC.at(i0) < fRange_lo || fEvtPsd->fADC.at(i0) > fRange_hi-Psd_Bgo_Gap || fEvtPsd->fADC.at(i1) < fRange_lo || fEvtPsd->fADC.at(i1) > fRange_hi-Psd_Bgo_Gap){
        continue;
      }
      fPsdMipsHist_Dy[gid_dy_s0]->Fill(fEvtPsd->fADC.at(i0));
      fPsdMipsHist_Dy[gid_dy_s1]->Fill(fEvtPsd->fADC.at(i1));
      fPsdMipsHist_Bar[DmpPsdBase::ConstructGlobalStripID(l,b)]->Fill(TMath::Sqrt(fEvtPsd->fADC.at(i0) * fEvtPsd->fADC.at(i1)));
    }
  }

  return true;
}

//-------------------------------------------------------------------
bool DmpAlgCalibrationMips::Finalize(){
  TF1 *lxg_f = gMyFunctions->GetLanGau();
  double sv[4]={20,500,5000,50};
  double pllo[4]={10,fRange_lo,1000,10};
  double plhi[4]={100,fRange_hi,100000,200};
  lxg_f->SetParameters(sv);
  for (int i=0; i<4; ++i) {
    lxg_f->SetParLimits(i, pllo[i], plhi[i]);
  }
  std::string histFileName = gRootIOSvc->GetOutputPath()+gRootIOSvc->GetInputStem()+"_MipsBarHist.root";
  TFile *histFile = gRootIOSvc->GetOutputRootFile();//new TFile(histFileName.c_str(),"RECREATE");

  // create output txtfile      BGO
  histFile->mkdir("Bgo");
  histFile->cd("Bgo");
  std::string name = "Bgo_"+gRootIOSvc->GetInputStem()+".mip1";
  o_MipData_BgoBar.open(name.c_str(),std::ios::out);
  o_MipData_BgoBar<<gRootIOSvc->GetInputFileName()<<std::endl;
  o_MipData_BgoBar<<DmpTimeConvertor::Second2Date(fFirstEvtTime)<<std::endl;
  o_MipData_BgoBar<<DmpTimeConvertor::Second2Date(fLastEvtTime)<<std::endl;
  o_MipData_BgoBar<<"globalBarID\tlayer\tbar\t\tWidth\tMP\tArea\tGSigma"<<std::endl;
  name = "Bgo_"+gRootIOSvc->GetInputStem()+".mip";
  o_MipData_BgoDy.open(name.c_str(),std::ios::out);
  o_MipData_BgoDy<<gRootIOSvc->GetInputFileName()<<std::endl;
  o_MipData_BgoDy<<DmpTimeConvertor::Second2Date(fFirstEvtTime)<<std::endl;
  o_MipData_BgoDy<<DmpTimeConvertor::Second2Date(fLastEvtTime)<<std::endl;
  o_MipData_BgoDy<<"globalDyID\tlayer\tbar\tside\t\tWidth\tMP\tArea\tGSigma"<<std::endl;
  lxg_f->SetRange(100,1500);
  short layerNo = DmpParameterBgo::kPlaneNo*2;
  short gid_bar = -1;
  short gid_dy = -1;
  for(short l=0;l<layerNo;++l){
    for(short b=0;b<DmpParameterBgo::kBarNo;++b){
      gid_bar = DmpBgoBase::ConstructGlobalBarID(l,b);
      o_MipData_BgoBar<<gid_bar<<"\t"<<l<<"\t"<<b<<"\t";
      fBgoMipsHist_Bar[gid_bar]->Fit(lxg_f,"RQB");
      for(int ip=0;ip<lxg_f->GetNumberFreeParameters();++ip){
        o_MipData_BgoBar<<"\t"<<lxg_f->GetParameter(ip);
      }
      o_MipData_BgoBar<<std::endl;
      fBgoMipsHist_Bar[gid_bar]->Write();
      delete fBgoMipsHist_Bar[gid_bar];

      for(short s = 0;s<DmpParameterBgo::kSideNo;++s){
        gid_dy = DmpBgoBase::ConstructGlobalDynodeID(l,b,s,8);
        o_MipData_BgoDy<<gid_dy<<"\t"<<l<<"\t"<<b<<"\t"<<s<<"\t";
        fBgoMipsHist_Dy[gid_dy]->Fit(lxg_f,"RQB");
        for(int ip=0;ip<lxg_f->GetNumberFreeParameters();++ip){
          o_MipData_BgoDy<<"\t"<<lxg_f->GetParameter(ip);
        }
        o_MipData_BgoDy<<std::endl;
        fBgoMipsHist_Dy[gid_dy]->Write();
        delete fBgoMipsHist_Dy[gid_dy];
      }
    }
  }
  o_MipData_BgoBar.close();
  o_MipData_BgoDy.close();

  // create output txtfile      PSD
  histFile->mkdir("Psd");
  histFile->cd("Psd");
  name = "Psd_"+gRootIOSvc->GetInputStem()+".mip1";
  o_MipData_PsdBar.open(name.c_str(),std::ios::out);
  o_MipData_PsdBar<<gRootIOSvc->GetInputFileName()<<std::endl;
  o_MipData_PsdBar<<DmpTimeConvertor::Second2Date(fFirstEvtTime)<<std::endl;
  o_MipData_PsdBar<<DmpTimeConvertor::Second2Date(fLastEvtTime)<<std::endl;
  o_MipData_PsdBar<<"globalStripID\tlayer\tstrip\t\tWidth\tMP\tArea\tGSigma"<<std::endl;
  name = "Psd_"+gRootIOSvc->GetInputStem()+".mip";
  o_MipData_PsdDy.open(name.c_str(),std::ios::out);
  o_MipData_PsdDy<<gRootIOSvc->GetInputFileName()<<std::endl;
  o_MipData_PsdDy<<DmpTimeConvertor::Second2Date(fFirstEvtTime)<<std::endl;
  o_MipData_PsdDy<<DmpTimeConvertor::Second2Date(fLastEvtTime)<<std::endl;
  o_MipData_PsdDy<<"globalDynodeID\tlayer\tstrip\tside\t\tWidth\tMP\tArea\tGSigma"<<std::endl;
  lxg_f->SetRange(100,1000);
  layerNo = DmpParameterPsd::kPlaneNo*2;
  for(short l=0;l<layerNo;++l){
    for(short b=0;b<DmpParameterPsd::kStripNo;++b){
      gid_bar = DmpPsdBase::ConstructGlobalStripID(l,b);
      o_MipData_PsdBar<<gid_bar<<"\t"<<l<<"\t"<<b<<"\t";
      fPsdMipsHist_Bar[gid_bar]->Fit(lxg_f,"RQB");
      for(int ip=0;ip<lxg_f->GetNumberFreeParameters();++ip){
        o_MipData_PsdBar<<"\t"<<lxg_f->GetParameter(ip);
      }
      o_MipData_PsdBar<<std::endl;
      fPsdMipsHist_Bar[gid_bar]->Write();
      delete fPsdMipsHist_Bar[gid_bar];

      for(short s = 0;s<DmpParameterPsd::kSideNo;++s){
        gid_dy = DmpPsdBase::ConstructGlobalDynodeID(l,b,s,8);
        o_MipData_PsdDy<<gid_dy<<"\t"<<l<<"\t"<<b<<"\t"<<s<<"\t";
        fPsdMipsHist_Dy[gid_dy]->Fit(lxg_f,"RQB");
        for(int ip=0;ip<lxg_f->GetNumberFreeParameters();++ip){
          o_MipData_PsdDy<<"\t"<<lxg_f->GetParameter(ip);
        }
        o_MipData_PsdDy<<std::endl;
        fPsdMipsHist_Dy[gid_dy]->Write();
        delete fPsdMipsHist_Dy[gid_dy];
      }
    }
  }
  o_MipData_PsdBar.close();
  o_MipData_PsdDy.close();

  delete histFile;
  return true;
}


