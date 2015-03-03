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
#include "TFile.h"
#include "TMath.h"

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

//-------------------------------------------------------------------
DmpAlgCalibrationMips::DmpAlgCalibrationMips()
 :DmpVAlg("Cal/Bgo/Ped"),
  fEvtHeader(0),
  fEvtBgo(0),
  fEvtPsd(0),
  fFirstEvtTime(-1),
  fLastEvtTime(-1),
  fEnableCaliBar(true),
  fEnableCaliDynode(true)
{
  gRootIOSvc->SetOutputKey("CalMip");
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
  // Bgo create Hist map
  if(fEnableCaliBar)
  {
  for(short l=0;l<DmpParameterBgo::kPlaneNo*2;++l){
    for(short b=0;b<DmpParameterBgo::kBarNo;++b){
          char name[50];
          short gid_bar = DmpBgoBase::ConstructGlobalBarID(l,b);
          snprintf(name,50,"BgoMip_%05d-L%02d_B%02d",gid_bar,l,b);
          fBgoMipsHist_Bar.insert(std::make_pair(gid_bar,new TH1D(name,name,750,0,1500)));
    }
  }
  }
  if(fEnableCaliDynode)
  {
  for(short l=0;l<DmpParameterBgo::kPlaneNo*2;++l){
    for(short b=0;b<DmpParameterBgo::kBarNo;++b){
      for(short s=0;s<DmpParameterBgo::kSideNo;++s){
          char name[50];
          short gid_dy = DmpBgoBase::ConstructGlobalDynodeID(l,b,s,8);
          snprintf(name,50,"BgoMip_%05d-L%02d_B%02d_Dy%02d",gid_dy,l,b,s*10+8);
          fBgoMipsHist_Dy.insert(std::make_pair(gid_dy,new TH1D(name,name,750,0,1500)));
      }
    }
  }
  }
    // Psd
  if(fEnableCaliBar)
  {
  for(short l=0;l<DmpParameterPsd::kPlaneNo*2;++l){
    for(short b=0;b<DmpParameterPsd::kStripNo;++b){
          char name[50];
          short gid_bar = DmpPsdBase::ConstructGlobalStripID(l,b);
          snprintf(name,50,"PsdMip_%05d-L%02d_S%02d",gid_bar,l,b);
          fPsdMipsHist_Bar.insert(std::make_pair(gid_bar,new TH1D(name,name,600,0,1200)));
    }
  }
  }
  if(fEnableCaliDynode)
  {
  for(short l=0;l<DmpParameterPsd::kPlaneNo*2;++l){
    for(short b=0;b<DmpParameterPsd::kStripNo;++b){
      for(short s=0;s<DmpParameterPsd::kSideNo;++s){
          char name[50];
          short gid_dy = DmpPsdBase::ConstructGlobalDynodeID(l,b,s,8);
          snprintf(name,50,"PsdMip_%05d-L%02d_S%02d_Dy%02d",gid_dy,l,b,s*10+8);
          fPsdMipsHist_Dy.insert(std::make_pair(gid_dy,new TH1D(name,name,600,0,1200)));
      }
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

  if(fEnableCaliBar)    {_calBarOrStrip();}
  if(fEnableCaliDynode) {_calDynode();}

  if(fFirstEvtTime == -1){
    fFirstEvtTime = fEvtHeader->GetSecond();
  }
  fLastEvtTime = fEvtHeader->GetSecond();
  return true;
}

//-------------------------------------------------------------------
bool DmpAlgCalibrationMips::Finalize(){
  TF1 *lxg_f = gMyFunctions->GetLanGau();
  std::string histFileName = gRootIOSvc->GetOutputPath()+gRootIOSvc->GetInputStem()+"_MipsBarHist.root";
  TFile *histFile = new TFile(histFileName.c_str(),"RECREATE");
  histFile->mkdir("Bgo");
  histFile->mkdir("Psd");

  if(fEnableCaliBar)
  {
  // create output txtfile      BGO
  histFile->cd("Bgo");
  std::string name = "Bgo_"+gRootIOSvc->GetInputStem()+".mip1";
  o_MipData_BgoBar.open(name.c_str(),std::ios::out);
  o_MipData_BgoBar<<gRootIOSvc->GetInputFileName()<<std::endl;
  o_MipData_BgoBar<<DmpTimeConvertor::Second2Date(fFirstEvtTime)<<std::endl;
  o_MipData_BgoBar<<DmpTimeConvertor::Second2Date(fLastEvtTime)<<std::endl;
  o_MipData_BgoBar<<"globalBarID\tlayer\tbar\t\tWidth\tMP\tArea\tGSigma"<<std::endl;
  lxg_f->SetRange(80,2000);
  for(std::map<short,TH1D*>::iterator aHist=fBgoMipsHist_Bar.begin();aHist!=fBgoMipsHist_Bar.end();++aHist){
    // Fit and save output data
      double mean = aHist->second->GetMean(), sigma = aHist->second->GetRMS();
      //aHist->second->Fit(lxg_f,"Q");
      o_MipData_BgoBar<<aHist->first<<"\t"<<DmpBgoBase::GetLayerID(aHist->first)<<"\t"<<DmpBgoBase::GetBarID(aHist->first)<<"\t";
      for(int ip=0;ip<lxg_f->GetNumberFreeParameters();++ip){
        o_MipData_BgoBar<<"\t"<<lxg_f->GetParameter(ip);
      }
      o_MipData_BgoBar<<std::endl;
      aHist->second->Write();
      delete aHist->second;
  }
  o_MipData_BgoBar.close();

  // create output txtfile      PSD
  histFile->cd("Psd");
  name = "Psd_"+gRootIOSvc->GetInputStem()+".mip1";
  o_MipData_PsdBar.open(name.c_str(),std::ios::out);
  o_MipData_PsdBar<<gRootIOSvc->GetInputFileName()<<std::endl;
  o_MipData_PsdBar<<DmpTimeConvertor::Second2Date(fFirstEvtTime)<<std::endl;
  o_MipData_PsdBar<<DmpTimeConvertor::Second2Date(fLastEvtTime)<<std::endl;
  o_MipData_PsdBar<<"globalStripID\tlayer\tstrip\t\tWidth\tMP\tArea\tGSigma"<<std::endl;
  lxg_f->SetRange(50,1200);
  for(std::map<short,TH1D*>::iterator aHist=fPsdMipsHist_Bar.begin();aHist!=fPsdMipsHist_Bar.end();++aHist){
    // Fit and save output data
      double mean = aHist->second->GetMean(), sigma = aHist->second->GetRMS();
      //aHist->second->Fit(lxg_f,"Q");
      o_MipData_PsdBar<<aHist->first<<"\t"<<DmpPsdBase::GetLayerID(aHist->first)<<"\t"<<DmpPsdBase::GetStripID(aHist->first)<<"\t";
      for(int ip=0;ip<lxg_f->GetNumberFreeParameters();++ip){
        o_MipData_PsdBar<<"\t"<<lxg_f->GetParameter(ip);
      }
      o_MipData_PsdBar<<std::endl;
      aHist->second->Write();
      delete aHist->second;
  }
  o_MipData_PsdBar.close();
  }

  if(fEnableCaliDynode)
  {
  // create output txtfile      BGO
  histFile->cd("Bgo");
  std::string name = "Bgo_"+gRootIOSvc->GetInputStem()+".mip";
  o_MipData_BgoDy.open(name.c_str(),std::ios::out);
  o_MipData_BgoDy<<gRootIOSvc->GetInputFileName()<<std::endl;
  o_MipData_BgoDy<<DmpTimeConvertor::Second2Date(fFirstEvtTime)<<std::endl;
  o_MipData_BgoDy<<DmpTimeConvertor::Second2Date(fLastEvtTime)<<std::endl;
  o_MipData_BgoDy<<"globalDyID\tlayer\tbar\tside\t\tWidth\tMP\tArea\tGSigma"<<std::endl;
  lxg_f->SetRange(80,2000);
  for(std::map<short,TH1D*>::iterator aHist=fBgoMipsHist_Dy.begin();aHist!=fBgoMipsHist_Dy.end();++aHist){
    // Fit and save output data
      double mean = aHist->second->GetMean(), sigma = aHist->second->GetRMS();
      //aHist->second->Fit(lxg_f,"Q");
      o_MipData_BgoDy<<aHist->first<<"\t"<<DmpBgoBase::GetLayerID(aHist->first)<<"\t"<<DmpBgoBase::GetBarID(aHist->first)<<"\t"<<DmpBgoBase::GetSideID(aHist->first)<<"\t";
      for(int ip=0;ip<lxg_f->GetNumberFreeParameters();++ip){
        o_MipData_BgoDy<<"\t"<<lxg_f->GetParameter(ip);
      }
      o_MipData_BgoDy<<std::endl;
      aHist->second->Write();
      delete aHist->second;
  }
  o_MipData_BgoDy.close();

  // create output txtfile      PSD
  histFile->cd("Psd");
  name = "Psd_"+gRootIOSvc->GetInputStem()+".mip";
  o_MipData_PsdDy.open(name.c_str(),std::ios::out);
  o_MipData_PsdDy<<gRootIOSvc->GetInputFileName()<<std::endl;
  o_MipData_PsdDy<<DmpTimeConvertor::Second2Date(fFirstEvtTime)<<std::endl;
  o_MipData_PsdDy<<DmpTimeConvertor::Second2Date(fLastEvtTime)<<std::endl;
  o_MipData_PsdDy<<"globalDynodeID\tlayer\tstrip\tside\t\tWidth\tMP\tArea\tGSigma"<<std::endl;
  lxg_f->SetRange(50,1200);
  for(std::map<short,TH1D*>::iterator aHist=fPsdMipsHist_Dy.begin();aHist!=fPsdMipsHist_Dy.end();++aHist){
    // Fit and save output data
      double mean = aHist->second->GetMean(), sigma = aHist->second->GetRMS();
      //aHist->second->Fit(lxg_f,"Q");
      o_MipData_PsdDy<<aHist->first<<"\t"<<DmpPsdBase::GetLayerID(aHist->first)<<"\t"<<DmpPsdBase::GetStripID(aHist->first)<<"\t"<<DmpPsdBase::GetSideID(aHist->first);
      for(int ip=0;ip<lxg_f->GetNumberFreeParameters();++ip){
        o_MipData_PsdDy<<"\t"<<lxg_f->GetParameter(ip);
      }
      o_MipData_PsdDy<<std::endl;
      aHist->second->Write();
      delete aHist->second;
  }
  o_MipData_PsdDy.close();
  }

  delete histFile;
  return true;
}

void DmpAlgCalibrationMips::_calBarOrStrip()
{
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
      if(fEvtBgo->fADC.at(i0) > 1500 || fEvtBgo->fADC.at(i1) > 1500){
        continue;
      }
      fBgoMipsHist_Bar[DmpBgoBase::ConstructGlobalBarID(l,b)]->Fill(TMath::Sqrt(fEvtBgo->fADC.at(i0) * fEvtBgo->fADC.at(i1)));
      fBgoMipsHist_Dy[gid_dy_s0]->Fill(fEvtBgo->fADC.at(i0));
      fBgoMipsHist_Dy[gid_dy_s1]->Fill(fEvtBgo->fADC.at(i1));
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
      if(fEvtPsd->fADC.at(i0) > 1500 || fEvtPsd->fADC.at(i1) > 1500){
        continue;
      }
      fPsdMipsHist_Bar[DmpPsdBase::ConstructGlobalStripID(l,b)]->Fill(TMath::Sqrt(fEvtPsd->fADC.at(i0) * fEvtPsd->fADC.at(i1)));
      fPsdMipsHist_Dy[gid_dy_s0]->Fill(fEvtPsd->fADC.at(i0));
      fPsdMipsHist_Dy[gid_dy_s1]->Fill(fEvtPsd->fADC.at(i1));
    }
  }
}

void DmpAlgCalibrationMips::_calDynode()
{
  return;
}

void DmpAlgCalibrationMips::SetCalibrationMode(int i)
{
  if(i = 1){
    fEnableCaliDynode = true;
  }else if(i = 2){
    fEnableCaliBar = true;
  }else{ // i <= 0
    fEnableCaliBar = true;
    fEnableCaliDynode = true;
  }
}


