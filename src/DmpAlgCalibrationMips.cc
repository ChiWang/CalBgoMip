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
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"

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

//-------------------------------------------------------------------
DmpAlgCalibrationMips::DmpAlgCalibrationMips()
 :DmpVAlg("DmpAlgCalibrationMips"),
  fEvtBgo(0),
  fEvtPsd(0),
  _fCutBarNo(2),
  _fCutLayerNo(10),
  _fCutEntries(300)
{
  gRootIOSvc->SetOutFileKey("CalMip");
  std::string root_path = (std::string)getenv("DMPSWSYS")+"/share/Calibration";
  gCore->GetJobOption()->SetOption(this->Name()+"/BgoPedestal",root_path+"/Bgo.ped");
  gCore->GetJobOption()->SetOption(this->Name()+"/PsdPedestal",root_path+"/Psd.ped");
  gCore->GetJobOption()->SetOption(this->Name()+"/BgoRelation",root_path+"/Bgo.rel");
  gCore->GetJobOption()->SetOption(this->Name()+"/PsdRelation",root_path+"/Psd.rel");
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
  // Bgo create Hist map
  for(short l=0;l<DmpParameterBgo::kPlaneNo*2;++l){
    for(short b=0;b<DmpParameterBgo::kBarNo;++b){
      for(short s=0;s<DmpParameterBgo::kSideNo;++s){
        snprintf(name,50,"BgoMip-L%02d_B%02d_S%d",l,b,s);
        fBgoMipsHist[l][b][s] = new TH1D(name,name,250,0,2000);
        fBgoMipsHist[l][b][s]->SetLineColor(s+8);
      }
      snprintf(name,50,"BgoMip-L%02d_B%02d",l,b);
      fBgoMipsHist[l][b][2] = new TH1D(name,name,250,0,2000);
      fBgoMipsHist[l][b][2]->SetLineColor(4);
    }
  }
  // Psd
  for(short l=0;l<DmpParameterPsd::kPlaneNo*2;++l){
    for(short b=0;b<DmpParameterPsd::kStripNo;++b){
      for(short s=0;s<DmpParameterPsd::kSideNo;++s){
        snprintf(name,50,"PsdMip-L%02d_B%02d_S%d",l,b,s);
        fPsdMipsHist[l][b][s] = new TH1D(name,name,250,0,2000);
        fPsdMipsHist[l][b][s]->SetLineColor(s+8);
      }
      snprintf(name,50,"PsdMip-L%02d_B%02d",l,b);
      fPsdMipsHist[l][b][2] = new TH1D(name,name,250,0,2000);
      fPsdMipsHist[l][b][2]->SetLineColor(4);
    }
  }


  DmpParameterSteering steering;
  std::string inputFile = gCore->GetJobOption()->GetValue(this->Name()+"/BgoPedestal");
  bool loadSta = DAMPE::LoadParameters(inputFile,fBgoPed,steering);
  if(not loadSta){
          return false;
  }else{
    gCore->GetJobOption()->SetOption(this->Name()+"/BgoPedestal/t0",steering["StartTime"]);
    gCore->GetJobOption()->SetOption(this->Name()+"/BgoPedestal/t1",steering["StopTime"]);
  }
  inputFile = gCore->GetJobOption()->GetValue(this->Name()+"/BgoRelation");
  loadSta = DAMPE::LoadParameters(inputFile,fBgoRel,steering);
  if(not loadSta){
          return false;
  }else{
    gCore->GetJobOption()->SetOption(this->Name()+"/BgoRelation/t0",steering["StartTime"]);
    gCore->GetJobOption()->SetOption(this->Name()+"/BgoRelation/t1",steering["StopTime"]);
  }

  inputFile = gCore->GetJobOption()->GetValue(this->Name()+"/PsdPedestal");
  loadSta = DAMPE::LoadParameters(inputFile,fPsdPed,steering);
  if(not loadSta){
          return false;
  }else{
    gCore->GetJobOption()->SetOption(this->Name()+"/PsdPedestal/t0",steering["StartTime"]);
    gCore->GetJobOption()->SetOption(this->Name()+"/PsdPedestal/t1",steering["StopTime"]);
  }
  inputFile = gCore->GetJobOption()->GetValue(this->Name()+"/PsdRelation");
  loadSta = DAMPE::LoadParameters(inputFile,fPsdRel,steering);
  if(not loadSta){
          return false;
  }else{
    gCore->GetJobOption()->SetOption(this->Name()+"/PsdRelation/t0",steering["StartTime"]);
    gCore->GetJobOption()->SetOption(this->Name()+"/PsdRelation/t1",steering["StopTime"]);
  }
  for(short l=0;l<DmpParameterBgo::kPlaneNo*2;++l){
    for(short b=0;b<DmpParameterBgo::kBarNo;++b){
      for(short s=0;s<2;++s){
        for(short d=0;d<2;++d){  // only for dy2,5
          short gid_s = DmpBgoBase::ConstructGlobalDynodeID(l,b,s,d*3+2);
          short gid_b = DmpBgoBase::ConstructGlobalDynodeID(l,b,s,(d+1)*3+2);
          //std::cout<<"DEBUG: "<<__FILE__<<"("<<__LINE__<<")\tl="<<l<<"\tb ="<<b<<"\ts="<<s<<"\td="<<d<<"\t\tped size ="<<fBgoPed[gid_s].size()<<"\t\t"<<fBgoPed[gid_b].size()<<"\t\t"<<fBgoRel[gid_s].size()<<std::endl;
          fTotSigmaBgo[gid_s] = TMath::Sqrt(TMath::Power(fBgoPed.at(gid_s).at(1),2) + TMath::Power(fBgoPed.at(gid_b).at(1)*fBgoRel.at(gid_s).at(1),2));
          //fTotSigmaBgo[gid_s] = TMath::Sqrt(TMath::Power(fBgoPed[gid_s].at(5),2) + TMath::Power(fBgoPed[gid_b].at(5)*fBgoRel[gid_s].at(5),2) + TMath::Power(fBgoRel[gid_s].at(4)/2,2));
        }
      }
    }
  }
  for(short l=0;l<DmpParameterPsd::kPlaneNo*2;++l){
    for(short b=0;b<DmpParameterPsd::kStripNo;++b){
      for(short s=0;s<2;++s){
        short gid_s = DmpPsdBase::ConstructGlobalDynodeID(l,b,s,5);
        short gid_b = DmpPsdBase::ConstructGlobalDynodeID(l,b,s,8);
          //std::cout<<"DEBUG: "<<__FILE__<<"("<<__LINE__<<")\tgid = "<<gid_s<<"\tl="<<l<<"\tb ="<<b<<"\ts="<<s<<"\t\tped size ="<<fPsdPed[gid_s].size()<<"\t\t"<<fPsdPed[gid_b].size()<<"\t\t"<<fPsdRel[gid_s].size()<<std::endl;
        fTotSigmaPsd[gid_s] = TMath::Sqrt(TMath::Power(fPsdPed.at(gid_s).at(1),2) + TMath::Power(fPsdPed.at(gid_b).at(1)*fPsdRel.at(gid_s).at(1),2));
      }
    }
  }

  return true;
}

//-------------------------------------------------------------------
bool DmpAlgCalibrationMips::ProcessThisEvent()
{
  if(gCore->GetEventHeader()->EnabledPeriodTrigger() && gCore->GetEventHeader()->GeneratedPeriodTrigger()){
    return false;
  }
  //-------------------------------------------------------------------
  // bgo Mip bar
  std::vector<short>   temp;
  std::map<short, std::vector<short> >  barIDs;
  static short dyid[2]={0,0};
  for(short l=0;l<DmpParameterBgo::kPlaneNo*2;++l){
    GetBarIDOfLayer_bgo(temp,l);
    if(temp.size() > _fCutBarNo) return false;
    barIDs.insert(std::make_pair(l,temp));
  }
  if(barIDs.size() < _fCutLayerNo){
    return false;
  }
  // bgo Mip bar
  for(std::map<short, std::vector<short> >::iterator it=barIDs.begin();it!=barIDs.end();++it){
    for(int i=0;i<it->second.size();++i){
      for(int s=0;s<2;++s){
        dyid[s] = DmpBgoBase::ConstructGlobalDynodeID(it->first,it->second.at(i),s,8);
        fBgoMipsHist[it->first][it->second.at(i)][s]->Fill(fEvtBgo->fADC[dyid[s]]);
      }
      fBgoMipsHist[it->first][it->second.at(i)][2]->Fill(TMath::Sqrt(fEvtBgo->fADC[dyid[0]] * fEvtBgo->fADC[dyid[1]]));
    }
  }
  // Psd Mip bar
  for(short l=0;l<DmpParameterPsd::kPlaneNo*2;++l){
    GetBarIDOfLayer_psd(temp,l);
    for(size_t ib = 0;ib<temp.size();++ib){
      for(int s=0;s<2;++s){
        dyid[s] = DmpPsdBase::ConstructGlobalDynodeID(l,temp[ib],s,8);
        fPsdMipsHist[l][temp[ib]][s]->Fill(fEvtPsd->fADC[dyid[s]]);
      }
      fPsdMipsHist[l][temp[ib]][2]->Fill(TMath::Sqrt(fEvtPsd->fADC[dyid[0]] * fEvtPsd->fADC[dyid[1]]));
    }
  }

  return true;
}

//-------------------------------------------------------------------
bool DmpAlgCalibrationMips::Finalize(){
  TF1 *lxg_f = gMyFunctions->GetLanGau();
  double sv[4]={20,500,5000,60};
  double pllo[4]={8,20,200,10};
  double plhi[4]={300,2000,100000,300};
  lxg_f->SetParameters(sv);
  for (int i=0; i<4; ++i) {
    lxg_f->SetParLimits(i, pllo[i], plhi[i]);
  }
  TFile *histFile = gRootIOSvc->GetOutputRootFile();//new TFile(histFileName.c_str(),"RECREATE");
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(0000);
  std::string epsFileName = gRootIOSvc->GetOutputPath()+gRootIOSvc->GetInputStem()+"_Mips.eps";
  TCanvas *c0 = new TCanvas("c0","c0",0,0,600,800);
  std::string epsFileName0 = epsFileName+"[";
  c0->Print(epsFileName0.c_str());

  // create output txtfile      BGO
  histFile->mkdir("Bgo");
  histFile->cd("Bgo");
  std::string name = "Bgo_"+gRootIOSvc->GetInputStem()+".mip";
  o_MipData_Bgo.open(name.c_str(),std::ios::out);
  o_MipData_Bgo<<Mark_S<<"\nFileName="<<gRootIOSvc->GetInputFileName()<<std::endl;
  o_MipData_Bgo<<"StartTime="<<gCore->GetTimeFirstOutput()<<"\nStopTime="<<gCore->GetTimeLastOutput()<<"\nDetector="<<DmpEDetectorID::kBgo<<"\nType=2"<<std::endl;
  o_MipData_Bgo<<Mark_D<<std::endl;
  TCanvas *c = new TCanvas("c1","",600,800);
  c->Divide(9,8);
  int fitStatus=0;
  for(short l=0;l<DmpParameterBgo::kPlaneNo*2;++l){
    for(short b=0;b<DmpParameterBgo::kBarNo;++b){
      for(short s = 0;s<3;++s){
        c->cd(b*3+s+1);
        gPad->SetGrid();
        if(s ==2){
          o_MipData_Bgo<<DmpBgoBase::ConstructGlobalBarID(l,b);
        }else{
          o_MipData_Bgo<<DmpBgoBase::ConstructGlobalDynodeID(l,b,s,8);
        }
        lxg_f->SetRange(fBgoMipsHist[l][b][s]->GetMean() * 0.3,fBgoMipsHist[l][b][s]->GetMean()*2.5);
        fitStatus = fBgoMipsHist[l][b][s]->Fit(lxg_f,"RQB");
        for(int ip=0;ip<lxg_f->GetNumberFreeParameters();++ip){
          o_MipData_Bgo<<"\t\t"<<lxg_f->GetParameter(ip);
        }
        fBgoMipsHist[l][b][s]->Draw();
        double mpv = lxg_f->GetParameter(1);
        o_MipData_Bgo<<"\t\t"<<lxg_f->GetMaximumX(0.8*mpv,1.2*mpv)<<"\t\t"<<lxg_f->GetChisquare() / lxg_f->GetNDF()<<"\t\t"<<fBgoMipsHist[l][b][s]->GetEntries()<<"\t\t"<<fitStatus<<std::endl;
        if(fitStatus == 0 && fBgoMipsHist[l][b][s]->GetEntries() > _fCutEntries){
          fBgoMipsHist[l][b][s]->Write();
        }
        //delete fBgoMipsHist_Bar[gid_bar];
      }
    }
    c->Print(epsFileName.c_str(),"Portrait");//Update();
  }
  o_MipData_Bgo<<Mark_N<<std::endl;
  o_MipData_Bgo.close();

  // create output txtfile      PSD
  histFile->mkdir("Psd");
  histFile->cd("Psd");
  name= "Psd_"+gRootIOSvc->GetInputStem()+".mip";
  o_MipData_Psd.open(name.c_str(),std::ios::out);
  o_MipData_Psd<<Mark_S<<"\nFileName="<<gRootIOSvc->GetInputFileName()<<std::endl;
  o_MipData_Psd<<"StartTime="<<gCore->GetTimeFirstOutput()<<"\nStopTime="<<gCore->GetTimeLastOutput()<<"\nDetector="<<DmpEDetectorID::kPsd<<"\nType=2"<<std::endl;
  o_MipData_Psd<<Mark_D<<std::endl;
  TCanvas *c2 = new TCanvas("c2","",600,800);
  c2->Divide(9,14);
  for(short l=0;l<DmpParameterPsd::kPlaneNo*2;++l){
    for(short b=0;b<DmpParameterPsd::kStripNo;++b){
      for(short s = 0;s<3;++s){
        c2->cd(b*3+s+1);
        gPad->SetGrid();
        if(s == 2){
          o_MipData_Psd<<DmpPsdBase::ConstructGlobalStripID(l,b);
        }else{
          o_MipData_Psd<<DmpPsdBase::ConstructGlobalDynodeID(l,b,s,8);
        }
        lxg_f->SetRange(fPsdMipsHist[l][b][s]->GetMean() * 0.3,fPsdMipsHist[l][b][s]->GetMean()*2.5);
        fitStatus = fPsdMipsHist[l][b][s]->Fit(lxg_f,"RQB");
        for(int ip=0;ip<lxg_f->GetNumberFreeParameters();++ip){
          o_MipData_Psd<<"\t\t"<<lxg_f->GetParameter(ip);
        }
        fPsdMipsHist[l][b][s]->Draw();
        double mpv = lxg_f->GetParameter(1);
        o_MipData_Psd<<"\t\t"<<lxg_f->GetMaximumX(0.8*mpv,1.2*mpv)<<"\t\t"<<lxg_f->GetChisquare() / lxg_f->GetNDF()<<"\t\t"<<fPsdMipsHist[l][b][s]->GetEntries()<<"\t\t"<<fitStatus<<std::endl;
        if(fitStatus == 0 && fPsdMipsHist[l][b][s]->GetEntries() > _fCutEntries){
          fPsdMipsHist[l][b][s]->Write();
        }
        //delete fBgoMipsHist_Bar[gid_bar];
      }
    }
    c2->Print(epsFileName.c_str(),"Portrait");//Update();
  }
  o_MipData_Psd<<Mark_N<<std::endl;
  o_MipData_Psd.close();

  TCanvas *c3 = new TCanvas("c3","c3");
  std::string epsFileName1 = epsFileName+"]";
  c3->Print(epsFileName1.c_str());

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
      if(fEvtBgo->fADC[gid_dy] > 5*fBgoPed.at(gid_dy).at(1)){
        gid_dy5 = DmpBgoBase::ConstructGlobalDynodeID(l,b,side,5);
        if(fEvtBgo->fADC.find(gid_dy5) != fEvtBgo->fADC.end()){
          if(fEvtBgo->fADC.at(gid_dy5) > 5 * fBgoPed.at(gid_dy5).at(1)){
            double deltaADC_small = fEvtBgo->fADC.at(gid_dy)*fBgoRel.at(gid_dy5).at(1) + fBgoRel.at(gid_dy5).at(0)  - fEvtBgo->fADC.at(gid_dy5);
            if(TMath::Abs(deltaADC_small) > 3*fTotSigmaBgo.at(gid_dy5)){
                    // affected by last event
              continue;
            }
          }
        }
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
  static short gid_dy5 =0;
  if(side == 0  || side == 1){
    for(short b=0;b<DmpParameterPsd::kStripNo;++b){
      gid_dy = DmpPsdBase::ConstructGlobalDynodeID(l,b,side,8);
      if(fEvtPsd->fADC.find(gid_dy) == fEvtPsd->fADC.end()){
        continue;
      }
      if(fEvtPsd->fADC.at(gid_dy) > 5* fPsdPed.at(gid_dy).at(1)){
        gid_dy5 = DmpPsdBase::ConstructGlobalDynodeID(l,b,side,5);
        if(fEvtPsd->fADC.find(gid_dy5) != fEvtPsd->fADC.end()){
          if(fEvtPsd->fADC.at(gid_dy5) > 5*fPsdPed.at(gid_dy5).at(1)){
            double deltaADC_small = fEvtPsd->fADC.at(gid_dy)*fPsdRel.at(gid_dy5).at(1) + fPsdRel.at(gid_dy5).at(0)  - fEvtPsd->fADC.at(gid_dy5);
            if(TMath::Abs(deltaADC_small) > 3*fTotSigmaPsd.at(gid_dy5)){
                    // affected by last event
              continue;
            }
          }
        }
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

void DmpAlgCalibrationMips::SetPedestalFile(std::string ID,std::string f)
{
        TString xx = f;
        if(not xx.Contains(".ped")){
                DmpLogWarning<<f<<"("<<ID<<") is not a pedestal file... will use the defaur one"<<DmpLogEndl;
                return;
        }
  TString Id = ID;
  Id.ToUpper();
  if(Id == "BGO"){
    gCore->GetJobOption()->SetOption(this->Name()+"/BgoPedestal",f);
  }else if(Id == "PSD"){
    gCore->GetJobOption()->SetOption(this->Name()+"/PsdPedestal",f);
  }
}

void DmpAlgCalibrationMips::SetRelationFile(std::string ID,std::string f)
{
        TString xx = f;
        if(not xx.Contains(".rel")){
                DmpLogWarning<<f<<"("<<ID<<") is not a relation file... will use the defaur one"<<DmpLogEndl;
                return;
        }
  TString Id = ID;
  Id.ToUpper();
  if(Id == "BGO"){
    gCore->GetJobOption()->SetOption(this->Name()+"/BgoRelation",f);
  }else if(Id == "PSD"){
    gCore->GetJobOption()->SetOption(this->Name()+"/PsdRelation",f);
  }
}


