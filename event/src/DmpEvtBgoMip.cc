/*
 *  $Id: DmpEvtBgoMip.cc, 2014-09-11 15:40:51 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 11/09/2014
*/

#include "DmpEvtBgoMip.h"

ClassImp(DmpEvtBgoMip)

//-------------------------------------------------------------------
DmpEvtBgoMip::DmpEvtBgoMip(){
  Reset();
}

//-------------------------------------------------------------------
DmpEvtBgoMip::DmpEvtBgoMip(const DmpEvtBgoMip &r){
  Reset();
  UsedFileName = r.UsedFileName;
  StartTime = r.StartTime;
  StopTime = r.StopTime;
  short n = GlobalDynodeID.size();
  for(size_t i = 0;i<n;++i){
    GlobalDynodeID.push_back(r.GlobalDynodeID[i]);
    Mean.push_back(r.Mean[i]);
    Sigma.push_back(r.Sigma[i]);
  }
}

//-------------------------------------------------------------------
DmpEvtBgoMip::DmpEvtBgoMip(const DmpEvtBgoMip *&r){
  Reset();
  UsedFileName = r->UsedFileName;
  StartTime = r->StartTime;
  StopTime = r->StopTime;
  short n = GlobalDynodeID.size();
  for(size_t i = 0;i<n;++i){
    GlobalDynodeID.push_back(r->GlobalDynodeID[i]);
    Mean.push_back(r->Mean[i]);
    Sigma.push_back(r->Sigma[i]);
  }
}

//-------------------------------------------------------------------
DmpEvtBgoMip::~DmpEvtBgoMip()
{
}

//-------------------------------------------------------------------
void DmpEvtBgoMip::Reset()
{
  UsedFileName="NO";
  StartTime = 0;
  StopTime  =0xafffffff;
  GlobalDynodeID.clear();
  Mean.clear();
  Sigma.clear();
}


