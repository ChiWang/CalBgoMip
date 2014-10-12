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
DmpEvtBgoMip& DmpEvtBgoMip::operator=(const DmpEvtBgoMip &r){
  Reset();
  GlobalPMTID = r.GlobalPMTID;
  Mean = r.Mean;
  Sigma = r.Sigma;
}

//-------------------------------------------------------------------
DmpEvtBgoMip::~DmpEvtBgoMip(){
  Reset();
}

//-------------------------------------------------------------------
void DmpEvtBgoMip::Reset(){
  GlobalPMTID.clear();
  Mean.clear();
  Sigma.clear();
}

//-------------------------------------------------------------------
void DmpEvtBgoMip::LoadFrom(DmpEvtBgoMip *r){
  Reset();
  GlobalPMTID = r->GlobalPMTID;
  Mean = r->Mean;
  Sigma = r->Sigma;
}

