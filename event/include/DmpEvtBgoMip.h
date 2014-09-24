/*
 *  $Id: DmpEvtBgoMip.h, 2014-09-11 15:41:06 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 11/09/2014
*/

#ifndef DmpEvtBgoMip_H
#define DmpEvtBgoMip_H

#include "TObject.h"

class DmpEvtBgoMip : public TObject{
/*
 *  DmpEvtBgoMip
 *
 */
public:
  DmpEvtBgoMip();
  DmpEvtBgoMip(const DmpEvtBgoMip &r);
  DmpEvtBgoMip(const DmpEvtBgoMip *&r);
  ~DmpEvtBgoMip();

  void Reset();
  //void SetUsedFileName(const std::string &n){fUsedFile = n;}
  //void SetStartTime(const int &t){fStartTime = t;}
  //void SetStopTime(const int &t){fStopTime = t;}
  //void LoadTimeRange(int &b,int &e)const{b=fStartTime;e=fStopTime;}
  //std::string UsedFile()const{return fUsedFile;}

//private:
public:
  std::string   UsedFileName;      // file name of raw data
  int   StartTime;     // the time of the first event used to cal. ped
  int   StopTime;      // 
  std::vector<short>    GlobalDynodeID;
  std::vector<double>   Mean;
  std::vector<double>   Sigma;
  
  ClassDef(DmpEvtBgoMip,1)

};

#endif

