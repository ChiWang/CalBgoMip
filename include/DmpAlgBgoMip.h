/*
 *  $Id: DmpAlgBgoMip.h, 2014-09-11 15:45:05 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 11/09/2014
*/

#ifndef DmpAlgBgoMip_H
#define DmpAlgBgoMip_H

#include "DmpVAlg.h"

class DmpEvtHeader;
class DmpEvtBgoRaw;
class DmpEvtBgoMip;
class TH1D;

class DmpAlgBgoMip : public DmpVAlg{
/*
 *  DmpAlgBgoMip
 *
 */
public:
  DmpAlgBgoMip();
  ~DmpAlgBgoMip();

  //void Set(const std::string &type,const std::string &value);
  bool Initialize();
  bool ProcessThisEvent();
  bool Finalize();

private:
  DmpEvtHeader          *fEvtHeader;
  DmpEvtBgoRaw          *fBgoRaw;
  DmpEvtBgoMip          *fBgoMip;
  std::map<short,TH1D*>  fMipHist;          // key is global dynode ID

};

#endif

