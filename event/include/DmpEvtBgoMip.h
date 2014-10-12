/*
 *  $Id: DmpEvtBgoMip.h, 2014-10-12 14:59:06 DAMPE $
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
  ~DmpEvtBgoMip();
  void Reset();
  DmpEvtBgoMip &operator=(const DmpEvtBgoMip &r);
  void LoadFrom(DmpEvtBgoMip *r);

public:
  std::vector<short>    GlobalPMTID;
  std::vector<double>   Mean;
  std::vector<double>   Sigma;
  
  ClassDef(DmpEvtBgoMip,1)

};

#endif

