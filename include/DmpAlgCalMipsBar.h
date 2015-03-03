/*
 *  $Id: DmpAlgCalMipsBar.h, 2015-03-03 09:06:15 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 19/07/2014
*/

#ifndef DmpAlgCalMipsBar_H
#define DmpAlgCalMipsBar_H

#include <map>
#include <fstream>
#include "DmpVAlg.h"
#include "DmpEvtPsdRaw.h"

class DmpEvtHeader;
class DmpEvtBgoRaw;
//class DmpEvtPsdRaw;
class TH1D;

class DmpAlgCalMipsBar : public DmpVAlg{
/*
 *  DmpAlgCalMipsBar
 *
 */
public:
  DmpAlgCalMipsBar();
  ~DmpAlgCalMipsBar();

  //void Set(const std::string &type,const std::string &value);
  bool Initialize();
  bool ProcessThisEvent();
  bool Finalize();

private:
  DmpEvtHeader          *fEvtHeader;
  DmpEvtBgoRaw          *fEvtBgo;   // without pedetal, pure signal
  DmpEvtPsdRaw          *fEvtPsd;

  int                   fFirstEvtTime;      // unit second
  int                   fLastEvtTime;       // unit second

  std::map<short,TH1D*>  fBgoMipHist;       // key is global bar ID
  std::ofstream         OutBgoMipData;      //
  std::map<short,TH1D*>  fPsdMipHist;       // key is global bar ID
  std::ofstream         OutPsdMipData;      //

};

#endif

