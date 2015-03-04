/*
 *  $Id: DmpAlgCalibrationMips.h, 2015-03-03 23:17:03 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 19/07/2014
*/

#ifndef DmpAlgCalibrationMips_H
#define DmpAlgCalibrationMips_H

#include <map>
#include <fstream>
#include "DmpVAlg.h"
#include "DmpEvtPsdRaw.h"

class DmpEvtHeader;
class DmpEvtBgoRaw;
//class DmpEvtPsdRaw;
class TH1D;

class DmpAlgCalibrationMips : public DmpVAlg{
/*
 *  DmpAlgCalibrationMips
 *
 */
public:
  DmpAlgCalibrationMips();
  ~DmpAlgCalibrationMips();

  bool Initialize();
  bool ProcessThisEvent();
  bool Finalize();

private:
  DmpEvtHeader          *fEvtHeader;
  DmpEvtBgoRaw          *fEvtBgo;   // without pedetal, pure signal
  DmpEvtPsdRaw          *fEvtPsd;

  int                   fFirstEvtTime;      // unit second
  int                   fLastEvtTime;       // unit second

  std::map<short,TH1D*>  fBgoMipsHist_Bar;       // key is global bar ID
  std::ofstream         o_MipData_BgoBar;      //

  std::map<short,TH1D*>  fBgoMipsHist_Dy;        // key is global dynode ID
  std::ofstream         o_MipData_BgoDy;      //

  std::map<short,TH1D*>  fPsdMipsHist_Bar;       // key is global bar ID
  std::ofstream         o_MipData_PsdBar;      //

  std::map<short,TH1D*>  fPsdMipsHist_Dy;        // key is global dynode ID
  std::ofstream         o_MipData_PsdDy;      //

};

#endif

