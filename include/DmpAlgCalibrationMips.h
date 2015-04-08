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
  DmpEvtBgoRaw          *fEvtBgo;   // without pedetal, pure signal
  DmpEvtPsdRaw          *fEvtPsd;

  std::map<short,TH1D*>  fBgoMipsHist_Bar;       // key is global bar ID
  std::map<short,TH1D*>  fBgoMipsHist_Dy;        // key is global dynode ID
  std::ofstream         o_MipData_Bgo;      //

  std::map<short,TH1D*>  fPsdMipsHist_Bar;       // key is global bar ID
  std::map<short,TH1D*>  fPsdMipsHist_Dy;        // key is global dynode ID
  std::ofstream         o_MipData_Psd;      //

private:
  std::vector<short>  GetBarIDOfLayer_bgo(short layerID,short side=-1)const;   // side = 0: check side_0;  side = 1: check side_1;  side = -1, check both
  std::vector<short>  GetBarIDOfLayer_psd(short layerID,short side=-1)const;   // side = 0: check side_0;  side = 1: check side_1;  side = -1, check both

};

#endif

