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
#include "DmpLoadParameters.h"

class DmpEvtBgoRaw;
class TH1D;

class DmpAlgCalibrationMips : public DmpVAlg{
/*
 *  DmpAlgCalibrationMips
 *
 */
public:
  DmpAlgCalibrationMips();
  ~DmpAlgCalibrationMips();

  void SetPedestalFile(std::string detectorID,std::string filename);
  void SetRelationFile(std::string detectorID,std::string filename);
  void SetCutBarNumber(int v){_fCutBarNo = v;}
  void SetCutLayerNumber(int v){_fCutLayerNo = v;}
  void SetCutEntries(int v){_fCutEntries = v;}

  bool Initialize();
  bool ProcessThisEvent();
  bool Finalize();

private:
  DmpEvtBgoRaw          *fEvtBgo;   // without pedetal, pure signal
  DmpEvtPsdRaw          *fEvtPsd;

  std::map<short,std::map<short,std::map<short, TH1D*> > >  fBgoMipsHist;       // key is l, b, side (side == 2 is combined)
  std::ofstream         o_MipData_Bgo;      //

  std::map<short,std::map<short,std::map<short, TH1D*> > >  fPsdMipsHist;       // key is l, b, side (side == 2 is combined)
  std::ofstream         o_MipData_Psd;      //

private:
  DmpParameterHolder    fBgoPed; 
  DmpParameterHolder    fPsdPed; 

  DmpParameterHolder    fBgoRel; 
  DmpParameterHolder    fPsdRel; 
  std::map<short,double>  fTotSigmaBgo;
  std::map<short,double>  fTotSigmaPsd;

private:
  int   _fCutBarNo;
  int   _fCutLayerNo;
  int   _fCutEntries;

private:
  void GetBarIDOfLayer_bgo(std::vector<short> &barIDs, short layerID,short side=-1)const;   // side = 0: check side_0;  side = 1: check side_1;  side = -1, check both
  void GetBarIDOfLayer_psd(std::vector<short> &barID,short layerID,short side=-1)const;   // side = 0: check side_0;  side = 1: check side_1;  side = -1, check both

};

#endif

