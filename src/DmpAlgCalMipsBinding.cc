/*
 *  $Id: DmpAlgCalMipsBinding.cc, 2015-03-03 18:36:32 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 03/09/2014
*/

#include <boost/python.hpp>
#include "DmpAlgCalibrationMips.h"

BOOST_PYTHON_MODULE(libDmpCalMip){
  using namespace boost::python;

  class_<DmpAlgCalibrationMips,boost::noncopyable,bases<DmpVAlg> >("DmpAlgCalibrationMips",init<>())
    .def("SetHistParameters",    &DmpAlgCalibrationMips::SetHistParameters)
    ;
}

