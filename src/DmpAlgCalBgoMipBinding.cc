/*
 *  $Id: DmpAlgCalBgoMipBinding.cc, 2014-09-03 11:35:31 DAMPE $
 *  Author(s):
 *    Chi WANG (chiwang@mail.ustc.edu.cn) 03/09/2014
*/

#include <boost/python.hpp>
#include "DmpAlgBgoMip.h"

BOOST_PYTHON_MODULE(libDmpBgoMip){
  using namespace boost::python;

  class_<DmpAlgBgoMip,boost::noncopyable,bases<DmpVAlg> >("DmpAlgBgoMip",init<>());
}

