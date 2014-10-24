#ifndef _mathieu_wrapper_h
#define _mathieu_wrapper_h

#include <boost/python/numeric.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python.hpp>
#include <mathieu.h>
#include <array>

namespace mathieu {
  namespace python = boost::python;

  class Mathieu;

  class MathieuWrapper {
  private:
    mathieu::Mathieu m;
    
  public:
    MathieuWrapper(double a, double q, unsigned int n_max=30);
    
    double get_a() const { return m.get_a(); }
    double get_q() const { return m.get_q(); }
    unsigned int get_n_max() const { return m.get_n_max(); }
    double  get_wronskian() const { return m.get_wronskian(); }
    double get_mu() const { return m.get_mu(); }

    
python::list eval(double t) const;

  };

}

#endif
