#ifndef __mathieu_h_
#define __mathieu_h_

#include <iostream>
#include <boost/python/tuple.hpp>
#include <array>
#include <string>
#include <map>

namespace mathieu {

  namespace python = boost::python;

  class Mathieu {
    double a, q;
    unsigned int n_max;
    double mu;
    std::vector<double> coeff_p, coeff_n;

    double wronskian;
    
  public:
    Mathieu(double a, double q, unsigned int n_max);
    
    inline double get_a() const { return a; }
    inline double get_q() const { return q; }
    inline unsigned int get_n_max() const { return n_max; }
    inline double get_wronskian() const { return wronskian; }
    inline double get_mu() const { return mu; }

    std::array<double, 4> eval(double x) const;

  private:
    void calculate_coeff();
    void calculate_misc();
  };

}

#endif
