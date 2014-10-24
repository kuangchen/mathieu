#include <mathieu_wrapper.h>
#include <array>

using namespace mathieu;

MathieuWrapper::MathieuWrapper(double a, double q, unsigned int n_max):
  m(a, q, n_max) 
{
}

python::list MathieuWrapper::eval(double t) const {
  python::list l;
  for (auto &i: m.eval(t)) l.append(i);

  return l;
}


