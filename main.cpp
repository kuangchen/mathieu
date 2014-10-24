#include <mathieu_wrapper.h>
#include <boost/python.hpp>

namespace mathieu {

  BOOST_PYTHON_MODULE(mathieu) 
  {
    python::class_<MathieuWrapper>("Mathieu", 
				   "Mathieu Function class",
				   python::init<double, double, python::optional<unsigned int> >(python::args("a", "q", "n_max")))
      .add_property("a", &MathieuWrapper::get_a)
      .add_property("q", &MathieuWrapper::get_q)
      .add_property("n_max", &MathieuWrapper::get_n_max)
      .add_property("mu", &MathieuWrapper::get_mu)
      .add_property("wronskian", &MathieuWrapper::get_wronskian)
      .def("eval", &MathieuWrapper::eval);
  }
}
