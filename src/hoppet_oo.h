#ifndef __HOPPET_OO__
#define __HOPPET_OO__
#include "hoppet.h"
#include <vector>

extern "C" {
  void * hoppet_cxx__grid_def__new(double dy, double ymax, int order, double eps);
  void   hoppet_cxx__grid_def__delete(void ** griddef);

  int    hoppet_cxx__grid_def__ny(void * griddef);
  void   hoppet_cxx__grid_def__y_values(void * griddef, double * yvals);
  void   hoppet_cxx__grid_def__x_values(void * griddef, double * xvals);
}

namespace hoppet {

class grid_def {
public:
  grid_def(double dy, double ymax, int order=-5, double eps=1e-7)
    : _ptr(hoppet_cxx__grid_def__new(dy, ymax, order, eps)) {}

  ~grid_def() {if (_ptr) hoppet_cxx__grid_def__delete(&_ptr); }

  int ny() const {return hoppet_cxx__grid_def__ny(_ptr); }

  std::vector<double> y_values() const {
    std::vector<double> yvals(ny()+1);
    hoppet_cxx__grid_def__y_values(_ptr, yvals.data());
    return yvals;
  }
  std::vector<double> x_values() const {
    std::vector<double> xvals(ny()+1);
    hoppet_cxx__grid_def__x_values(_ptr, xvals.data());
    return xvals;
  }

  void * ptr() const { return _ptr; }
private:
  void * _ptr = nullptr;
};  

} // end namespace hoppet
#endif // __HOPPET_OO__