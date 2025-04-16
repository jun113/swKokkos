#ifndef SRC_AXPY_FUNCTOR_HPP_
#define SRC_AXPY_FUNCTOR_HPP_

#include "Kokkos_Core.hpp"

#include <stdio.h>

template<typename T>
class FunctorAXPY {
 public:
  using KokView1D = Kokkos::View<T *>;
  FunctorAXPY (const T &alpha,
    const KokView1D &v_x, const KokView1D &v_y)
    : alpha_(alpha), v_x_(v_x), v_y_(v_y) {}
  KOKKOS_INLINE_FUNCTION 
  void operator () (const int &i) const {
    v_y_(i) = alpha_ * v_x_(i) + v_y_(i);
    return ;
  }
 private:
  const T alpha_;
  const KokView1D v_x_, v_y_;
};

KOKKOS_REGISTER_FOR_1D(axpy_f, FunctorAXPY<float>)
KOKKOS_REGISTER_FOR_1D(axpy_d, FunctorAXPY<double>)

#endif // SRC_AXPY_FUNCTOR_HPP_