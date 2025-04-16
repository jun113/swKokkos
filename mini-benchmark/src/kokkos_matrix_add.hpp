#ifndef SRC_FUNCTOR_MATRIX_ADD_HPP_
#define SRC_FUNCTOR_MATRIX_ADD_HPP_

#include "Kokkos_Core.hpp"

template<typename T>
class FunctorMatrixAdd {
 public:
  using KokView2D = Kokkos::View<T **>;
  FunctorMatrixAdd (const KokView2D &v_A, 
      const KokView2D &v_B, const KokView2D &v_C)
    : v_A_(v_A), v_B_(v_B), v_C_(v_C) {}
  KOKKOS_INLINE_FUNCTION void operator () (
      const uint64_t &i0, const uint64_t &i1) const {
    v_C_(i0, i1) = v_A_(i0, i1) + v_B_(i0, i1);
    return ;
  }
 private:
  const KokView2D v_A_, v_B_, v_C_;
};

KOKKOS_REGISTER_FOR_2D(my_matx_add_f, FunctorMatrixAdd<float>)
KOKKOS_REGISTER_FOR_2D(my_matx_add_d, FunctorMatrixAdd<double>)

#endif // SRC_FUNCTOR_MATRIX_ADD_HPP_