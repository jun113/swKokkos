#ifndef SRC_FUNCTOR_STENCIL_HPP_
#define SRC_FUNCTOR_STENCIL_HPP_

#include "Kokkos_Core.hpp"

template<typename T>
class Functor5PointStencil {
 public:
  using KokView2D = Kokkos::View<T **>;
  Functor5PointStencil (const KokView2D &v_A, const KokView2D &v_B)
    : v_A_(v_A), v_B_(v_B) {}
  KOKKOS_INLINE_FUNCTION void operator () (
      const uint64_t &i0, const uint64_t &i1) const {
	  v_B_(i0, i1) += v_A_((i0  ), (i1  ))
                  + v_A_((i0-1), (i1  ))
                  + v_A_((i0+1), (i1  ))
                  + v_A_((i0  ), (i1-1))
                  + v_A_((i0  ), (i1+1));
    return ;
  }
 private:
  const KokView2D v_A_, v_B_;
};

KOKKOS_REGISTER_FOR_2D(my_5_point_stencil_f, Functor5PointStencil<float>)
KOKKOS_REGISTER_FOR_2D(my_5_point_stencil_d, Functor5PointStencil<double>)

#endif // SRC_FUNCTOR_STENCIL_HPP_