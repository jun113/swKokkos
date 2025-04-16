#ifndef SRC_KOKKOS_DOT_PRODUCT_HPP_
#define SRC_KOKKOS_DOT_PRODUCT_HPP_

#include "Kokkos_Core.hpp"

template<typename T>
class FunctorDotPro1D {
 public:
  using KokView = Kokkos::View<T **>;
	FunctorDotPro1D (const KokView &v_a, const KokView &v_b)
	  : v_a_(v_a), v_b_(v_b) {}
	KOKKOS_INLINE_FUNCTION
	void operator () (const uint64_t &idx, T &result) const {
    result += v_a_.data()[idx] * v_b_.data()[idx];
		return ;
	}
 private:
  KokView v_a_, v_b_;
};

template<typename T>
class FunctorDotPro2D {
 public:
  using KokView = Kokkos::View<T **>;
	FunctorDotPro2D (const KokView &v_a, const KokView &v_b)
	  : v_a_(v_a), v_b_(v_b) {}
	KOKKOS_INLINE_FUNCTION
	void operator () (
			const uint64_t &i0, const uint64_t &i1, T &result) const {
    result += v_a_(i0, i1) * v_b_(i0, i1);
		return ;
	}
 private:
  KokView v_a_, v_b_;
};

// KOKKOS_REGISTER_REDUCE_2D(my_dot_pro_2d_f, FunctorDotPro2D<float>)
KOKKOS_REGISTER_REDUCE_2D(my_dot_pro_2d_d, FunctorDotPro2D<double>)
KOKKOS_REGISTER_REDUCE_1D(my_dot_pro_1d_d, FunctorDotPro1D<double>)

#endif  // SRC_KOKKOS_DOT_PRODUCT_HPP_