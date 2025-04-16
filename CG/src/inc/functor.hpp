#ifndef PCG_SRC_INC_FUNCTOR_HPP_
#define PCG_SRC_INC_FUNCTOR_HPP_

#include <Kokkos_Core.hpp>

using Layout = Kokkos::LayoutRight;
using ViewInt1D = Kokkos::View<int*,       Layout, Kokkos::MemoryUnmanaged>;
using ViewFloat1D = Kokkos::View<float*,   Layout, Kokkos::MemoryUnmanaged>;
using ViewDouble1D = Kokkos::View<double*, Layout, Kokkos::MemoryUnmanaged>;

class FunctorSpmv {
 public:
  FunctorSpmv (
      const ViewInt1D &v_row_off,
      const ViewInt1D &v_cols,
      const ViewDouble1D &v_p,
      const ViewDouble1D &v_Ap,
      const ViewDouble1D &v_matrx_data) 
    : v_row_off_(v_row_off), v_cols_(v_cols), 
        v_p_(v_p), v_Ap_(v_Ap), v_matrx_data_(v_matrx_data) {}
  KOKKOS_INLINE_FUNCTION
  void operator () (const int i) const {
		const int row_start = v_row_off_(i);
		const int row_end   = v_row_off_(i + 1);
    v_Ap_(i) = 0.0;
		for (int j = row_start; j < row_end; ++j) {
			v_Ap_(i) += v_matrx_data_(j) * v_p_(v_cols_(j));
		}
    return ;
  }
 private:
  const ViewInt1D v_row_off_, v_cols_;
  const ViewDouble1D v_p_, v_Ap_, v_matrx_data_;
};

template<typename T>
class FunctorAXPY {
  using View1D = Kokkos::View<T *, Layout, Kokkos::MemoryUnmanaged>;
 public:
  FunctorAXPY (const T &alpha,
    const View1D &v_x, const View1D &v_y, const View1D &v_z)
    : alpha_(alpha), v_x_(v_x), v_y_(v_y), v_z_(v_z) {}
  KOKKOS_INLINE_FUNCTION 
  void operator () (const int &i) const {
    v_x_(i) = alpha_ * v_y_(i) + v_z_(i);
    return ;
  }
 private:
  const T alpha_;
  const View1D v_x_, v_y_, v_z_;
};

class MyReduce {
 public:
  MyReduce (const ViewDouble1D &v) : v_a_(v), v_b_(v) {}
  MyReduce (const ViewDouble1D &v_a, const ViewDouble1D &v_b) 
      : v_a_(v_a), v_b_(v_b) {}
  KOKKOS_INLINE_FUNCTION 
  void operator () (const int &i, double &result) const {
    result += v_a_(i) * v_b_(i);
    return ;
  }
 private:
  const ViewDouble1D v_a_, v_b_;
};

KOKKOS_REGISTER_FOR_1D(spmv, FunctorSpmv)
KOKKOS_REGISTER_FOR_1D(axpy_f, FunctorAXPY<float>)
KOKKOS_REGISTER_FOR_1D(axpy_d, FunctorAXPY<double>)
KOKKOS_REGISTER_REDUCE_1D(my_reduce, MyReduce)

#endif // PCG_SRC_INC_FUNCTOR_HPP_