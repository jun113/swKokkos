#ifndef SRC_FUNCTOR_3D_MULT_HPP_
#define SRC_FUNCTOR_3D_MULT_HPP_

template<typename T>
class Functor3DMult {
 public:
  using KokView = Kokkos::View<T ***>;
  Functor3DMult (const KokView &v_A, const KokView &v_B, const KokView &v_C)
    : v_A_(v_A), v_B_(v_B), v_C_(v_C) {}
  KOKKOS_INLINE_FUNCTION void operator () (
      const uint64_t &i0, const uint64_t &i1, const uint64_t &i2) const {
    v_C_(i0, i1, i2) = v_A_(i0, i1, i2) * v_B_(i0, i1, i2);
    return ;
  }
 private:
  const KokView v_A_, v_B_, v_C_;
};

KOKKOS_REGISTER_FOR_3D(my_3d_mutl_f, Functor3DMult<float>)
KOKKOS_REGISTER_FOR_3D(my_3d_mutl_d, Functor3DMult<double>)

#endif // SRC_FUNCTOR_3D_MULT_HPP_