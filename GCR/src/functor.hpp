#ifndef GCR_SRC_FUNCTOR_HPP_
#define GCR_SRC_FUNCTOR_HPP_

#include "Kokkos_Core.hpp"

using Layout = Kokkos::LayoutRight;
using ViewInt1D = Kokkos::View<int*,         Layout, Kokkos::MemoryUnmanaged>;
using ViewFloat1D = Kokkos::View<float*,     Layout, Kokkos::MemoryUnmanaged>;
using ViewFloat2D = Kokkos::View<float**,    Layout, Kokkos::MemoryUnmanaged>;
using ViewFloat3D = Kokkos::View<float***,   Layout, Kokkos::MemoryUnmanaged>;
using ViewDouble1D = Kokkos::View<double*,   Layout, Kokkos::MemoryUnmanaged>;
using ViewDouble2D = Kokkos::View<double**,  Layout, Kokkos::MemoryUnmanaged>;
using ViewDouble3D = Kokkos::View<double***, Layout, Kokkos::MemoryUnmanaged>;

class Functor1 {
 public:
  Functor1 (const int &its, const int &kts, const int &jts, const int &kte,
            const int &ims, const int &ime, const int &jms,
            const int &NX, const int &NY, const int &NZ,
            const ViewFloat1D &v_p, const ViewDouble1D &v_x0)
      : its_(its), kts_(kts), jts_(jts), kte_(kte), ims_(ims), ime_(ime), 
      jms_(jms), NX_(NX), NY_(NY), NZ_(NZ), v_p_(v_p), v_x0_(v_x0) {}
  KOKKOS_INLINE_FUNCTION
  void operator () (const int &j, const int &k, const int &i) const {
    v_p_(index4b(i, k, j, 0)) = v_x0_(index_x(i, k, j));
   return ; 
  }
 private:
  const int its_, kts_, jts_, kte_;
  const int ims_, ime_, jms_;
  const int NX_, NY_, NZ_;
  const ViewFloat1D v_p_;
  const ViewDouble1D v_x0_;
// #define index4b(i,k,j,m) (i+(1-its) + ((k)+(1-kts))*NX + ((j)+(1-jts))*NX*NY + (m) * NY*NX*NZ)
  KOKKOS_INLINE_FUNCTION 
  int index4b (const int &i, const int &k, const int &j, const int &m) const {
    return i+1-its_ + (k+1-kts_)*NX_ + (j+1-jts_)*NX_*NY_ + m*NY_*NX_*NZ_;
  }
// #define index_x(i,k,j) (i-ims + (k-kts+1)*(ime-ims+1) + (j-jms)*(ime-ims+1)*(kte-kts+3))
  KOKKOS_INLINE_FUNCTION 
  int index_x (const int &i, const int &k, const int &j) const {
  return i-ims_ + (k-kts_+1)*(ime_-ims_+1) + (j-jms_)*(ime_-ims_+1)*(kte_-kts_+3);
  }
};

class Functor2 {
 public:
  Functor2 (const int &its, const int &kts, const int &jts, const int &ite,
            const int &NX, const int &NY,
            const ViewFloat1D &v_r, const ViewFloat1D &v_f0)
      : its_(its), kts_(kts), jts_(jts), ite_(ite),
          NX_(NX), NY_(NY), v_r_(v_r), v_f0_(v_f0) {}
  KOKKOS_INLINE_FUNCTION
  void operator () (const int &j, const int &k, const int &i) const {
    v_r_(index3(i,k,j)) = v_f0_(index3(i,k,j)) - v_r_(index3(i,k,j));
   return ; 
  }
 private:
  const int its_, kts_, jts_, ite_;
  const int NX_, NY_;
  const ViewFloat1D v_r_, v_f0_;
  // #define index3(i,k,j) ( (i)-its+(k-kts+1)*(ite-its+1)+(j-jts)*(ite-its+1)*NY)
  KOKKOS_INLINE_FUNCTION 
  int index3 (const int &i, const int &k, const int &j) const {
    return i-its_ + (k-kts_+1)*(ite_-its_+1) + (j-jts_)*(ite_-its_+1)*NY_;
  }
};

class Functor3 {
 public:
  Functor3 (const int &m, const int &its, const int &kts, const int &jts, const int &ite,
            const int &NX, const int &NY, const int &NZ,
            const ViewFloat1D &v_p, const ViewFloat1D &v_r)
      : m_(m), its_(its), kts_(kts), jts_(jts), ite_(ite), 
          NX_(NX), NY_(NY), NZ_(NZ), v_p_(v_p), v_r_(v_r) {}
  KOKKOS_INLINE_FUNCTION
  void operator () (const int &j, const int &k, const int &i) const {
    v_p_(index4b(i,k,j,m_)) = v_r_(index3(i,k,j));
    return ; 
  }
 private:
  const int m_, its_, kts_, jts_, ite_;
  const int NX_, NY_, NZ_;
  const ViewFloat1D v_p_, v_r_;
// #define index4b(i,k,j,m) (i+(1-its) + ((k)+(1-kts))*NX + ((j)+(1-jts))*NX*NY + (m) * NY*NX*NZ)
  KOKKOS_INLINE_FUNCTION 
  int index4b (const int &i, const int &k, const int &j, const int &m) const {
    return i+1-its_ + (k+1-kts_)*NX_ + (j+1-jts_)*NX_*NY_ + m*NY_*NX_*NZ_;
  }
  // #define index3(i,k,j) ( (i)-its+(k-kts+1)*(ite-its+1)+(j-jts)*(ite-its+1)*NY)
  KOKKOS_INLINE_FUNCTION 
  int index3 (const int &i, const int &k, const int &j) const {
    return i-its_ + (k-kts_+1)*(ite_-its_+1) + (j-jts_)*(ite_-its_+1)*NY_;
  }
};

class Functor4 {
 public:
  Functor4 (const int &its, const int &kts, const int &jts,
            const int &NX, const int &NY, const int &NZ,
            const int &iter_max,
            const ViewFloat1D &v_p)
      : its_(its), kts_(kts), jts_(jts), NX_(NX), NY_(NY), NZ_(NZ), 
          iter_max_(iter_max), v_p_(v_p) {}
  KOKKOS_INLINE_FUNCTION
  void operator () (const int &j, const int &k, const int &i) const {
    v_p_(index4b(i, k, j, 0)) = v_p_(index4b(i, k, j, iter_max_-1));
   return ; 
  }
 private:
  const int its_, kts_, jts_;
  const int NX_, NY_, NZ_, iter_max_;
  const ViewFloat1D v_p_;
// #define index4b(i,k,j,m) (i+(1-its) + ((k)+(1-kts))*NX + ((j)+(1-jts))*NX*NY + (m) * NY*NX*NZ)
  KOKKOS_INLINE_FUNCTION 
  int index4b (const int &i, const int &k, const int &j, const int &m) const {
    return i+1-its_ + (k+1-kts_)*NX_ + (j+1-jts_)*NX_*NY_ + m*NY_*NX_*NZ_;
  }
};

class Functor5 {
 public:
  Functor5 (const int &its, const int &kts, const int &jts, const int &kte,
            const int &ims, const int &ime, const int &jms,
            const int &NX, const int &NY, const int &NZ,
            const int &m, const double &ac,
            const ViewFloat1D &v_p, const ViewDouble1D &v_x0)
      : its_(its), kts_(kts), jts_(jts), kte_(kte), 
          ims_(ims), ime_(ime), jms_(jms), NX_(NX), NY_(NY), NZ_(NZ), 
              ac_(ac), m_(m), v_p_(v_p), v_x0_(v_x0) {}
  KOKKOS_INLINE_FUNCTION
  void operator () (const int &j, const int &k, const int &i) const {
    v_x0_(index_x(i,k,j)) = v_x0_(index_x(i,k,j)) + ac_ * v_p_(index4b(i,k,j,m_));
   return ; 
  }
 private:
  const int its_, kts_, jts_, kte_;
  const int ims_, ime_, jms_;
  const int NX_, NY_, NZ_, m_;
  const double ac_;
  const ViewFloat1D v_p_;
  const ViewDouble1D v_x0_;
  // #define index4b(i,k,j,m) (i+(1-its) + ((k)+(1-kts))*NX + ((j)+(1-jts))*NX*NY + (m) * NY*NX*NZ)
  KOKKOS_INLINE_FUNCTION 
  int index4b (const int &i, const int &k, const int &j, const int &m) const {
    return i+1-its_ + (k+1-kts_)*NX_ + (j+1-jts_)*NX_*NY_ + m*NY_*NX_*NZ_;
  }
  // #define index_x(i,k,j) (i-ims + (k-kts+1)*(ime-ims+1) + (j-jms)*(ime-ims+1)*(kte-kts+3))
  KOKKOS_INLINE_FUNCTION 
  int index_x (const int &i, const int &k, const int &j) const {
  return i-ims_ + (k-kts_+1)*(ime_-ims_+1) + (j-jms_)*(ime_-ims_+1)*(kte_-kts_+3);
  }
};

class Functor6 {
 public:
  Functor6 (const int &its, const int &kts, const int &jts, 
            const int &ite, const int &jte, const int &NY,
            const int &m, const double &ac,
            const ViewFloat1D &v_r, const ViewFloat1D &v_ap)
      : its_(its), kts_(kts), jts_(jts), ite_(ite), jte_(jte), 
          NY_(NY), ac_(ac), m_(m), v_r_(v_r), v_ap_(v_ap) {}
  KOKKOS_INLINE_FUNCTION
  void operator () (const int &j, const int &k, const int &i) const {
    v_r_(index3(i,k,j)) = v_r_(index3(i,k,j)) - ac_ * v_ap_(index4(i,k,j,m_));
   return ; 
  }
 private:
  const int its_, kts_, jts_, ite_, jte_, NY_, m_;
  const double ac_;
  const ViewFloat1D v_r_, v_ap_;
  // #define index3(i,k,j) ( (i)-its+(k-kts+1)*(ite-its+1)+(j-jts)*(ite-its+1)*NY)
  KOKKOS_INLINE_FUNCTION 
  int index3 (const int &i, const int &k, const int &j) const {
    return i-its_ + (k-kts_+1)*(ite_-its_+1) + (j-jts_)*(ite_-its_+1)*NY_;
  }
  // #define index4(i,k,j,m) (i-its + ((k)+(1-kts))*(ite-its+1) + (j-jts)*NY*(ite-its+1) + (m) * NY*(ite-its+1)*(jte-jts+1))
  KOKKOS_INLINE_FUNCTION 
  int index4 (const int &i, const int &k, const int &j, const int &m) 
      const {
    return i-its_ + (k+(1-kts_))*(ite_-its_+1) + (j-jts_)*NY_*(ite_-its_+1) + m*NY_*(ite_-its_+1)*(jte_-jts_+1);
  }
};

class Functor7 {
 public:
  Functor7 (const ViewDouble1D &v_b, 
            const ViewDouble1D &v_c2,
            const ViewDouble1D &v_aps)
      : v_b_(v_b), v_c2_(v_c2), v_aps_(v_aps) {}
  KOKKOS_INLINE_FUNCTION
  void operator () (const int &l) const {
    v_b_(l) = - v_c2_(l) / v_aps_(l);
   return ; 
  }
 private:
  const ViewDouble1D v_b_, v_c2_, v_aps_;
};

// class Functor8 {
//  public:
//   Functor8 (const int &m, const int &its, const int &kts, const int &jts,
//             const int &NX, const int &NY, const int &NZ,
//             const ViewFloat1D &v_p, const ViewDouble1D &v_b)
//       : m_(m), its_(its), kts_(kts), jts_(jts), 
//           NX_(NX), NY_(NY), NZ_(NZ), v_p_(v_p), v_b_(v_b) {}
//   KOKKOS_INLINE_FUNCTION
//   void operator () (const int &j, const int &l, const int &k, const int &i) 
//       const {
//     v_p_(index4b(i,k,j,m_+1)) += v_b_(l) * v_p_(index4b(i,k,j,l));
//    return ; 
//   }
//  private:
//   const int m_, its_, kts_, jts_;
//   const int NX_, NY_, NZ_;
//   const ViewFloat1D v_p_;
//   const ViewDouble1D v_b_;
//   // #define index4b(i,k,j,m) (i+(1-its) + ((k)+(1-kts))*NX + ((j)+(1-jts))*NX*NY + (m) * NY*NX*NZ)
//   KOKKOS_INLINE_FUNCTION 
//   int index4b (const int &i, const int &k, const int &j, const int &m) const {
//     return i+1-its_ + (k+1-kts_)*NX_ + (j+1-jts_)*NX_*NY_ + m*NY_*NX_*NZ_;
//   }
// };

class Functor8 {
 public:
  Functor8 (const int &m, const int &l, const int &its, const int &kts, const int &jts,
            const int &NX, const int &NY, const int &NZ,
            const ViewFloat1D &v_p, const ViewDouble1D &v_b)
      : m_(m), l_(l), its_(its), kts_(kts), jts_(jts),
          NX_(NX), NY_(NY), NZ_(NZ), v_p_(v_p), v_b_(v_b) {}
  KOKKOS_INLINE_FUNCTION
  void operator () (const int &j, const int &k, const int &i) 
      const {
    v_p_(index4b(i,k,j,m_+1)) += v_b_(l_) * v_p_(index4b(i,k,j,l_));
   return ; 
  }
 private:
  const int m_, l_, its_, kts_, jts_;
  const int NX_, NY_, NZ_;
  const ViewFloat1D v_p_;
  const ViewDouble1D v_b_;
  // #define index4b(i,k,j,m) (i+(1-its) + ((k)+(1-kts))*NX + ((j)+(1-jts))*NX*NY + (m) * NY*NX*NZ)
  KOKKOS_INLINE_FUNCTION 
  int index4b (const int &i, const int &k, const int &j, const int &m) const {
    return i+1-its_ + (k+1-kts_)*NX_ + (j+1-jts_)*NX_*NY_ + m*NY_*NX_*NZ_;
  }
};

class Functor9 {
 public:
  Functor9 (const int &its, const int &kts, const int &jts, const int &ite,
            const int &NX, const int &NY,
            const ViewFloat1D &v_r, const ViewFloat1D &v_f0, const ViewFloat1D &v_ap)
      : its_(its), kts_(kts), jts_(jts), ite_(ite),
          NX_(NX), NY_(NY), v_r_(v_r), v_f0_(v_f0), v_ap_(v_ap) {}
  KOKKOS_INLINE_FUNCTION
  void operator () (const int &j, const int &k, const int &i) const {
    v_r_(index3(i,k,j)) = v_f0_(index3(i,k,j)) - v_ap_(index3(i,k,j));
   return ; 
  }
 private:
  const int its_, kts_, jts_, ite_;
  const int NX_, NY_;
  const ViewFloat1D v_r_, v_f0_, v_ap_;
  // #define index3(i,k,j) ( (i)-its+(k-kts+1)*(ite-its+1)+(j-jts)*(ite-its+1)*NY)
  KOKKOS_INLINE_FUNCTION 
  int index3 (const int &i, const int &k, const int &j) const {
    return i-its_ + (k-kts_+1)*(ite_-its_+1) + (j-jts_)*(ite_-its_+1)*NY_;
  }
};

class MatrixProKtsKte {
 public:
  MatrixProKtsKte (const int &its, const int &kts, const int &jts, 
      const int &ite, const int &NX, const int &NY, 
      const int &offset_a, const int &offset_b, const int &offset_c,
      const ViewFloat1D &v_a, const ViewFloat1D &v_b, const ViewFloat1D &v_c) 
    : its_(its), kts_(kts), jts_(jts), ite_(ite), NX_(NX), NY_(NY), 
        offset_a_(offset_a), offset_b_(offset_b), offset_c_(offset_c),
            v_a_(v_a), v_b_(v_b), v_c_(v_c) {}
  KOKKOS_INLINE_FUNCTION
  void operator () (const int &j, const int &k, const int &i) const {
    v_c_(index3(i,k,j)+offset_c_) = + v_a_(index_a(1 ,i,k,j)+offset_a_) * v_b_(index3b(i  ,k  ,j  )+offset_b_)
                                    + v_a_(index_a(2 ,i,k,j)+offset_a_) * v_b_(index3b(i-1,k  ,j  )+offset_b_)
                                    + v_a_(index_a(3 ,i,k,j)+offset_a_) * v_b_(index3b(i+1,k  ,j  )+offset_b_)
                                    + v_a_(index_a(4 ,i,k,j)+offset_a_) * v_b_(index3b(i  ,k  ,j-1)+offset_b_)
                                    + v_a_(index_a(5 ,i,k,j)+offset_a_) * v_b_(index3b(i  ,k  ,j+1)+offset_b_)
                                    + v_a_(index_a(6 ,i,k,j)+offset_a_) * v_b_(index3b(i+1,k  ,j+1)+offset_b_)
                                    + v_a_(index_a(7 ,i,k,j)+offset_a_) * v_b_(index3b(i+1,k  ,j-1)+offset_b_)
                                    + v_a_(index_a(8 ,i,k,j)+offset_a_) * v_b_(index3b(i-1,k  ,j-1)+offset_b_)
                                    + v_a_(index_a(9 ,i,k,j)+offset_a_) * v_b_(index3b(i-1,k  ,j+1)+offset_b_)
                                    + v_a_(index_a(10,i,k,j)+offset_a_) * v_b_(index3b(i  ,k-1,j  )+offset_b_)
                                    + v_a_(index_a(11,i,k,j)+offset_a_) * v_b_(index3b(i-1,k-1,j  )+offset_b_)
                                    + v_a_(index_a(12,i,k,j)+offset_a_) * v_b_(index3b(i+1,k-1,j  )+offset_b_)
                                    + v_a_(index_a(13,i,k,j)+offset_a_) * v_b_(index3b(i  ,k-1,j-1)+offset_b_)
                                    + v_a_(index_a(14,i,k,j)+offset_a_) * v_b_(index3b(i  ,k-1,j+1)+offset_b_)
                                    + v_a_(index_a(15,i,k,j)+offset_a_) * v_b_(index3b(i  ,k+1,j  )+offset_b_)
                                    + v_a_(index_a(16,i,k,j)+offset_a_) * v_b_(index3b(i-1,k+1,j  )+offset_b_)
                                    + v_a_(index_a(17,i,k,j)+offset_a_) * v_b_(index3b(i+1,k+1,j  )+offset_b_)
                                    + v_a_(index_a(18,i,k,j)+offset_a_) * v_b_(index3b(i  ,k+1,j-1)+offset_b_)
                                    + v_a_(index_a(19,i,k,j)+offset_a_) * v_b_(index3b(i  ,k+1,j+1)+offset_b_);
   return ; 
  }
 private:
  const int its_, kts_, jts_, ite_;
  const int NX_, NY_;
  const int offset_a_, offset_b_, offset_c_;
  const ViewFloat1D v_a_, v_b_, v_c_;

  // #define index3(i,k,j) ( (i)-its+(k-kts+1)*(ite-its+1)+(j-jts)*(ite-its+1)*NY)
  KOKKOS_INLINE_FUNCTION 
  int index3 (const int &i, const int &k, const int &j) const {
    return i-its_ + (k-kts_+1)*(ite_-its_+1) + (j-jts_)*(ite_-its_+1)*NY_;
  }
  // #define index_a(m,i,k,j) (m-1+(i-its)*19+ (k-kts+1)*19*(ite-its+1) + (j-jts)*19*(ite-its+1)*NY)
  KOKKOS_INLINE_FUNCTION 
  int index_a (const int &m, const int &i, const int &k, const int &j) const {
    return m-1 + (i-its_)*19+ (k-kts_+1)*19*(ite_-its_+1) + (j-jts_)*19*(ite_-its_+1)*NY_;
  }
  // #define index3b(i,k,j) (i-its+1 + (k-kts+1)*NX + (j-jts+1)*NX*NY)
  KOKKOS_INLINE_FUNCTION 
  int index3b (const int &i, const int &k, const int &j) const {
    return i-its_+1 + (k-kts_+1)*NX_ + (j-jts_+1)*NX_*NY_;
  }
};

class MatrixProKts {
 public:
  MatrixProKts (const int &its, const int &kts, const int &jts, 
      const int &ite, const int &NX, const int &NY, 
      const int &offset_a, const int &offset_b, const int &offset_c,
      const ViewFloat1D &v_a, const ViewFloat1D &v_b, const ViewFloat1D &v_c) 
    : its_(its), kts_(kts), jts_(jts), ite_(ite), NX_(NX), NY_(NY), 
        offset_a_(offset_a), offset_b_(offset_b), offset_c_(offset_c),
            v_a_(v_a), v_b_(v_b), v_c_(v_c) {}
  KOKKOS_INLINE_FUNCTION
  void operator () (const int &j, const int &i) const {
    const int k = kts_ - 1;
    v_c_(index3(i,k,j)+offset_c_) = + v_a_(index_a(1 ,i,k,j)+offset_a_) * v_b_(index3b(i  ,k  ,j  )+offset_b_)
                                    + v_a_(index_a(2 ,i,k,j)+offset_a_) * v_b_(index3b(i-1,k  ,j  )+offset_b_)
                                    + v_a_(index_a(3 ,i,k,j)+offset_a_) * v_b_(index3b(i+1,k  ,j  )+offset_b_)
                                    + v_a_(index_a(4 ,i,k,j)+offset_a_) * v_b_(index3b(i  ,k  ,j-1)+offset_b_)
                                    + v_a_(index_a(5 ,i,k,j)+offset_a_) * v_b_(index3b(i  ,k  ,j+1)+offset_b_)
                                    + v_a_(index_a(6 ,i,k,j)+offset_a_) * v_b_(index3b(i+1,k  ,j+1)+offset_b_)
                                    + v_a_(index_a(7 ,i,k,j)+offset_a_) * v_b_(index3b(i+1,k  ,j-1)+offset_b_)
                                    + v_a_(index_a(8 ,i,k,j)+offset_a_) * v_b_(index3b(i-1,k  ,j-1)+offset_b_)
                                    + v_a_(index_a(9 ,i,k,j)+offset_a_) * v_b_(index3b(i-1,k  ,j+1)+offset_b_)
                                    + v_a_(index_a(15,i,k,j)+offset_a_) * v_b_(index3b(i  ,k+1,j  )+offset_b_)
                                    + v_a_(index_a(16,i,k,j)+offset_a_) * v_b_(index3b(i-1,k+1,j  )+offset_b_)
                                    + v_a_(index_a(17,i,k,j)+offset_a_) * v_b_(index3b(i+1,k+1,j  )+offset_b_)
                                    + v_a_(index_a(18,i,k,j)+offset_a_) * v_b_(index3b(i  ,k+1,j-1)+offset_b_)
                                    + v_a_(index_a(19,i,k,j)+offset_a_) * v_b_(index3b(i  ,k+1,j+1)+offset_b_);
   return ; 
  }
 private:
  const int its_, kts_, jts_, ite_;
  const int NX_, NY_;
  const int offset_a_, offset_b_, offset_c_;
  const ViewFloat1D v_a_, v_b_, v_c_;

  // #define index3(i,k,j) ( (i)-its+(k-kts+1)*(ite-its+1)+(j-jts)*(ite-its+1)*NY)
  KOKKOS_INLINE_FUNCTION 
  int index3 (const int &i, const int &k, const int &j) const {
    return i-its_ + (k-kts_+1)*(ite_-its_+1) + (j-jts_)*(ite_-its_+1)*NY_;
  }
  // #define index_a(m,i,k,j) (m-1+(i-its)*19+ (k-kts+1)*19*(ite-its+1) + (j-jts)*19*(ite-its+1)*NY)
  KOKKOS_INLINE_FUNCTION 
  int index_a (const int &m, const int &i, const int &k, const int &j) const {
    return m-1 + (i-its_)*19+ (k-kts_+1)*19*(ite_-its_+1) + (j-jts_)*19*(ite_-its_+1)*NY_;
  }
  // #define index3b(i,k,j) (i-its+1 + (k-kts+1)*NX + (j-jts+1)*NX*NY)
  KOKKOS_INLINE_FUNCTION 
  int index3b (const int &i, const int &k, const int &j) const {
    return i-its_+1 + (k-kts_+1)*NX_ + (j-jts_+1)*NX_*NY_;
  }
};

class MatrixProKte {
 public:
  MatrixProKte (const int &its, const int &kts, const int &jts, 
      const int &ite, const int &kte, const int &NX, const int &NY, 
      const int &offset_a, const int &offset_b, const int &offset_c,
      const ViewFloat1D &v_a, const ViewFloat1D &v_b, const ViewFloat1D &v_c) 
    : its_(its), kts_(kts), jts_(jts), ite_(ite), kte_(kte), NX_(NX), NY_(NY), 
        offset_a_(offset_a), offset_b_(offset_b), offset_c_(offset_c),
            v_a_(v_a), v_b_(v_b), v_c_(v_c) {}
  KOKKOS_INLINE_FUNCTION
  void operator () (const int &j, const int &i) const {
    const int k = kte_ + 1;
    v_c_(index3(i,k,j)+offset_c_) = + v_a_(index_a(1 ,i,k,j)+offset_a_) * v_b_(index3b(i  ,k  ,j  )+offset_b_)
                                    + v_a_(index_a(2 ,i,k,j)+offset_a_) * v_b_(index3b(i-1,k  ,j  )+offset_b_)
                                    + v_a_(index_a(3 ,i,k,j)+offset_a_) * v_b_(index3b(i+1,k  ,j  )+offset_b_)
                                    + v_a_(index_a(4 ,i,k,j)+offset_a_) * v_b_(index3b(i  ,k  ,j-1)+offset_b_)
                                    + v_a_(index_a(5 ,i,k,j)+offset_a_) * v_b_(index3b(i  ,k  ,j+1)+offset_b_)
                                    + v_a_(index_a(6 ,i,k,j)+offset_a_) * v_b_(index3b(i+1,k  ,j+1)+offset_b_)
                                    + v_a_(index_a(7 ,i,k,j)+offset_a_) * v_b_(index3b(i+1,k  ,j-1)+offset_b_)
                                    + v_a_(index_a(8 ,i,k,j)+offset_a_) * v_b_(index3b(i-1,k  ,j-1)+offset_b_)
                                    + v_a_(index_a(9 ,i,k,j)+offset_a_) * v_b_(index3b(i-1,k  ,j+1)+offset_b_)
                                    + v_a_(index_a(10,i,k,j)+offset_a_) * v_b_(index3b(i  ,k-1,j  )+offset_b_)
                                    + v_a_(index_a(11,i,k,j)+offset_a_) * v_b_(index3b(i-1,k-1,j  )+offset_b_)
                                    + v_a_(index_a(12,i,k,j)+offset_a_) * v_b_(index3b(i+1,k-1,j  )+offset_b_)
                                    + v_a_(index_a(13,i,k,j)+offset_a_) * v_b_(index3b(i  ,k-1,j-1)+offset_b_)
                                    + v_a_(index_a(14,i,k,j)+offset_a_) * v_b_(index3b(i  ,k-1,j+1)+offset_b_);
   return ; 
  }
 private:
  const int its_, kts_, jts_, ite_, kte_;
  const int NX_, NY_;
  const int offset_a_, offset_b_, offset_c_;
  const ViewFloat1D v_a_, v_b_, v_c_;

  // #define index3(i,k,j) ( (i)-its+(k-kts+1)*(ite-its+1)+(j-jts)*(ite-its+1)*NY)
  KOKKOS_INLINE_FUNCTION 
  int index3 (const int &i, const int &k, const int &j) const {
    return i-its_ + (k-kts_+1)*(ite_-its_+1) + (j-jts_)*(ite_-its_+1)*NY_;
  }
  // #define index_a(m,i,k,j) (m-1+(i-its)*19+ (k-kts+1)*19*(ite-its+1) + (j-jts)*19*(ite-its+1)*NY)
  KOKKOS_INLINE_FUNCTION 
  int index_a (const int &m, const int &i, const int &k, const int &j) const {
    return m-1 + (i-its_)*19+ (k-kts_+1)*19*(ite_-its_+1) + (j-jts_)*19*(ite_-its_+1)*NY_;
  }
  // #define index3b(i,k,j) (i-its+1 + (k-kts+1)*NX + (j-jts+1)*NX*NY)
  KOKKOS_INLINE_FUNCTION 
  int index3b (const int &i, const int &k, const int &j) const {
    return i-its_+1 + (k-kts_+1)*NX_ + (j-jts_+1)*NX_*NY_;
  }
};

class MatrixProKtsKteD {
 public:
  MatrixProKtsKteD (const int &its, const int &kts, const int &jts, 
      const int &ite, const int &NX, const int &NY, 
      const int &offset_a, const int &offset_b, const int &offset_c,
      const ViewFloat1D &v_a, const ViewDouble1D &v_b, const ViewFloat1D &v_c) 
    : its_(its), kts_(kts), jts_(jts), ite_(ite), NX_(NX), NY_(NY), 
        offset_a_(offset_a), offset_b_(offset_b), offset_c_(offset_c),
            v_a_(v_a), v_b_(v_b), v_c_(v_c) {}
  KOKKOS_INLINE_FUNCTION
  void operator () (const int &j, const int &k, const int &i) const {
    v_c_(index3(i,k,j)+offset_c_) = + v_a_(index_a(1 ,i,k,j)+offset_a_) * v_b_(index3b(i  ,k  ,j  )+offset_b_)
                                    + v_a_(index_a(2 ,i,k,j)+offset_a_) * v_b_(index3b(i-1,k  ,j  )+offset_b_)
                                    + v_a_(index_a(3 ,i,k,j)+offset_a_) * v_b_(index3b(i+1,k  ,j  )+offset_b_)
                                    + v_a_(index_a(4 ,i,k,j)+offset_a_) * v_b_(index3b(i  ,k  ,j-1)+offset_b_)
                                    + v_a_(index_a(5 ,i,k,j)+offset_a_) * v_b_(index3b(i  ,k  ,j+1)+offset_b_)
                                    + v_a_(index_a(6 ,i,k,j)+offset_a_) * v_b_(index3b(i+1,k  ,j+1)+offset_b_)
                                    + v_a_(index_a(7 ,i,k,j)+offset_a_) * v_b_(index3b(i+1,k  ,j-1)+offset_b_)
                                    + v_a_(index_a(8 ,i,k,j)+offset_a_) * v_b_(index3b(i-1,k  ,j-1)+offset_b_)
                                    + v_a_(index_a(9 ,i,k,j)+offset_a_) * v_b_(index3b(i-1,k  ,j+1)+offset_b_)
                                    + v_a_(index_a(10,i,k,j)+offset_a_) * v_b_(index3b(i  ,k-1,j  )+offset_b_)
                                    + v_a_(index_a(11,i,k,j)+offset_a_) * v_b_(index3b(i-1,k-1,j  )+offset_b_)
                                    + v_a_(index_a(12,i,k,j)+offset_a_) * v_b_(index3b(i+1,k-1,j  )+offset_b_)
                                    + v_a_(index_a(13,i,k,j)+offset_a_) * v_b_(index3b(i  ,k-1,j-1)+offset_b_)
                                    + v_a_(index_a(14,i,k,j)+offset_a_) * v_b_(index3b(i  ,k-1,j+1)+offset_b_)
                                    + v_a_(index_a(15,i,k,j)+offset_a_) * v_b_(index3b(i  ,k+1,j  )+offset_b_)
                                    + v_a_(index_a(16,i,k,j)+offset_a_) * v_b_(index3b(i-1,k+1,j  )+offset_b_)
                                    + v_a_(index_a(17,i,k,j)+offset_a_) * v_b_(index3b(i+1,k+1,j  )+offset_b_)
                                    + v_a_(index_a(18,i,k,j)+offset_a_) * v_b_(index3b(i  ,k+1,j-1)+offset_b_)
                                    + v_a_(index_a(19,i,k,j)+offset_a_) * v_b_(index3b(i  ,k+1,j+1)+offset_b_);
   return ; 
  }
 private:
  const int its_, kts_, jts_, ite_;
  const int NX_, NY_;
  const int offset_a_, offset_b_, offset_c_;
  const ViewFloat1D v_a_;
  const ViewDouble1D v_b_;
  const ViewFloat1D v_c_;

  // #define index3(i,k,j) ( (i)-its+(k-kts+1)*(ite-its+1)+(j-jts)*(ite-its+1)*NY)
  KOKKOS_INLINE_FUNCTION 
  int index3 (const int &i, const int &k, const int &j) const {
    return i-its_ + (k-kts_+1)*(ite_-its_+1) + (j-jts_)*(ite_-its_+1)*NY_;
  }
  // #define index_a(m,i,k,j) (m-1+(i-its)*19+ (k-kts+1)*19*(ite-its+1) + (j-jts)*19*(ite-its+1)*NY)
  KOKKOS_INLINE_FUNCTION 
  int index_a (const int &m, const int &i, const int &k, const int &j) const {
    return m-1 + (i-its_)*19+ (k-kts_+1)*19*(ite_-its_+1) + (j-jts_)*19*(ite_-its_+1)*NY_;
  }
  // #define index3b(i,k,j) (i-its+1 + (k-kts+1)*NX + (j-jts+1)*NX*NY)
  KOKKOS_INLINE_FUNCTION 
  int index3b (const int &i, const int &k, const int &j) const {
    return i-its_+1 + (k-kts_+1)*NX_ + (j-jts_+1)*NX_*NY_;
  }
};

class MatrixProKtsD {
 public:
  MatrixProKtsD (const int &its, const int &kts, const int &jts, 
      const int &ite, const int &NX, const int &NY, 
      const int &offset_a, const int &offset_b, const int &offset_c,
      const ViewFloat1D &v_a, const ViewDouble1D &v_b, const ViewFloat1D &v_c) 
    : its_(its), kts_(kts), jts_(jts), ite_(ite), NX_(NX), NY_(NY), 
        offset_a_(offset_a), offset_b_(offset_b), offset_c_(offset_c),
            v_a_(v_a), v_b_(v_b), v_c_(v_c) {}
  KOKKOS_INLINE_FUNCTION
  void operator () (const int &j, const int &i) const {
    const int k = kts_ - 1;
    v_c_(index3(i,k,j)+offset_c_) = + v_a_(index_a(1 ,i,k,j)+offset_a_) * v_b_(index3b(i  ,k  ,j  )+offset_b_)
                                    + v_a_(index_a(2 ,i,k,j)+offset_a_) * v_b_(index3b(i-1,k  ,j  )+offset_b_)
                                    + v_a_(index_a(3 ,i,k,j)+offset_a_) * v_b_(index3b(i+1,k  ,j  )+offset_b_)
                                    + v_a_(index_a(4 ,i,k,j)+offset_a_) * v_b_(index3b(i  ,k  ,j-1)+offset_b_)
                                    + v_a_(index_a(5 ,i,k,j)+offset_a_) * v_b_(index3b(i  ,k  ,j+1)+offset_b_)
                                    + v_a_(index_a(6 ,i,k,j)+offset_a_) * v_b_(index3b(i+1,k  ,j+1)+offset_b_)
                                    + v_a_(index_a(7 ,i,k,j)+offset_a_) * v_b_(index3b(i+1,k  ,j-1)+offset_b_)
                                    + v_a_(index_a(8 ,i,k,j)+offset_a_) * v_b_(index3b(i-1,k  ,j-1)+offset_b_)
                                    + v_a_(index_a(9 ,i,k,j)+offset_a_) * v_b_(index3b(i-1,k  ,j+1)+offset_b_)
                                    + v_a_(index_a(15,i,k,j)+offset_a_) * v_b_(index3b(i  ,k+1,j  )+offset_b_)
                                    + v_a_(index_a(16,i,k,j)+offset_a_) * v_b_(index3b(i-1,k+1,j  )+offset_b_)
                                    + v_a_(index_a(17,i,k,j)+offset_a_) * v_b_(index3b(i+1,k+1,j  )+offset_b_)
                                    + v_a_(index_a(18,i,k,j)+offset_a_) * v_b_(index3b(i  ,k+1,j-1)+offset_b_)
                                    + v_a_(index_a(19,i,k,j)+offset_a_) * v_b_(index3b(i  ,k+1,j+1)+offset_b_);
   return ; 
  }
 private:
  const int its_, kts_, jts_, ite_;
  const int NX_, NY_;
  const int offset_a_, offset_b_, offset_c_;
  const ViewFloat1D v_a_;
  const ViewDouble1D v_b_;
  const ViewFloat1D v_c_;

  // #define index3(i,k,j) ( (i)-its+(k-kts+1)*(ite-its+1)+(j-jts)*(ite-its+1)*NY)
  KOKKOS_INLINE_FUNCTION 
  int index3 (const int &i, const int &k, const int &j) const {
    return i-its_ + (k-kts_+1)*(ite_-its_+1) + (j-jts_)*(ite_-its_+1)*NY_;
  }
  // #define index_a(m,i,k,j) (m-1+(i-its)*19+ (k-kts+1)*19*(ite-its+1) + (j-jts)*19*(ite-its+1)*NY)
  KOKKOS_INLINE_FUNCTION 
  int index_a (const int &m, const int &i, const int &k, const int &j) const {
    return m-1 + (i-its_)*19+ (k-kts_+1)*19*(ite_-its_+1) + (j-jts_)*19*(ite_-its_+1)*NY_;
  }
  // #define index3b(i,k,j) (i-its+1 + (k-kts+1)*NX + (j-jts+1)*NX*NY)
  KOKKOS_INLINE_FUNCTION 
  int index3b (const int &i, const int &k, const int &j) const {
    return i-its_+1 + (k-kts_+1)*NX_ + (j-jts_+1)*NX_*NY_;
  }
};

class MatrixProKteD {
 public:
  MatrixProKteD (const int &its, const int &kts, const int &jts, 
      const int &ite, const int &kte, const int &NX, const int &NY, 
      const int &offset_a, const int &offset_b, const int &offset_c,
      const ViewFloat1D &v_a, const ViewDouble1D &v_b, const ViewFloat1D &v_c) 
    : its_(its), kts_(kts), jts_(jts), ite_(ite), kte_(kte), NX_(NX), NY_(NY), 
        offset_a_(offset_a), offset_b_(offset_b), offset_c_(offset_c),
            v_a_(v_a), v_b_(v_b), v_c_(v_c) {}
  KOKKOS_INLINE_FUNCTION
  void operator () (const int &j, const int &i) const {
    const int k = kte_ + 1;
    v_c_(index3(i,k,j)+offset_c_) = + v_a_(index_a(1 ,i,k,j)+offset_a_) * v_b_(index3b(i  ,k  ,j  )+offset_b_)
                                    + v_a_(index_a(2 ,i,k,j)+offset_a_) * v_b_(index3b(i-1,k  ,j  )+offset_b_)
                                    + v_a_(index_a(3 ,i,k,j)+offset_a_) * v_b_(index3b(i+1,k  ,j  )+offset_b_)
                                    + v_a_(index_a(4 ,i,k,j)+offset_a_) * v_b_(index3b(i  ,k  ,j-1)+offset_b_)
                                    + v_a_(index_a(5 ,i,k,j)+offset_a_) * v_b_(index3b(i  ,k  ,j+1)+offset_b_)
                                    + v_a_(index_a(6 ,i,k,j)+offset_a_) * v_b_(index3b(i+1,k  ,j+1)+offset_b_)
                                    + v_a_(index_a(7 ,i,k,j)+offset_a_) * v_b_(index3b(i+1,k  ,j-1)+offset_b_)
                                    + v_a_(index_a(8 ,i,k,j)+offset_a_) * v_b_(index3b(i-1,k  ,j-1)+offset_b_)
                                    + v_a_(index_a(9 ,i,k,j)+offset_a_) * v_b_(index3b(i-1,k  ,j+1)+offset_b_)
                                    + v_a_(index_a(10,i,k,j)+offset_a_) * v_b_(index3b(i  ,k-1,j  )+offset_b_)
                                    + v_a_(index_a(11,i,k,j)+offset_a_) * v_b_(index3b(i-1,k-1,j  )+offset_b_)
                                    + v_a_(index_a(12,i,k,j)+offset_a_) * v_b_(index3b(i+1,k-1,j  )+offset_b_)
                                    + v_a_(index_a(13,i,k,j)+offset_a_) * v_b_(index3b(i  ,k-1,j-1)+offset_b_)
                                    + v_a_(index_a(14,i,k,j)+offset_a_) * v_b_(index3b(i  ,k-1,j+1)+offset_b_);
   return ; 
  }
 private:
  const int its_, kts_, jts_, ite_, kte_;
  const int NX_, NY_;
  const int offset_a_, offset_b_, offset_c_;
  const ViewFloat1D v_a_;
  const ViewDouble1D v_b_;
  const ViewFloat1D v_c_;

  // #define index3(i,k,j) ( (i)-its+(k-kts+1)*(ite-its+1)+(j-jts)*(ite-its+1)*NY)
  KOKKOS_INLINE_FUNCTION 
  int index3 (const int &i, const int &k, const int &j) const {
    return i-its_ + (k-kts_+1)*(ite_-its_+1) + (j-jts_)*(ite_-its_+1)*NY_;
  }
  // #define index_a(m,i,k,j) (m-1+(i-its)*19+ (k-kts+1)*19*(ite-its+1) + (j-jts)*19*(ite-its+1)*NY)
  KOKKOS_INLINE_FUNCTION 
  int index_a (const int &m, const int &i, const int &k, const int &j) const {
    return m-1 + (i-its_)*19+ (k-kts_+1)*19*(ite_-its_+1) + (j-jts_)*19*(ite_-its_+1)*NY_;
  }
  // #define index3b(i,k,j) (i-its+1 + (k-kts+1)*NX + (j-jts+1)*NX*NY)
  KOKKOS_INLINE_FUNCTION 
  int index3b (const int &i, const int &k, const int &j) const {
    return i-its_+1 + (k-kts_+1)*NX_ + (j-jts_+1)*NX_*NY_;
  }
};

class Svrasr {
 public:
  Svrasr (const int &m, const int &ni, const int &nk, const int &nj,
      const ViewFloat1D &v_yl, const ViewFloat1D &v_b) 
    : m_(m), ni_(ni), nk_(nk), nj_(nj), v_yl_(v_yl), v_b_(v_b) {}
  KOKKOS_INLINE_FUNCTION
  void operator () (const int &j) const {
    int i, k;
    for (i = 2; i <= ni_; i++) {
      v_yl_(index_yl(i,1,j,m_)) = v_yl_(index_yl(i,1,j,m_)) 
          - v_b_(index_b(i,1,j,2)) * v_yl_(index_yl(i-1,1,j,m_));
    }
    for (k = 2; k <= nk_; k++) {
      for (i = 1; i <= ni_; i++) {
        v_yl_(index_yl(i,k,j,m_)) = v_yl_(index_yl(i,k,j,m_)) 
            - v_b_(index_b(i,k,j,4)) * v_yl_(index_yl(i,k-1,j,m_));
      }
      for (i = 2; i <= ni_; i++) {
        v_yl_(index_yl(i,k,j,m_)) = v_yl_(index_yl(i,k,j,m_)) 
            - v_b_(index_b(i,k,j,2)) * v_yl_(index_yl(i-1,k,j,m_));
      }
    }
    if (j == nj_) {
      v_yl_(index_yl(ni_,nk_,nj_,m_)) = v_yl_(index_yl(ni_,nk_,nj_,m_))
          * v_b_(index_b(ni_,nk_,nj_,1));
    }
    for (i = ni_-1; i >= 1; i -= 1) {
      v_yl_(index_yl(i,nk_,j,m_)) = (v_yl_[index_yl(i,nk_,j,m_)] - v_b_(index_b(i,nk_,j,3)) 
          * v_yl_(index_yl(i+1,nk_,j,m_))) * v_b_(index_b(i,nk_,j,1));
    }
    for (k = nk_-1; k >= 1; k -= 1) {
      for (i = 1; i <= ni_; i++) {
        v_yl_(index_yl(i,k,j,m_)) = v_yl_(index_yl(i,k,j,m_)) - v_b_(index_b(i,k,j,5))
            * v_yl_(index_yl(i,k+1,j,m_));
      }
      v_yl_(index_yl(ni_,k,j,m_)) = v_yl_(index_yl(ni_,k,j,m_)) * v_b_(index_b(ni_,k,j,1));
      for (i = ni_-1; i >= 1; i -= 1) {
        v_yl_(index_yl(i,k,j,m_)) = (v_yl_(index_yl(i,k,j,m_)) - v_b_(index_b(i,k,j,3))
            * v_yl_(index_yl(i+1,k,j,m_))) * v_b_(index_b(i,k,j,1));
      }
    }
    return ;
  }
 private:
  const int m_, ni_, nk_, nj_;
  const ViewFloat1D v_yl_, v_b_;
  // #define index_b(i,k,j,m) (i-1+(k-1)*ni+(j-1)*ni*nk+(m-1)*ni*nk*nj)
  KOKKOS_INLINE_FUNCTION
  int index_b (const int &i, const int &k, const int &j, const int &m) const {
    return i-1 + (k-1)*ni_ + (j-1)*ni_*nk_ + (m-1)*ni_*nk_*nj_;
  }
  // #define index_yl(i,k,j,m) (i-1+(k-1)*ni+(j-1)*ni*nk +m*(ni*nk*nj))
  KOKKOS_INLINE_FUNCTION
  int index_yl (const int &i, const int &k, const int &j, const int &m) const {
    return i-1 + (k-1)*ni_ + (j-1)*ni_*nk_ + m*(ni_*nk_*nj_);
  }
};

class FunctorReduceD {
 public:
  FunctorReduceD (const int &its, const int &kts,
      const int &jts, const int &ite, const int &NY, const ViewFloat1D &v_r) 
    : its_(its), kts_(kts), jts_(jts), ite_(ite), NY_(NY), 
        v_r_(v_r) {}
  KOKKOS_INLINE_FUNCTION
  void operator () (const int &j, const int &k, const int &i,
      double &d) const {
    d = d + v_r_(index3(i,k,j)) * v_r_(index3(i,k,j));
    return ;
  }
 private:
  const int its_, kts_, jts_, ite_, NY_;
  const ViewFloat1D v_r_;
// #define index3(i,k,j) ( (i)-its+(k-kts+1)*(ite-its+1)+(j-jts)*(ite-its+1)*NY)
  KOKKOS_INLINE_FUNCTION
  int index3 (const int &i, const int &k, const int &j) const {
    return i-its_ + (k-kts_+1)*(ite_-its_+1) + (j-jts_)*(ite_-its_+1)*NY_;
  }
};

class FunctorReduceC11 {
 public:
  FunctorReduceC11 (const int &m, const int &its, const int &kts,
      const int &jts, const int &ite, const int &jte, const int &NY, 
      const ViewFloat1D &v_r, const ViewFloat1D &v_ap) 
    : m_(m), its_(its), kts_(kts), jts_(jts), ite_(ite), jte_(jte), NY_(NY), 
        v_r_(v_r), v_ap_(v_ap) {}
  KOKKOS_INLINE_FUNCTION
  void operator () (const int &j, const int &k, const int &i,
      double &c11) const {
    c11 = c11 + v_r_(index3(i,k,j)) * v_ap_(index4(i,k,j,m_));
    return ;
  }
 private:
  const int m_, its_, kts_, jts_, ite_, jte_, NY_;
  const ViewFloat1D v_r_, v_ap_;
// #define index3(i,k,j) ( (i)-its+(k-kts+1)*(ite-its+1)+(j-jts)*(ite-its+1)*NY)
  KOKKOS_INLINE_FUNCTION
  int index3 (const int &i, const int &k, const int &j) const {
    return i-its_ + (k-kts_+1)*(ite_-its_+1) + (j-jts_)*(ite_-its_+1)*NY_;
  }
// #define index4(i,k,j,m) (i-its + ((k)+(1-kts))*(ite-its+1) + (j-jts)*NY*(ite-its+1) + (m) * NY*(ite-its+1)*(jte-jts+1))
  KOKKOS_INLINE_FUNCTION
  int index4 (const int &i, const int &k, const int &j, const int &m) const {
    return i-its_ + (k+(1-kts_))*(ite_-its_+1) + (j-jts_)*NY_*(ite_-its_+1) + m*NY_*(ite_-its_+1)*(jte_-jts_+1);
  }
};

class FunctorReduceC12 {
 public:
  FunctorReduceC12 (const int &m, const int &its, const int &kts,
      const int &jts, const int &ite, const int &jte, const int &NY, 
      const ViewFloat1D &v_ap) 
    : m_(m), its_(its), kts_(kts), jts_(jts), ite_(ite), jte_(jte), NY_(NY), v_ap_(v_ap) {}
  KOKKOS_INLINE_FUNCTION
  void operator () (const int &j, const int &k, const int &i,
      double &c12) const {
    c12 = c12 + v_ap_(index4(i,k,j,m_)) * v_ap_(index4(i,k,j,m_));
    return ;
  }
 private:
  const int m_, its_, kts_, jts_, ite_, jte_, NY_;
  const ViewFloat1D v_ap_;
// #define index4(i,k,j,m) (i-its + ((k)+(1-kts))*(ite-its+1) + (j-jts)*NY*(ite-its+1) + (m) * NY*(ite-its+1)*(jte-jts+1))
  KOKKOS_INLINE_FUNCTION
  int index4 (const int &i, const int &k, const int &j, const int &m) const {
    return i-its_ + (k+(1-kts_))*(ite_-its_+1) + (j-jts_)*NY_*(ite_-its_+1) + m*NY_*(ite_-its_+1)*(jte_-jts_+1);
  }
};

class FunctorReduceCl {
 public:
  FunctorReduceCl (const int &m, const int &its, const int &kts,
      const int &jts, const int &ite, const int &jte, const int &NY, 
      const ViewFloat1D &v_ar, const ViewFloat1D &v_ap) 
    : m_(m), its_(its), kts_(kts), jts_(jts), ite_(ite), jte_(jte), NY_(NY), 
        v_ar_(v_ar), v_ap_(v_ap) {}
  KOKKOS_INLINE_FUNCTION
  void operator () (const int &j, const int &k, const int &i,
      double &cl) const {
    cl = cl + v_ar_(index3(i,k,j)) * v_ap_(index4(i,k,j,m_));
    return ;
  }
 private:
  const int m_, its_, kts_, jts_, ite_, jte_, NY_;
  const ViewFloat1D v_ar_, v_ap_;
  // #define index3(i,k,j) ( (i)-its+(k-kts+1)*(ite-its+1)+(j-jts)*(ite-its+1)*NY)
  KOKKOS_INLINE_FUNCTION
  int index3 (const int &i, const int &k, const int &j) const {
    return i-its_ + (k-kts_+1)*(ite_-its_+1) + (j-jts_)*(ite_-its_+1)*NY_;
  }
  // #define index4(i,k,j,m) (i-its + ((k)+(1-kts))*(ite-its+1) + (j-jts)*NY*(ite-its+1) + (m) * NY*(ite-its+1)*(jte-jts+1))
  KOKKOS_INLINE_FUNCTION
  int index4 (const int &i, const int &k, const int &j, const int &m) const {
    return i-its_ + (k+(1-kts_))*(ite_-its_+1) + (j-jts_)*NY_*(ite_-its_+1) + m*NY_*(ite_-its_+1)*(jte_-jts_+1);
  }
};

class FunctorNorm2 {
 public:
  FunctorNorm2 (const ViewFloat1D &v_r) : v_r_(v_r) {}
  KOKKOS_INLINE_FUNCTION
  void operator () (const int &i, double &c) const {
    c += v_r_(i) * v_r_(i);
    return ;
  }
 private:
  const ViewFloat1D v_r_;
};

KOKKOS_REGISTER_FOR_3D(functor1, Functor1)
KOKKOS_REGISTER_FOR_3D(functor2, Functor2)
KOKKOS_REGISTER_FOR_3D(functor3, Functor3)
KOKKOS_REGISTER_FOR_3D(functor4, Functor4)
KOKKOS_REGISTER_FOR_3D(functor5, Functor5)
KOKKOS_REGISTER_FOR_3D(functor6, Functor6)
KOKKOS_REGISTER_FOR_1D(functor7, Functor7)
KOKKOS_REGISTER_FOR_3D(functor8, Functor8)
KOKKOS_REGISTER_FOR_3D(functor9, Functor9)
KOKKOS_REGISTER_FOR_3D(MatrixProKtsKte, MatrixProKtsKte)
KOKKOS_REGISTER_FOR_2D(MatrixProKts, MatrixProKts)
KOKKOS_REGISTER_FOR_2D(MatrixProKte, MatrixProKte)
KOKKOS_REGISTER_FOR_3D(MatrixProKtsKteD, MatrixProKtsKteD)
KOKKOS_REGISTER_FOR_2D(MatrixProKtsD, MatrixProKtsD)
KOKKOS_REGISTER_FOR_2D(MatrixProKteD, MatrixProKteD)
KOKKOS_REGISTER_FOR_1D(Svrasr, Svrasr)
KOKKOS_REGISTER_REDUCE_3D(FunctorReduceD, FunctorReduceD)
KOKKOS_REGISTER_REDUCE_3D(FunctorReduceC11, FunctorReduceC11)
KOKKOS_REGISTER_REDUCE_3D(FunctorReduceC12, FunctorReduceC12)
KOKKOS_REGISTER_REDUCE_3D(FunctorReduceCl, FunctorReduceCl)
KOKKOS_REGISTER_REDUCE_1D(FunctorNorm2, FunctorNorm2)

#endif // GCR_SRC_FUNCTOR_HPP_