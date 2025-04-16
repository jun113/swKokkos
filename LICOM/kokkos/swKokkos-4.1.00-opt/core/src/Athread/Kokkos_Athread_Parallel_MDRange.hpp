//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKO_ATHREAD_PARALLEL_MDRANGE_HPP
#define KOKKO_ATHREAD_PARALLEL_MDRANGE_HPP

#include "Kokkos_Athread_ParamWrap.h"
#include "Kokkos_Athread_KernelLaunch.h"
#include "Kokkos_Athread_FastThreadSpawn.h"

#include <athread.h>
#include <Kokkos_Parallel.hpp>
#include <KokkosExp_MDRangePolicy.hpp>

#include <iostream>
#include <type_traits>

namespace Kokkos {
namespace Impl {

// template <class FunctorType, class... Traits>
// class ParallelFor<FunctorType, Kokkos::MDRangePolicy<Traits...>,
//                   Kokkos::Athread> {
//  private:
//   using MDRangePolicy = Kokkos::MDRangePolicy<Traits...>;
//   using Policy        = typename MDRangePolicy::impl_range_policy;

//   using iterate_type = typename Kokkos::Impl::HostIterateTile<
//       MDRangePolicy, FunctorType, typename MDRangePolicy::work_tag, void>;

//   const iterate_type m_iter;

//   void exec() const {
//     const typename Policy::member_type e = m_iter.m_rp.m_num_tiles;
//     for (typename Policy::member_type i = 0; i < e; ++i) {
//       m_iter(i);
//     }
//   }

//  public:
//   inline void execute() const { this->exec(); }
//   template <typename Policy, typename Functor>
//   static int max_tile_size_product(const Policy&, const Functor&) {
//     /**
//      * 1024 here is just our guess for a reasonable max tile size,
//      * it isn't a hardware constraint. If people see a use for larger
//      * tile size products, we're happy to change this.
//      */
//     return 1024;
//   }
//   inline ParallelFor(const FunctorType& arg_functor,
//                      const MDRangePolicy& arg_policy)
//       : m_iter(arg_policy, arg_functor) {}
// };

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::MDRangePolicy<Traits...>,
                  Kokkos::Athread> {
 private:
  using MDPolicy = Kokkos::MDRangePolicy<Traits...>;

  using iterate_type = typename Kokkos::Impl::HostIterateTile<
      MDPolicy, FunctorType, typename MDPolicy::work_tag, void>;

  using array_index_type = typename MDPolicy::array_index_type;
  using index_type       = typename MDPolicy::index_type;
  using LaunchBounds     = typename MDPolicy::launch_bounds;

  const FunctorType m_functor;
  const MDPolicy    m_policy;

  void exec() const {
    static AthreadParamWrap param;
    // 2D
    if constexpr (std::is_invocable<FunctorType, 
        index_type, index_type>::value) {
      athread_get_key<FunctorType>(2);
   
      param.range[0][0] = m_policy.m_lower[0];
      param.range[0][1] = m_policy.m_upper[0];
      param.range[1][0] = m_policy.m_lower[1];
      param.range[1][1] = m_policy.m_upper[1];
      param.tile[0]     = m_policy.m_tile[0];
      param.tile[1]     = m_policy.m_tile[1];

      param.total_tiles[1] =  (param.range[1][1] - param.range[1][0] + param.tile[1] - 1
          ) / param.tile[1];
      param.total_tiles[0] = ((param.range[0][1] - param.range[0][0] + param.tile[0] - 1
          ) / param.tile[0]) * param.total_tiles[1];
   
      param.hash_value = kokkos_athread_get_hash_host (g_athread_functor_key);
      param.functor = reinterpret_cast<const void *>(&m_functor);
   
      // CPEs
      kokkos_athread_parallel_for_launch (&param);
   
      int index_tile0, index_tile1;
      int start0, start1, end0, end1;
      // MPE
      for (int index_tile = 64; index_tile < param.total_tiles[0];
          index_tile += 65) {
        index_tile0 = index_tile / param.total_tiles[1];
        index_tile1 = index_tile % param.total_tiles[1];
        start0 = param.range[0][0] + index_tile0 * param.tile[0];
        end0   = ((start0 + param.tile[0]) < param.range[0][1]) ?
            (start0 + param.tile[0]) : param.range[0][1];
        start1 = param.range[1][0] + index_tile1 * param.tile[1];
        end1   = ((start1 + param.tile[1]) < param.range[1][1]) ?
            (start1 + param.tile[1]) : param.range[1][1];
        _Pragma("ivdep")
        for (int i0 = start0; i0 < end0; ++ i0) {
          _Pragma("ivdep")
          for (int i1 = start1; i1 < end1; ++ i1) {
            m_functor (i0, i1);
          }
        }
      }
    // 3D
    } else if constexpr (std::is_invocable<FunctorType, 
        index_type, index_type, index_type>::value) {
      athread_get_key<FunctorType>(3);

      param.range[0][0] = m_policy.m_lower[0];
      param.range[0][1] = m_policy.m_upper[0];
      param.range[1][0] = m_policy.m_lower[1];
      param.range[1][1] = m_policy.m_upper[1];
      param.range[2][0] = m_policy.m_lower[2];
      param.range[2][1] = m_policy.m_upper[2];
      param.tile[0]     = m_policy.m_tile[0];
      param.tile[1]     = m_policy.m_tile[1];
      param.tile[2]     = m_policy.m_tile[2];

      // printf ("lower: %d %d %d\n", param.range[0][0],param.range[1][0],param.range[2][0]);
      // printf ("upper: %d %d %d\n", param.range[0][1],param.range[1][1],param.range[2][1]);
      // printf ("tile: %d %d %d\n", m_policy.m_tile[0],m_policy.m_tile[1],m_policy.m_tile[2]);

      param.total_tiles[2] =  (param.range[2][1] - param.range[2][0] + param.tile[2] - 1) 
          / param.tile[2];
      param.total_tiles[1] = ((param.range[1][1] - param.range[1][0] + param.tile[1] - 1) 
          / param.tile[1]) * param.total_tiles[2];
      param.total_tiles[0] = ((param.range[0][1] - param.range[0][0] + param.tile[0] - 1) 
          / param.tile[0]) * param.total_tiles[1];

      param.hash_value = kokkos_athread_get_hash_host (g_athread_functor_key);
      param.functor = reinterpret_cast<const void *>(&m_functor);

      // CPEs
      kokkos_athread_parallel_for_launch (&param);

      int index_tile0, index_tile1, index_tile2;
      int tile12;
      int start0, start1, start2;
      int end0, end1, end2;
      // MPE
      for (int index_tile = 64; index_tile < param.total_tiles[0]; index_tile += 65) {
        index_tile0 = index_tile / param.total_tiles[1];
        tile12      = index_tile % param.total_tiles[1];
        index_tile1 = tile12     / param.total_tiles[2];
        index_tile2 = tile12     % param.total_tiles[2];
        start0 = param.range[0][0] + index_tile0 * param.tile[0];
        end0   = ((start0 + param.tile[0]) < param.range[0][1]) 
            ? (start0 + param.tile[0]) : param.range[0][1];
        start1 = param.range[1][0] + index_tile1 * param.tile[1];
        end1   = ((start1 + param.tile[1]) < param.range[1][1]) 
            ? (start1 + param.tile[1]) : param.range[1][1];
        start2 = param.range[2][0] + index_tile2 * param.tile[2];
        end2   = ((start2 + param.tile[2]) < param.range[2][1]) 
            ? (start2 + param.tile[2]) : param.range[2][1];
        _Pragma("ivdep")
        for (int i0 = start0; i0 < end0; ++ i0) {
          _Pragma("ivdep")
          for (int i1 = start1; i1 < end1; ++ i1) {
            _Pragma("ivdep")
            for (int i2 = start2; i2 < end2; ++ i2) {
              m_functor (i0, i1, i2);
            }
          }
        }
      }
    // 4D
    } else if constexpr (std::is_invocable<FunctorType, 
        index_type, index_type, index_type, index_type>::value) {
      athread_get_key<FunctorType>(4);

      //   _Pragma("ivdep")
      //   for (int i0 = m_policy.m_lower[0]; i0 < m_policy.m_upper[0]; ++ i0) {
      //     _Pragma("ivdep")
      //     for (int i1 = m_policy.m_lower[1]; i1 < m_policy.m_upper[1]; ++ i1) {
      //       _Pragma("ivdep")
      //       for (int i2 = m_policy.m_lower[2]; i2 < m_policy.m_upper[2]; ++ i2) {
      //         _Pragma("ivdep")
      //         for (int i3 = m_policy.m_lower[3]; i3 < m_policy.m_upper[3]; ++ i3) {
      //           m_functor (i0, i1, i2, i3);
      //         }
      //       }
      //     }
      //   }
      // return ;

      param.range[0][0] = m_policy.m_lower[0];
      param.range[0][1] = m_policy.m_upper[0];
      param.range[1][0] = m_policy.m_lower[1];
      param.range[1][1] = m_policy.m_upper[1];
      param.range[2][0] = m_policy.m_lower[2];
      param.range[2][1] = m_policy.m_upper[2];
      param.range[3][0] = m_policy.m_lower[3];
      param.range[3][1] = m_policy.m_upper[3];
      param.tile[0]     = m_policy.m_tile[0];
      param.tile[1]     = m_policy.m_tile[1];
      param.tile[2]     = m_policy.m_tile[2];
      param.tile[3]     = m_policy.m_tile[3];

      param.total_tiles[3] =  (param.range[3][1] - param.range[3][0] + param.tile[3] - 1) 
          / param.tile[3];
      param.total_tiles[2] = ((param.range[2][1] - param.range[2][0] + param.tile[2] - 1) 
          / param.tile[2]) * param.total_tiles[3];
      param.total_tiles[1] = ((param.range[1][1] - param.range[1][0] + param.tile[1] - 1) 
          / param.tile[1]) * param.total_tiles[2];
      param.total_tiles[0] = ((param.range[0][1] - param.range[0][0] + param.tile[0] - 1) 
          / param.tile[0]) * param.total_tiles[1];

      param.hash_value = kokkos_athread_get_hash_host (g_athread_functor_key);
      param.functor = reinterpret_cast<const void *>(&m_functor);

      // CPEs
      kokkos_athread_parallel_for_launch (&param);

      int index_tile0, index_tile1, index_tile2, index_tile3;
      int tile123, tile23;
      int start0, start1, start2, start3;
      int end0, end1, end2, end3;
      // MPE
      for (int index_tile = 64; index_tile < param.total_tiles[0]; 
          index_tile += 65) {
        index_tile0 = index_tile / param.total_tiles[1];
        tile123     = index_tile % param.total_tiles[1];
        index_tile1 = tile123    / param.total_tiles[2];
        tile23      = index_tile % param.total_tiles[2];
        index_tile2 = tile23     / param.total_tiles[3];
        index_tile3 = tile23     % param.total_tiles[3];
        start0 = param.range[0][0] + index_tile0 * param.tile[0];
        end0   = ((start0 + param.tile[0]) < param.range[0][1]) 
            ? (start0 + param.tile[0]) : param.range[0][1];
        start1 = param.range[1][0] + index_tile1 * param.tile[1];
        end1   = ((start1 + param.tile[1]) < param.range[1][1]) 
            ? (start1 + param.tile[1]) : param.range[1][1];
        start2 = param.range[2][0] + index_tile2 * param.tile[2];
        end2   = ((start2 + param.tile[2]) < param.range[2][1]) 
            ? (start2 + param.tile[2]) : param.range[2][1];
        start3 = param.range[3][0] + index_tile3 * param.tile[3];
        end3   = ((start3 + param.tile[3]) < param.range[3][1]) 
            ? (start3 + param.tile[3]) : param.range[3][1];
        _Pragma("ivdep")
        for (int i0 = start0; i0 < end0; ++ i0) {
          _Pragma("ivdep")
          for (int i1 = start1; i1 < end1; ++ i1) {
            _Pragma("ivdep")
            for (int i2 = start2; i2 < end2; ++ i2) {
              _Pragma("ivdep")
              for (int i3 = start3; i3 < end3; ++ i3) {
                m_functor (i0, i1, i2, i3);
              }
            }
          }
        }
      }
    // 5D
    } else if constexpr (std::is_invocable<FunctorType, 
        index_type, index_type, index_type, index_type, index_type>::value) {
      athread_get_key<FunctorType>(5);
      std::cout<<"TODO ParallelFor MDRangePolicy::rank == 5"<<std::endl;
      exit (EXIT_FAILURE);
    } else {
      std::cout<<"error in ParallelReduce->exec: "<<__PRETTY_FUNCTION__<<std::endl;
      exit (EXIT_FAILURE);
    }
    spawn_proxy_join();
    // athread_join ();
    return ;
  }

 public:
  inline void execute() const { this->exec(); }
  template <typename Policy, typename Functor>
  static int max_tile_size_product(const MDPolicy&, const Functor&) {
    /**
     * 1024 here is just our guess for a reasonable max tile size,
     * it isn't a hardware constraint. If people see a use for larger
     * tile size products, we're happy to change this.
     */
    return 1024;
  }
  inline ParallelFor (const FunctorType& arg_functor,
                      const MDPolicy& arg_policy)
      : m_functor (arg_functor), m_policy (arg_policy) {}
};

template <class CombinedFunctorReducerType, class... Traits>
class ParallelReduce<CombinedFunctorReducerType,
                     Kokkos::MDRangePolicy<Traits...>, Kokkos::Athread> {
 private:
  using MDRangePolicy = Kokkos::MDRangePolicy<Traits...>;
  using Policy        = typename MDRangePolicy::impl_range_policy;
  using index_type    = typename MDRangePolicy::index_type;
  using FunctorType   = typename CombinedFunctorReducerType::functor_type;
  using ReducerType   = typename CombinedFunctorReducerType::reducer_type;

  using WorkTag = typename MDRangePolicy::work_tag;

  using pointer_type   = typename ReducerType::pointer_type;
  using value_type     = typename ReducerType::value_type;
  using reference_type = typename ReducerType::reference_type;

  const CombinedFunctorReducerType m_functor_reducer;
  const MDRangePolicy m_policy;
  const pointer_type m_result_ptr;

  inline void exec(reference_type update) const {
    static AthreadParamWrap param;
    // 2D
    if constexpr (std::is_invocable<FunctorType, 
        index_type, index_type, reference_type>::value) {
      athread_get_key<FunctorType>(2);
   
      param.range[0][0] = m_policy.m_lower[0];
      param.range[0][1] = m_policy.m_upper[0];
      param.range[1][0] = m_policy.m_lower[1];
      param.range[1][1] = m_policy.m_upper[1];
      param.tile[0]     = m_policy.m_tile[0];
      param.tile[1]     = m_policy.m_tile[1];

      param.total_tiles[1] =  (param.range[1][1] - param.range[1][0] + param.tile[1] - 1
          ) / param.tile[1];
      param.total_tiles[0] = ((param.range[0][1] - param.range[0][0] + param.tile[0] - 1
          ) / param.tile[0]) * param.total_tiles[1];
   
      param.hash_value = kokkos_athread_get_hash_host (g_athread_functor_key);
      param.functor = reinterpret_cast<const void *>(&(m_functor_reducer.get_functor()));
   
      // TODO
      //param.reduce_result = update;
      param.reduce_double = update;
      // CPEs
      kokkos_athread_parallel_reduce_launch (&param);
   
      int index_tile0, index_tile1;
      int start0, start1, end0, end1;
      // MPE
      for (int index_tile = 64; index_tile < param.total_tiles[0];
          index_tile += 65) {
        index_tile0 = index_tile / param.total_tiles[1];
        index_tile1 = index_tile % param.total_tiles[1];
        start0 = param.range[0][0] + index_tile0 * param.tile[0];
        end0   = ((start0 + param.tile[0]) < param.range[0][1]) ?
            (start0 + param.tile[0]) : param.range[0][1];
        start1 = param.range[1][0] + index_tile1 * param.tile[1];
        end1   = ((start1 + param.tile[1]) < param.range[1][1]) ?
            (start1 + param.tile[1]) : param.range[1][1];
        _Pragma("ivdep")
        for (int i0 = start0; i0 < end0; ++ i0) {
          _Pragma("ivdep")
          for (int i1 = start1; i1 < end1; ++ i1) {
            m_functor_reducer.get_functor()(i0, i1, update);
          }
        }
      }
      spawn_proxy_join();
      // athread_join ();
      update += static_cast<reference_type>(param.reduce_double);
      return ;
    // 3D
    } else if constexpr (std::is_invocable<FunctorType, 
        index_type, index_type, index_type, reference_type>::value) {
      athread_get_key<FunctorType>(3);

      param.range[0][0] = m_policy.m_lower[0];
      param.range[0][1] = m_policy.m_upper[0];
      param.range[1][0] = m_policy.m_lower[1];
      param.range[1][1] = m_policy.m_upper[1];
      param.range[2][0] = m_policy.m_lower[2];
      param.range[2][1] = m_policy.m_upper[2];
      param.tile[0]     = m_policy.m_tile[0];
      param.tile[1]     = m_policy.m_tile[1];
      param.tile[2]     = m_policy.m_tile[2];

      param.total_tiles[2] =  (param.range[2][1] - param.range[2][0] + param.tile[2] - 1) 
          / param.tile[2];
      param.total_tiles[1] = ((param.range[1][1] - param.range[1][0] + param.tile[1] - 1) 
          / param.tile[1]) * param.total_tiles[2];
      param.total_tiles[0] = ((param.range[0][1] - param.range[0][0] + param.tile[0] - 1) 
          / param.tile[0]) * param.total_tiles[1];

      param.hash_value = kokkos_athread_get_hash_host (g_athread_functor_key);
      param.functor = reinterpret_cast<const void *>(&(m_functor_reducer.get_functor()));

      param.reduce_double = update;

      // CPEs
      kokkos_athread_parallel_reduce_launch (&param);

      int index_tile0, index_tile1, index_tile2;
      int tile12;
      int start0, start1, start2;
      int end0, end1, end2;
      // MPE
      for (int index_tile = 64; index_tile < param.total_tiles[0]; index_tile += 65) {
        index_tile0 = index_tile / param.total_tiles[1];
        tile12      = index_tile % param.total_tiles[1];
        index_tile1 = tile12     / param.total_tiles[2];
        index_tile2 = tile12     % param.total_tiles[2];
        start0 = param.range[0][0] + index_tile0 * param.tile[0];
        end0   = ((start0 + param.tile[0]) < param.range[0][1]) 
            ? (start0 + param.tile[0]) : param.range[0][1];
        start1 = param.range[1][0] + index_tile1 * param.tile[1];
        end1   = ((start1 + param.tile[1]) < param.range[1][1]) 
            ? (start1 + param.tile[1]) : param.range[1][1];
        start2 = param.range[2][0] + index_tile2 * param.tile[2];
        end2   = ((start2 + param.tile[2]) < param.range[2][1]) 
            ? (start2 + param.tile[2]) : param.range[2][1];
        _Pragma("ivdep")
        for (int i0 = start0; i0 < end0; ++ i0) {
          _Pragma("ivdep")
          for (int i1 = start1; i1 < end1; ++ i1) {
            _Pragma("ivdep")
            for (int i2 = start2; i2 < end2; ++ i2) {
              m_functor_reducer.get_functor() (i0, i1, i2, update);
            }
          }
        }
      }
      spawn_proxy_join();
      // athread_join ();
      update += static_cast<reference_type>(param.reduce_double);
      // update = static_cast<reference_type>(param.reduce_double);
      return ;
    // 4D
    } else if constexpr (std::is_invocable<FunctorType, 
        index_type, index_type, index_type, index_type, reference_type>::value) {
      athread_get_key<FunctorType>(4);
      std::cout<<"TODO ParallelReduce MDRangePolicy::rank == 4"<<std::endl;
      exit (EXIT_FAILURE);
    // 5D
    } else if constexpr (std::is_invocable<FunctorType, 
        index_type, index_type, index_type, index_type, index_type, reference_type>::value) {
      athread_get_key<FunctorType>(5);
      std::cout<<"TODO ParallelReduce MDRangePolicy::rank == 5"<<std::endl;
      exit (EXIT_FAILURE);
    } else {
      std::cout<<"error in ParallelReduce->exec: "<<__PRETTY_FUNCTION__<<std::endl;
      exit (EXIT_FAILURE);
    }
  }

 public:
  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy&, const Functor&) {
    /**
     * 1024 here is just our guess for a reasonable max tile size,
     * it isn't a hardware constraint. If people see a use for larger
     * tile size products, we're happy to change this.
     */
    return 1024;
  }
  inline void execute() const {
    const ReducerType& reducer = m_functor_reducer.get_reducer();
    if (m_result_ptr == nullptr) {
      printf("error in ParallelReduce, execute(), m_result_ptr == nullptr\n");
      exit(0);
    }
    pointer_type ptr = m_result_ptr;
    reference_type update = reducer.init(ptr);
    this->exec(update);
    reducer.final(ptr);
  }

  template <class ViewType>
  ParallelReduce(const CombinedFunctorReducerType& arg_functor_reducer,
                 const MDRangePolicy& arg_policy,
                 const ViewType& arg_result_view)
      : m_functor_reducer(arg_functor_reducer),
        m_policy(arg_policy),
        m_result_ptr(arg_result_view.data()) {
    static_assert(Kokkos::is_view<ViewType>::value,
                  "Kokkos::Athread reduce result must be a View");

    static_assert(
        Kokkos::Impl::MemorySpaceAccess<typename ViewType::memory_space,
                                        Kokkos::HostSpace>::accessible,
        "Kokkos::Athread reduce result must be a View accessible from "
        "HostSpace");
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif
