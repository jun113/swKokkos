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

#ifndef KOKKO_ATHREAD_PARALLEL_RANGE_HPP
#define KOKKO_ATHREAD_PARALLEL_RANGE_HPP

#include "Kokkos_Athread_ParamWrap.h"
#include "Kokkos_Athread_KernelLaunch.h"
#include "Kokkos_Athread_FastThreadSpawn.h"

#include <Kokkos_Parallel.hpp>

#include "athread.h"

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::RangePolicy<Traits...>, Kokkos::Athread> {
 private:
  using Policy = Kokkos::RangePolicy<Traits...>;

  const FunctorType m_functor;
  const Policy m_policy;

  template <class TagType>
  std::enable_if_t<std::is_void<TagType>::value> exec() const {
    static AthreadParamWrap param;

    const int len = m_policy.end() - m_policy.begin();
    // MPE + 64 CPEs = 65 cores
    const int tile = (len + 64) / 65;

    athread_get_key<FunctorType>(1);

    param.range[0][0] = m_policy.begin();
    param.range[0][1] = m_policy.end();
    param.tile[0]     = tile;
    param.functor = reinterpret_cast<const void *>(&m_functor);
    param.hash_value = kokkos_athread_get_hash_host (g_athread_functor_key);

    kokkos_athread_parallel_for_launch (&param);

    const int start = m_policy.begin() + 64 * tile;

    _Pragma("ivdep")
    for (int i = start; i < m_policy.end(); ++i) {
      m_functor(i);
    }
    spawn_proxy_join();
    return ;
  }

  template <class TagType>
  std::enable_if_t<!std::is_void<TagType>::value> exec() const {
    const TagType t{};
    const typename Policy::member_type e = m_policy.end();
    for (typename Policy::member_type i = m_policy.begin(); i < e; ++i) {
      m_functor(t, i);
    }
  }

 public:
  inline void execute() const {
    this->template exec<typename Policy::work_tag>();
  }

  inline ParallelFor(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
};

/*--------------------------------------------------------------------------*/

template <class CombinedFunctorReducerType, class... Traits>
class ParallelReduce<CombinedFunctorReducerType, Kokkos::RangePolicy<Traits...>,
                     Kokkos::Athread> {
 private:
  using Policy      = Kokkos::RangePolicy<Traits...>;
  // using WorkTag     = typename Policy::work_tag;
  const Policy m_policy;

  using FunctorType   = typename CombinedFunctorReducerType::functor_type;
  using ReducerType   = typename CombinedFunctorReducerType::reducer_type;


  using pointer_type   = typename ReducerType::pointer_type;
  using value_type     = typename ReducerType::value_type;
  using reference_type = typename ReducerType::reference_type;

  const CombinedFunctorReducerType m_functor_reducer;
  const pointer_type m_result_ptr;

  // template <class TagType>
  // inline std::enable_if_t<std::is_void<TagType>::value> exec(
  //     reference_type update) const {
  //   const typename Policy::member_type e = m_policy.end();
  //   for (typename Policy::member_type i = m_policy.begin(); i < e; ++i) {
  //     m_functor_reducer.get_functor()(i, update);
  //   }
  // }

  // template <class TagType>
  // inline std::enable_if_t<!std::is_void<TagType>::value> exec(
  //     reference_type update) const {
  //   const TagType t{};

  //   const typename Policy::member_type e = m_policy.end();
  //   for (typename Policy::member_type i = m_policy.begin(); i < e; ++i) {
  //     m_functor_reducer.get_functor()(t, i, update);
  //   }
  // }
  inline void exec(reference_type update) const {
    static AthreadParamWrap param;

    const int len = m_policy.end() - m_policy.begin();
    // MPE + 64 CPEs = 65 cores
    const int tile = (len + 64) / 65;

    athread_get_key<FunctorType>(1);

    param.range[0][0] = m_policy.begin();
    param.range[0][1] = m_policy.end();
    param.tile[0]     = tile;
    param.functor = reinterpret_cast<const void *>(&(m_functor_reducer.get_functor()));
    param.hash_value = kokkos_athread_get_hash_host (g_athread_functor_key);

    param.reduce_result = static_cast<double>(update);
    // CPEs
    kokkos_athread_parallel_reduce_launch (&param);

    const int start = m_policy.begin() + 64 * tile;

    _Pragma("ivdep")
    for (int i = start; i < m_policy.end(); ++i) {
      m_functor_reducer.get_functor()(i, update);
    }
    spawn_proxy_join();
    update += static_cast<reference_type>(param.reduce_result);
    // std::cout<<"TODO ParallelReduce RangePolicy"<<std::endl;
    // exit (EXIT_FAILURE);

    // param.range[0][0] = m_policy.begin();
    // param.range[0][1] = m_policy.end();

    // // athread_get_key<FunctorType>(1, param.num_intv16);
    // param.reduce_result = static_cast<double>(update);
    // param.functor = reinterpret_cast<const void *>(&(m_functor_reducer.get_functor()));
    // // athread_parallel_reduce_launch_1D(&param);

    // update = static_cast<reference_type>(param.reduce_result);

    return ;
  }

 public:
  inline void execute() const {
    if (m_result_ptr == nullptr) {
      printf("error in ParallelReduce, execute(), m_result_ptr == nullptr\n");
      exit(0);
    }

    pointer_type ptr = m_result_ptr;
    reference_type update = m_functor_reducer.get_reducer().init(ptr);

    // this->template exec<WorkTag>(update);
    this->exec(update);

    m_functor_reducer.get_reducer().final(ptr);
    return ;
  }

  template <class ViewType>
  ParallelReduce(const CombinedFunctorReducerType& arg_functor_reducer,
                 const Policy& arg_policy,
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

/*--------------------------------------------------------------------------*/

template <class FunctorType, class... Traits>
class ParallelScan<FunctorType, Kokkos::RangePolicy<Traits...>,
                   Kokkos::Athread> {
 private:
  using Policy  = Kokkos::RangePolicy<Traits...>;
  using WorkTag = typename Policy::work_tag;

  using Analysis =
      FunctorAnalysis<FunctorPatternInterface::SCAN, Policy, FunctorType, void>;

  using pointer_type   = typename Analysis::pointer_type;
  using reference_type = typename Analysis::reference_type;

  const CombinedFunctorReducer<FunctorType, typename Analysis::Reducer>
      m_functor_reducer;
  const Policy m_policy;

  template <class TagType>
  inline std::enable_if_t<std::is_void<TagType>::value> exec(
      reference_type update) const {
    const typename Policy::member_type e = m_policy.end();
    for (typename Policy::member_type i = m_policy.begin(); i < e; ++i) {
      m_functor_reducer.get_functor()(i, update, true);
    }
  }

  template <class TagType>
  inline std::enable_if_t<!std::is_void<TagType>::value> exec(
      reference_type update) const {
    const TagType t{};
    const typename Policy::member_type e = m_policy.end();
    for (typename Policy::member_type i = m_policy.begin(); i < e; ++i) {
      m_functor_reducer.get_functor()(t, i, update, true);
    }
  }

 public:
  inline void execute() const {
    const typename Analysis::Reducer& final_reducer =
        m_functor_reducer.get_reducer();
    const size_t pool_reduce_size  = final_reducer.value_size();
    const size_t team_reduce_size  = 0;  // Never shrinks
    const size_t team_shared_size  = 0;  // Never shrinks
    const size_t thread_local_size = 0;  // Never shrinks

    // Need to lock resize_thread_team_data
    auto* internal_instance = m_policy.space().impl_internal_space_instance();
#if defined(__sw_host__)
    std::lock_guard<std::mutex> lock(
        internal_instance->m_thread_team_data_mutex);
#endif
    internal_instance->resize_thread_team_data(
        pool_reduce_size, team_reduce_size, team_shared_size,
        thread_local_size);

    reference_type update = final_reducer.init(pointer_type(
        internal_instance->m_thread_team_data.pool_reduce_local()));

    this->template exec<WorkTag>(update);
  }

  inline ParallelScan(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_functor_reducer(arg_functor, typename Analysis::Reducer{arg_functor}),
        m_policy(arg_policy) {}
};

/*--------------------------------------------------------------------------*/
template <class FunctorType, class ReturnType, class... Traits>
class ParallelScanWithTotal<FunctorType, Kokkos::RangePolicy<Traits...>,
                            ReturnType, Kokkos::Athread> {
 private:
  using Policy  = Kokkos::RangePolicy<Traits...>;
  using WorkTag = typename Policy::work_tag;

  using Analysis = FunctorAnalysis<FunctorPatternInterface::SCAN, Policy,
                                   FunctorType, ReturnType>;

  using value_type     = typename Analysis::value_type;
  using pointer_type   = typename Analysis::pointer_type;
  using reference_type = typename Analysis::reference_type;

  const CombinedFunctorReducer<FunctorType, typename Analysis::Reducer>
      m_functor_reducer;
  const Policy m_policy;
  const pointer_type m_result_ptr;

  template <class TagType>
  inline std::enable_if_t<std::is_void<TagType>::value> exec(
      reference_type update) const {
    const typename Policy::member_type e = m_policy.end();
    for (typename Policy::member_type i = m_policy.begin(); i < e; ++i) {
      m_functor_reducer.get_functor()(i, update, true);
    }
  }

  template <class TagType>
  inline std::enable_if_t<!std::is_void<TagType>::value> exec(
      reference_type update) const {
    const TagType t{};
    const typename Policy::member_type e = m_policy.end();
    for (typename Policy::member_type i = m_policy.begin(); i < e; ++i) {
      m_functor_reducer.get_functor()(t, i, update, true);
    }
  }

 public:
  inline void execute() {
    const size_t pool_reduce_size =
        m_functor_reducer.get_reducer().value_size();
    const size_t team_reduce_size  = 0;  // Never shrinks
    const size_t team_shared_size  = 0;  // Never shrinks
    const size_t thread_local_size = 0;  // Never shrinks

    // Need to lock resize_thread_team_data
    auto* internal_instance = m_policy.space().impl_internal_space_instance();
#if defined(__sw_host__)
    std::lock_guard<std::mutex> lock(
        internal_instance->m_thread_team_data_mutex);
#endif // defined(__sw_host__)
    internal_instance->resize_thread_team_data(
        pool_reduce_size, team_reduce_size, team_shared_size,
        thread_local_size);

    const typename Analysis::Reducer& final_reducer =
        m_functor_reducer.get_reducer();

    reference_type update = final_reducer.init(pointer_type(
        internal_instance->m_thread_team_data.pool_reduce_local()));

    this->template exec<WorkTag>(update);

    *m_result_ptr = update;
  }

  template <class ViewType,
            class Enable = std::enable_if_t<Kokkos::is_view<ViewType>::value>>
  ParallelScanWithTotal(const FunctorType& arg_functor,
                        const Policy& arg_policy,
                        const ViewType& arg_result_view)
      : m_functor_reducer(arg_functor, typename Analysis::Reducer{arg_functor}),
        m_policy(arg_policy),
        m_result_ptr(arg_result_view.data()) {
    static_assert(
        Kokkos::Impl::MemorySpaceAccess<typename ViewType::memory_space,
                                        Kokkos::HostSpace>::accessible,
        "Kokkos::Athread parallel_scan result must be host-accessible!");
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif
