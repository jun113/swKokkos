// __INCLUDE_FUNCTOR_HPP_START__
#include "/home/export/online1/mdt00/shisuan/swiap/wjl/swKokkos-paper/swKokkos-AD/LICOM/src/kokkos-impl/kokkos_readyc.hpp"
#include "/home/export/online1/mdt00/shisuan/swiap/wjl/swKokkos-paper/swKokkos-AD/LICOM/src/kokkos-impl/kokkos_bclinc.hpp"
#include "/home/export/online1/mdt00/shisuan/swiap/wjl/swKokkos-paper/swKokkos-AD/LICOM/src/kokkos-impl/kokkos_icesnow.hpp"
#include "/home/export/online1/mdt00/shisuan/swiap/wjl/swKokkos-paper/swKokkos-AD/LICOM/src/kokkos-impl/kokkos_convadj.hpp"
#include "/home/export/online1/mdt00/shisuan/swiap/wjl/swKokkos-paper/swKokkos-AD/LICOM/src/kokkos-impl/kokkos_barotr.hpp"
#include "/home/export/online1/mdt00/shisuan/swiap/wjl/swKokkos-paper/swKokkos-AD/LICOM/src/kokkos-impl/kokkos_readyt.hpp"
#include "/home/export/online1/mdt00/shisuan/swiap/wjl/swKokkos-paper/swKokkos-AD/LICOM/src/kokkos-impl/kokkos_jra_daily.hpp"
#include "/home/export/online1/mdt00/shisuan/swiap/wjl/swKokkos-paper/swKokkos-AD/LICOM/src/kokkos-impl/kokkos_nextstep.hpp"
#include "/home/export/online1/mdt00/shisuan/swiap/wjl/swKokkos-paper/swKokkos-AD/LICOM/src/kokkos-impl/kokkos_tracer.hpp"
#include "/home/export/online1/mdt00/shisuan/swiap/wjl/swKokkos-paper/swKokkos-AD/LICOM/src/util/pop_haloupdate.hpp"
// __INCLUDE_FUNCTOR_HPP_END__

#include "Athread/Kokkos_Athread_Utils_Slave.h"
#include "Athread/Kokkos_Athread_HashTable.hpp"

#include "Athread/Kokkos_Athread_RegisterFunction.hpp"
#include "Athread/Kokkos_Athread_ParamWrap.h"
#include "Athread/Kokkos_Athread_FastThreadSpawn.h"

#include "simd.h"
#include "slave.h"

#include <string>

#include <time.h>

//ATHREAD DMA transfer definition 
__thread_local unsigned int D_COUNT  = 0;
__thread_local crts_rply_t  dma_rply = 0;
__thread_local crts_rply_t  l_rply   = 0;
__thread_local crts_rply_t  r_rply   = 0;

__thread_local_fix AthreadParamWrap kokkos_athread_local_param __attribute__ ((aligned(64)));

// Buffer for the operator of a reduction between slave cores
__thread_local double buf_reduce[64] __attribute__ ((aligned(64)));
__thread_local int buf_reduce_int __attribute__ ((aligned(64)));


extern "C" void register_kernel() {

  kokkos_athread_hash_table_init ();
  // __REGISTER_START__
  AthreadHashNode* func0 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func0->key = arr_char_to_arr_int (
      "FunctorReadyc12");
  func0->fp = FunctorReadyc1_2D;
  kokkos_athread_hash_table_insert (func0);
  AthreadHashNode* func1 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func1->key = arr_char_to_arr_int (
      "FunctorReadyc23");
  func1->fp = FunctorReadyc2_3D;
  kokkos_athread_hash_table_insert (func1);
  AthreadHashNode* func2 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func2->key = arr_char_to_arr_int (
      "FunctorReadyc33");
  func2->fp = FunctorReadyc3_3D;
  kokkos_athread_hash_table_insert (func2);
  AthreadHashNode* func3 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func3->key = arr_char_to_arr_int (
      "FunctorReadyc43");
  func3->fp = FunctorReadyc4_3D;
  kokkos_athread_hash_table_insert (func3);
  AthreadHashNode* func4 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func4->key = arr_char_to_arr_int (
      "FunctorReadyc513");
  func4->fp = FunctorReadyc51_3D;
  kokkos_athread_hash_table_insert (func4);
  AthreadHashNode* func5 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func5->key = arr_char_to_arr_int (
      "FunctorReadyc522");
  func5->fp = FunctorReadyc52_2D;
  kokkos_athread_hash_table_insert (func5);
  AthreadHashNode* func6 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func6->key = arr_char_to_arr_int (
      "FunctorReadyc533");
  func6->fp = FunctorReadyc53_3D;
  kokkos_athread_hash_table_insert (func6);
  AthreadHashNode* func7 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func7->key = arr_char_to_arr_int (
      "FunctorReadyc542");
  func7->fp = FunctorReadyc54_2D;
  kokkos_athread_hash_table_insert (func7);
  AthreadHashNode* func8 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func8->key = arr_char_to_arr_int (
      "FunctorReadyc553");
  func8->fp = FunctorReadyc55_3D;
  kokkos_athread_hash_table_insert (func8);
  AthreadHashNode* func9 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func9->key = arr_char_to_arr_int (
      "FunctorReadyc63");
  func9->fp = FunctorReadyc6_3D;
  kokkos_athread_hash_table_insert (func9);
  AthreadHashNode* func10 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func10->key = arr_char_to_arr_int (
      "FunctorReadyc72");
  func10->fp = FunctorReadyc7_2D;
  kokkos_athread_hash_table_insert (func10);
  AthreadHashNode* func11 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func11->key = arr_char_to_arr_int (
      "FunctorReadyc83");
  func11->fp = FunctorReadyc8_3D;
  kokkos_athread_hash_table_insert (func11);
  AthreadHashNode* func12 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func12->key = arr_char_to_arr_int (
      "FunctorReadyc93");
  func12->fp = FunctorReadyc9_3D;
  kokkos_athread_hash_table_insert (func12);
  AthreadHashNode* func13 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func13->key = arr_char_to_arr_int (
      "FunctorReadyc113");
  func13->fp = FunctorReadyc11_3D;
  kokkos_athread_hash_table_insert (func13);
  AthreadHashNode* func14 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func14->key = arr_char_to_arr_int (
      "FuncAdvMomCen13");
  func14->fp = FuncAdvMomCen1_3D;
  kokkos_athread_hash_table_insert (func14);
  AthreadHashNode* func15 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func15->key = arr_char_to_arr_int (
      "FuncAdvMomFlu13");
  func15->fp = FuncAdvMomFlu1_3D;
  kokkos_athread_hash_table_insert (func15);
  AthreadHashNode* func16 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func16->key = arr_char_to_arr_int (
      "FuncAdvMomCen23");
  func16->fp = FuncAdvMomCen2_3D;
  kokkos_athread_hash_table_insert (func16);
  AthreadHashNode* func17 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func17->key = arr_char_to_arr_int (
      "FuncAdvMomFlu23");
  func17->fp = FuncAdvMomFlu2_3D;
  kokkos_athread_hash_table_insert (func17);
  AthreadHashNode* func18 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func18->key = arr_char_to_arr_int (
      "FunctorReadyc143");
  func18->fp = FunctorReadyc14_3D;
  kokkos_athread_hash_table_insert (func18);
  AthreadHashNode* func19 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func19->key = arr_char_to_arr_int (
      "FunctorReadyc153");
  func19->fp = FunctorReadyc15_3D;
  kokkos_athread_hash_table_insert (func19);
  AthreadHashNode* func20 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func20->key = arr_char_to_arr_int (
      "FunctorReadyc163");
  func20->fp = FunctorReadyc16_3D;
  kokkos_athread_hash_table_insert (func20);
  AthreadHashNode* func21 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func21->key = arr_char_to_arr_int (
      "FunctorReadyc173");
  func21->fp = FunctorReadyc17_3D;
  kokkos_athread_hash_table_insert (func21);
  AthreadHashNode* func22 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func22->key = arr_char_to_arr_int (
      "FunctorReadyc202");
  func22->fp = FunctorReadyc20_2D;
  kokkos_athread_hash_table_insert (func22);
  AthreadHashNode* func23 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func23->key = arr_char_to_arr_int (
      "FunctorReadyc212");
  func23->fp = FunctorReadyc21_2D;
  kokkos_athread_hash_table_insert (func23);
  AthreadHashNode* func24 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func24->key = arr_char_to_arr_int (
      "FunctorReadyc222");
  func24->fp = FunctorReadyc22_2D;
  kokkos_athread_hash_table_insert (func24);
  AthreadHashNode* func25 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func25->key = arr_char_to_arr_int (
      "FunctorReadyc233");
  func25->fp = FunctorReadyc23_3D;
  kokkos_athread_hash_table_insert (func25);
  AthreadHashNode* func26 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func26->key = arr_char_to_arr_int (
      "FunctorReadyc243");
  func26->fp = FunctorReadyc24_3D;
  kokkos_athread_hash_table_insert (func26);
  AthreadHashNode* func27 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func27->key = arr_char_to_arr_int (
      "FunctorReadyc253");
  func27->fp = FunctorReadyc25_3D;
  kokkos_athread_hash_table_insert (func27);
  AthreadHashNode* func28 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func28->key = arr_char_to_arr_int (
      "FunctorBclinc12");
  func28->fp = FunctorBclinc1_2D;
  kokkos_athread_hash_table_insert (func28);
  AthreadHashNode* func29 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func29->key = arr_char_to_arr_int (
      "FunctorBclinc22");
  func29->fp = FunctorBclinc2_2D;
  kokkos_athread_hash_table_insert (func29);
  AthreadHashNode* func30 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func30->key = arr_char_to_arr_int (
      "FunctorBclinc33");
  func30->fp = FunctorBclinc3_3D;
  kokkos_athread_hash_table_insert (func30);
  AthreadHashNode* func31 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func31->key = arr_char_to_arr_int (
      "FunctorBclinc53");
  func31->fp = FunctorBclinc5_3D;
  kokkos_athread_hash_table_insert (func31);
  AthreadHashNode* func32 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func32->key = arr_char_to_arr_int (
      "FunctorBclinc63");
  func32->fp = FunctorBclinc6_3D;
  kokkos_athread_hash_table_insert (func32);
  AthreadHashNode* func33 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func33->key = arr_char_to_arr_int (
      "FunctorBclinc73");
  func33->fp = FunctorBclinc7_3D;
  kokkos_athread_hash_table_insert (func33);
  AthreadHashNode* func34 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func34->key = arr_char_to_arr_int (
      "FunctorBclinc83");
  func34->fp = FunctorBclinc8_3D;
  kokkos_athread_hash_table_insert (func34);
  AthreadHashNode* func35 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func35->key = arr_char_to_arr_int (
      "FunctorBclinc93");
  func35->fp = FunctorBclinc9_3D;
  kokkos_athread_hash_table_insert (func35);
  AthreadHashNode* func36 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func36->key = arr_char_to_arr_int (
      "FunctorBclinc122");
  func36->fp = FunctorBclinc12_2D;
  kokkos_athread_hash_table_insert (func36);
  AthreadHashNode* func37 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func37->key = arr_char_to_arr_int (
      "FunctorBclinc133");
  func37->fp = FunctorBclinc13_3D;
  kokkos_athread_hash_table_insert (func37);
  AthreadHashNode* func38 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func38->key = arr_char_to_arr_int (
      "FunctorBclinc163");
  func38->fp = FunctorBclinc16_3D;
  kokkos_athread_hash_table_insert (func38);
  AthreadHashNode* func39 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func39->key = arr_char_to_arr_int (
      "FunctorBclinc193");
  func39->fp = FunctorBclinc19_3D;
  kokkos_athread_hash_table_insert (func39);
  AthreadHashNode* func40 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func40->key = arr_char_to_arr_int (
      "FunctorIcesnow12");
  func40->fp = FunctorIcesnow1_2D;
  kokkos_athread_hash_table_insert (func40);
  AthreadHashNode* func41 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func41->key = arr_char_to_arr_int (
      "FunctorIcesnow23");
  func41->fp = FunctorIcesnow2_3D;
  kokkos_athread_hash_table_insert (func41);
  AthreadHashNode* func42 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func42->key = arr_char_to_arr_int (
      "FunctorConvadj12");
  func42->fp = FunctorConvadj1_2D;
  kokkos_athread_hash_table_insert (func42);
  AthreadHashNode* func43 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func43->key = arr_char_to_arr_int (
      "FunctorConvadj24");
  func43->fp = FunctorConvadj2_4D;
  kokkos_athread_hash_table_insert (func43);
  AthreadHashNode* func44 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func44->key = arr_char_to_arr_int (
      "FunctorConvadj32");
  func44->fp = FunctorConvadj3_2D;
  kokkos_athread_hash_table_insert (func44);
  AthreadHashNode* func45 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func45->key = arr_char_to_arr_int (
      "FunctorBarotr13");
  func45->fp = FunctorBarotr1_3D;
  kokkos_athread_hash_table_insert (func45);
  AthreadHashNode* func46 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func46->key = arr_char_to_arr_int (
      "FunctorBarotr22");
  func46->fp = FunctorBarotr2_2D;
  kokkos_athread_hash_table_insert (func46);
  AthreadHashNode* func47 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func47->key = arr_char_to_arr_int (
      "FunctorBarotr32");
  func47->fp = FunctorBarotr3_2D;
  kokkos_athread_hash_table_insert (func47);
  AthreadHashNode* func48 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func48->key = arr_char_to_arr_int (
      "FunctorBarotr42");
  func48->fp = FunctorBarotr4_2D;
  kokkos_athread_hash_table_insert (func48);
  AthreadHashNode* func49 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func49->key = arr_char_to_arr_int (
      "FunctorBarotr52");
  func49->fp = FunctorBarotr5_2D;
  kokkos_athread_hash_table_insert (func49);
  AthreadHashNode* func50 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func50->key = arr_char_to_arr_int (
      "FunctorBarotr62");
  func50->fp = FunctorBarotr6_2D;
  kokkos_athread_hash_table_insert (func50);
  AthreadHashNode* func51 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func51->key = arr_char_to_arr_int (
      "FunctorBarotr72");
  func51->fp = FunctorBarotr7_2D;
  kokkos_athread_hash_table_insert (func51);
  AthreadHashNode* func52 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func52->key = arr_char_to_arr_int (
      "FunctorBarotr82");
  func52->fp = FunctorBarotr8_2D;
  kokkos_athread_hash_table_insert (func52);
  AthreadHashNode* func53 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func53->key = arr_char_to_arr_int (
      "FunctorBarotr92");
  func53->fp = FunctorBarotr9_2D;
  kokkos_athread_hash_table_insert (func53);
  AthreadHashNode* func54 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func54->key = arr_char_to_arr_int (
      "FunctorBarotr112");
  func54->fp = FunctorBarotr11_2D;
  kokkos_athread_hash_table_insert (func54);
  AthreadHashNode* func55 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func55->key = arr_char_to_arr_int (
      "FunctorBarotr122");
  func55->fp = FunctorBarotr12_2D;
  kokkos_athread_hash_table_insert (func55);
  AthreadHashNode* func56 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func56->key = arr_char_to_arr_int (
      "FunctorBarotr132");
  func56->fp = FunctorBarotr13_2D;
  kokkos_athread_hash_table_insert (func56);
  AthreadHashNode* func57 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func57->key = arr_char_to_arr_int (
      "FunctorBarotr142");
  func57->fp = FunctorBarotr14_2D;
  kokkos_athread_hash_table_insert (func57);
  AthreadHashNode* func58 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func58->key = arr_char_to_arr_int (
      "FunctorBarotr152");
  func58->fp = FunctorBarotr15_2D;
  kokkos_athread_hash_table_insert (func58);
  AthreadHashNode* func59 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func59->key = arr_char_to_arr_int (
      "FunctorBarotr162");
  func59->fp = FunctorBarotr16_2D;
  kokkos_athread_hash_table_insert (func59);
  AthreadHashNode* func60 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func60->key = arr_char_to_arr_int (
      "FunctorBarotr172");
  func60->fp = FunctorBarotr17_2D;
  kokkos_athread_hash_table_insert (func60);
  AthreadHashNode* func61 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func61->key = arr_char_to_arr_int (
      "FunctorBarotr182");
  func61->fp = FunctorBarotr18_2D;
  kokkos_athread_hash_table_insert (func61);
  AthreadHashNode* func62 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func62->key = arr_char_to_arr_int (
      "FunctorBarotr192");
  func62->fp = FunctorBarotr19_2D;
  kokkos_athread_hash_table_insert (func62);
  AthreadHashNode* func63 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func63->key = arr_char_to_arr_int (
      "FunctorBarotr202");
  func63->fp = FunctorBarotr20_2D;
  kokkos_athread_hash_table_insert (func63);
  AthreadHashNode* func64 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func64->key = arr_char_to_arr_int (
      "FunctorBarotr212");
  func64->fp = FunctorBarotr21_2D;
  kokkos_athread_hash_table_insert (func64);
  AthreadHashNode* func65 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func65->key = arr_char_to_arr_int (
      "FunctorReadyt13");
  func65->fp = FunctorReadyt1_3D;
  kokkos_athread_hash_table_insert (func65);
  AthreadHashNode* func66 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func66->key = arr_char_to_arr_int (
      "FunctorReadyt23");
  func66->fp = FunctorReadyt2_3D;
  kokkos_athread_hash_table_insert (func66);
  AthreadHashNode* func67 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func67->key = arr_char_to_arr_int (
      "FunctorReadyt33");
  func67->fp = FunctorReadyt3_3D;
  kokkos_athread_hash_table_insert (func67);
  AthreadHashNode* func68 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func68->key = arr_char_to_arr_int (
      "FunctorReadyt53");
  func68->fp = FunctorReadyt5_3D;
  kokkos_athread_hash_table_insert (func68);
  AthreadHashNode* func69 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func69->key = arr_char_to_arr_int (
      "FunctorReadyt63");
  func69->fp = FunctorReadyt6_3D;
  kokkos_athread_hash_table_insert (func69);
  AthreadHashNode* func70 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func70->key = arr_char_to_arr_int (
      "FunctorReadyt73");
  func70->fp = FunctorReadyt7_3D;
  kokkos_athread_hash_table_insert (func70);
  AthreadHashNode* func71 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func71->key = arr_char_to_arr_int (
      "FunctorReadyt82");
  func71->fp = FunctorReadyt8_2D;
  kokkos_athread_hash_table_insert (func71);
  AthreadHashNode* func72 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func72->key = arr_char_to_arr_int (
      "FunctorReadyt102");
  func72->fp = FunctorReadyt10_2D;
  kokkos_athread_hash_table_insert (func72);
  AthreadHashNode* func73 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func73->key = arr_char_to_arr_int (
      "FunctorReadyt112");
  func73->fp = FunctorReadyt11_2D;
  kokkos_athread_hash_table_insert (func73);
  AthreadHashNode* func74 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func74->key = arr_char_to_arr_int (
      "FunctorReadyt122");
  func74->fp = FunctorReadyt12_2D;
  kokkos_athread_hash_table_insert (func74);
  AthreadHashNode* func75 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func75->key = arr_char_to_arr_int (
      "FunctorReadyt132");
  func75->fp = FunctorReadyt13_2D;
  kokkos_athread_hash_table_insert (func75);
  AthreadHashNode* func76 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func76->key = arr_char_to_arr_int (
      "FunctorInterplationNearest2");
  func76->fp = FunctorInterplationNearest_2D;
  kokkos_athread_hash_table_insert (func76);
  AthreadHashNode* func77 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func77->key = arr_char_to_arr_int (
      "FunctorNear12");
  func77->fp = FunctorNear1_2D;
  kokkos_athread_hash_table_insert (func77);
  AthreadHashNode* func78 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func78->key = arr_char_to_arr_int (
      "FunctorNear21");
  func78->fp = FunctorNear2_1D;
  kokkos_athread_hash_table_insert (func78);
  AthreadHashNode* func79 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func79->key = arr_char_to_arr_int (
      "FunctorJRADaily12");
  func79->fp = FunctorJRADaily1_2D;
  kokkos_athread_hash_table_insert (func79);
  AthreadHashNode* func80 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func80->key = arr_char_to_arr_int (
      "FunctorJRADaily22");
  func80->fp = FunctorJRADaily2_2D;
  kokkos_athread_hash_table_insert (func80);
  AthreadHashNode* func81 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func81->key = arr_char_to_arr_int (
      "FunctorJRADaily32");
  func81->fp = FunctorJRADaily3_2D;
  kokkos_athread_hash_table_insert (func81);
  AthreadHashNode* func82 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func82->key = arr_char_to_arr_int (
      "FunctorNcarOceanFluxesJra2");
  func82->fp = FunctorNcarOceanFluxesJra_2D;
  kokkos_athread_hash_table_insert (func82);
  AthreadHashNode* func83 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func83->key = arr_char_to_arr_int (
      "FunctorNextStep13");
  func83->fp = FunctorNextStep1_3D;
  kokkos_athread_hash_table_insert (func83);
  AthreadHashNode* func84 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func84->key = arr_char_to_arr_int (
      "FunctorNextStep22");
  func84->fp = FunctorNextStep2_2D;
  kokkos_athread_hash_table_insert (func84);
  AthreadHashNode* func85 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func85->key = arr_char_to_arr_int (
      "FunctorTracer12");
  func85->fp = FunctorTracer1_2D;
  kokkos_athread_hash_table_insert (func85);
  AthreadHashNode* func86 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func86->key = arr_char_to_arr_int (
      "FunctorTracer23");
  func86->fp = FunctorTracer2_3D;
  kokkos_athread_hash_table_insert (func86);
  AthreadHashNode* func87 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func87->key = arr_char_to_arr_int (
      "FunctorTracer32");
  func87->fp = FunctorTracer3_2D;
  kokkos_athread_hash_table_insert (func87);
  AthreadHashNode* func88 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func88->key = arr_char_to_arr_int (
      "FunctorTracer43");
  func88->fp = FunctorTracer4_3D;
  kokkos_athread_hash_table_insert (func88);
  AthreadHashNode* func89 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func89->key = arr_char_to_arr_int (
      "FunctorTracer53");
  func89->fp = FunctorTracer5_3D;
  kokkos_athread_hash_table_insert (func89);
  AthreadHashNode* func90 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func90->key = arr_char_to_arr_int (
      "FunctorTracer73");
  func90->fp = FunctorTracer7_3D;
  kokkos_athread_hash_table_insert (func90);
  AthreadHashNode* func91 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func91->key = arr_char_to_arr_int (
      "FunctorTracer83");
  func91->fp = FunctorTracer8_3D;
  kokkos_athread_hash_table_insert (func91);
  AthreadHashNode* func92 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func92->key = arr_char_to_arr_int (
      "FunctorTracer153");
  func92->fp = FunctorTracer15_3D;
  kokkos_athread_hash_table_insert (func92);
  AthreadHashNode* func93 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func93->key = arr_char_to_arr_int (
      "FunctorTracer163");
  func93->fp = FunctorTracer16_3D;
  kokkos_athread_hash_table_insert (func93);
  AthreadHashNode* func94 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func94->key = arr_char_to_arr_int (
      "FunctorTracer173");
  func94->fp = FunctorTracer17_3D;
  kokkos_athread_hash_table_insert (func94);
  AthreadHashNode* func95 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func95->key = arr_char_to_arr_int (
      "FunctorTracer183");
  func95->fp = FunctorTracer18_3D;
  kokkos_athread_hash_table_insert (func95);
  AthreadHashNode* func96 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func96->key = arr_char_to_arr_int (
      "FunctorTracer192");
  func96->fp = FunctorTracer19_2D;
  kokkos_athread_hash_table_insert (func96);
  AthreadHashNode* func97 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func97->key = arr_char_to_arr_int (
      "FunctorTracer203");
  func97->fp = FunctorTracer20_3D;
  kokkos_athread_hash_table_insert (func97);
  AthreadHashNode* func98 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func98->key = arr_char_to_arr_int (
      "FunctorTracer212");
  func98->fp = FunctorTracer21_2D;
  kokkos_athread_hash_table_insert (func98);
  AthreadHashNode* func99 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func99->key = arr_char_to_arr_int (
      "FunctorTracer222");
  func99->fp = FunctorTracer22_2D;
  kokkos_athread_hash_table_insert (func99);
  AthreadHashNode* func100 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func100->key = arr_char_to_arr_int (
      "FunctorTracer233");
  func100->fp = FunctorTracer23_3D;
  kokkos_athread_hash_table_insert (func100);
  AthreadHashNode* func101 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func101->key = arr_char_to_arr_int (
      "FunctorTracer242");
  func101->fp = FunctorTracer24_2D;
  kokkos_athread_hash_table_insert (func101);
  AthreadHashNode* func102 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func102->key = arr_char_to_arr_int (
      "FunctorTracer252");
  func102->fp = FunctorTracer25_2D;
  kokkos_athread_hash_table_insert (func102);
  AthreadHashNode* func103 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func103->key = arr_char_to_arr_int (
      "FunctorTracer272");
  func103->fp = FunctorTracer27_2D;
  kokkos_athread_hash_table_insert (func103);
  AthreadHashNode* func104 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func104->key = arr_char_to_arr_int (
      "FunctorTracer282");
  func104->fp = FunctorTracer28_2D;
  kokkos_athread_hash_table_insert (func104);
  AthreadHashNode* func105 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func105->key = arr_char_to_arr_int (
      "FunctorTracer292");
  func105->fp = FunctorTracer29_2D;
  kokkos_athread_hash_table_insert (func105);
  AthreadHashNode* func106 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func106->key = arr_char_to_arr_int (
      "FunctorTracer303");
  func106->fp = FunctorTracer30_3D;
  kokkos_athread_hash_table_insert (func106);
  AthreadHashNode* func107 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func107->key = arr_char_to_arr_int (
      "FunctorTracer313");
  func107->fp = FunctorTracer31_3D;
  kokkos_athread_hash_table_insert (func107);
  AthreadHashNode* func108 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func108->key = arr_char_to_arr_int (
      "FunctorTracer322");
  func108->fp = FunctorTracer32_2D;
  kokkos_athread_hash_table_insert (func108);
  AthreadHashNode* func109 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func109->key = arr_char_to_arr_int (
      "FunctorTracer332");
  func109->fp = FunctorTracer33_2D;
  kokkos_athread_hash_table_insert (func109);
  AthreadHashNode* func110 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func110->key = arr_char_to_arr_int (
      "FunctorTracer343");
  func110->fp = FunctorTracer34_3D;
  kokkos_athread_hash_table_insert (func110);
  AthreadHashNode* func111 = (AthreadHashNode*)
      libc_aligned_malloc(sizeof(AthreadHashNode));
  func111->key = arr_char_to_arr_int (
      "FunctorTracer353");
  func111->fp = FunctorTracer35_3D;
  kokkos_athread_hash_table_insert (func111);
  // __REGISTER_END__
  
  /*
  if (athread_tid == 0) {
    kokkos_athread_hash_table_print (); 
  }
  */
  return ;
}

inline void check_function_pointer (const void* fp, const char* func_name) {
  if (fp == nullptr) {
    if (athread_tid == 0) {
      // printf ("Error in \"%s\" -> \"", __PRETTY_FUNCTION__);
      printf ("Error in \"%s\" -> \"", func_name);
      for (int i = 0; i < KOKKOS_ATHREAD_KEY_LEN; ++i) {
        printf ("%c", g_athread_functor_key[i]);
      }
      printf ("\". Hash value: %d\n", kokkos_athread_local_param.hash_value);
    }
    athread_ssync_array(); 
    exit (EXIT_FAILURE);
  }
  return ;
}

extern "C" void kokkos_athread_parallel_for (AthreadParamWrap* host_param) {
  athread_dma_iget(&kokkos_athread_local_param, host_param, 
      sizeof(AthreadParamWrap), &dma_rply);
  D_COUNT++;
  athread_dma_wait_value(&dma_rply, D_COUNT);
  const FunctionPointer fp __attribute__ ((aligned(64))) = 
      kokkos_athread_hash_table_lookup(&kokkos_athread_local_param)->fp;

  check_function_pointer ((void*)fp, __PRETTY_FUNCTION__);

  fp(&kokkos_athread_local_param, athread_tid);
  // if (athread_tid == 0) {
  //   printf("launch fp: %p\n", fp);
  // }
  return ;
}

extern "C" void kokkos_athread_parallel_reduce (AthreadParamWrap* host_param) {
  athread_dma_iget(&kokkos_athread_local_param, host_param, 
      sizeof(AthreadParamWrap), &dma_rply);
  D_COUNT++;
  athread_dma_wait_value(&dma_rply, D_COUNT);
  const FunctionPointer fp __attribute__ ((aligned(64))) = 
      kokkos_athread_hash_table_lookup(&kokkos_athread_local_param)->fp;
  check_function_pointer ((void*)fp, __PRETTY_FUNCTION__);
  fp(&kokkos_athread_local_param, athread_tid);
  // reduce from all slave cores
  slave_cores_reduce(&(kokkos_athread_local_param.reduce_double),
      &(kokkos_athread_local_param.reduce_double), 1, athread_double, OP_add, &buf_reduce, 64);
  // slave_cores_reduce(&(kokkos_athread_local_param.reduce_int),
  //     &(kokkos_athread_local_param.reduce_int), 1, athread_int, OP_and, &buf_reduce_int, 1);
  if (athread_tid == 0) {
    // host_param->reduce_int    = kokkos_athread_local_param.reduce_int;
    host_param->reduce_double = kokkos_athread_local_param.reduce_double;
  }
  return ;
}

// extern "C" void parallel_reduce_1D(AthreadParamWrap* host_param) {

//   athread_dma_iget(&kokkos_athread_local_param, host_param, sizeof(AthreadParamWrap), &dma_rply);
//   D_COUNT++;
//   athread_dma_wait_value(&dma_rply, D_COUNT);

//   const FunctionPointer fp __attribute__ ((aligned(64))) = lookup_fp();

//   const int ATHREAD_SLAVE_CORES = 64;

//   const int start = kokkos_athread_local_param.range[0][0];
//   const int end   = kokkos_athread_local_param.range[0][1];
//   const int len   = end - start;
//   const int times = (len + ATHREAD_SLAVE_CORES - 1) / ATHREAD_SLAVE_CORES;

//   for (int i = 0; i < times; ++i) {
//     int index = start + athread_tid * times + i;
//     if (index < end) {
//       kokkos_athread_local_param.index[0] = index;
//       fp(&kokkos_athread_local_param);
//     }
//   }
//   // reduce from all slave cores
//   slave_cores_reduce(&(kokkos_athread_local_param.reduce_result),
//       &(kokkos_athread_local_param.reduce_result), 1, athread_double, OP_add, &buf_reduce, 1);
//   if (athread_tid == 0) {
//     host_param->reduce_result = kokkos_athread_local_param.reduce_result;
//   }
//   return ;
// }

// extern "C" void parallel_reduce_2D(AthreadParamWrap* host_param) {

//   athread_dma_iget(&kokkos_athread_local_param, host_param, sizeof(AthreadParamWrap), &dma_rply);
//   D_COUNT++;
//   athread_dma_wait_value(&dma_rply, D_COUNT);

//   const FunctionPointer fp __attribute__ ((aligned(64))) = lookup_fp();

//   const int ATHREAD_SLAVE_CORES = 64;

//   const int start0 = kokkos_athread_local_param.range[0][0];
//   const int end0   = kokkos_athread_local_param.range[0][1];
//   const int start1 = kokkos_athread_local_param.range[1][0];
//   const int end1   = kokkos_athread_local_param.range[1][1];
//   const int tile0  = kokkos_athread_local_param.tile[0];
//   const int tile1  = kokkos_athread_local_param.tile[1];

//   const int num_tiles0 = (end0 - start0 + tile0 - 1) / tile0;
//   const int num_tiles1 = (end1 - start1 + tile1 - 1) / tile1;

//   const int num_tiles = num_tiles0 * num_tiles1;
  
//   for (int index_tile = athread_tid; index_tile < num_tiles; 
//       index_tile += ATHREAD_SLAVE_CORES) {
//     const int index_tile0 = index_tile / num_tiles1;
//     const int index_tile1 = index_tile % num_tiles1;

//     for (int index00 = 0; index00 < tile0; ++ index00) {
//       const int index0 = start0 + index_tile0 * tile0 + index00;
//       if (index0 >= end0) { break; }
//       for (int index_11 = 0; index_11 < tile1; ++ index_11) {
//         const int index1 = start1 + index_tile1 * tile1 + index_11;
//         if (index1 >= end1) { break; }
//         kokkos_athread_local_param.index[0] = index0;
//         kokkos_athread_local_param.index[1] = index1;
//         fp(&kokkos_athread_local_param);
//       }
//     }
//   }
//   // reduce from all slave cores
//   slave_cores_reduce(&(kokkos_athread_local_param.reduce_result),
//       &(kokkos_athread_local_param.reduce_result), 1, athread_double, OP_add, &buf_reduce, 1);
//   if (athread_tid == 0) {
//     host_param->reduce_result = kokkos_athread_local_param.reduce_result;
//   }
//   return ;
// }

// extern "C" void parallel_reduce_3D(AthreadParamWrap* host_param) {

//   athread_dma_iget(&kokkos_athread_local_param, host_param, sizeof(AthreadParamWrap), &dma_rply);
//   D_COUNT++;
//   athread_dma_wait_value(&dma_rply, D_COUNT);

//   const FunctionPointer fp __attribute__ ((aligned(64))) = lookup_fp();

//   const int ATHREAD_SLAVE_CORES = 64;

//   const int start[3] = {kokkos_athread_local_param.range[0][0], kokkos_athread_local_param.range[1][0], kokkos_athread_local_param.range[2][0]};
//   const int end[3]   = {kokkos_athread_local_param.range[0][1], kokkos_athread_local_param.range[1][1], kokkos_athread_local_param.range[2][1]};
//   const int tile[3]  = {kokkos_athread_local_param.tile[0],     kokkos_athread_local_param.tile[1],     kokkos_athread_local_param.tile[2]};

//   const int num_tiles0 = (end[0] - start[0] + tile[0] - 1) / tile[0];
//   const int num_tiles1 = (end[1] - start[1] + tile[1] - 1) / tile[1];
//   const int num_tiles2 = (end[2] - start[2] + tile[2] - 1) / tile[2];

//   const int tmp0 = num_tiles1 * num_tiles2;
//   const int num_tiles = num_tiles0 * tmp0;
  
//   for (int index_tile = athread_tid; index_tile < num_tiles; index_tile += ATHREAD_SLAVE_CORES) {

//     const int index_tile0 = index_tile / tmp0;
//     const int tmp1        = index_tile % tmp0;
//     const int index_tile1 = tmp1      / num_tiles2;
//     const int index_tile2 = tmp1      % num_tiles2;

//     for (int index_00 = 0; index_00 < tile[0]; ++ index_00) {
//       const int index0 = start[0] + index_tile0 * tile[0] + index_00;
//       if (index0 >= end[0]) { break; }
//       for (int index_11 = 0; index_11 < tile[1]; ++ index_11) {
//         const int index1 = start[1] + index_tile1 * tile[1] + index_11;
//         if (index1 >= end[1]) { break; }
//         for (int index_22 = 0; index_22 < tile[2]; ++ index_22) {
//           const int index2 = start[2] + index_tile2 * tile[2] + index_22;
//           if (index2 >= end[2]) { break; }

//           kokkos_athread_local_param.index[0] = index0;
//           kokkos_athread_local_param.index[1] = index1;
//           kokkos_athread_local_param.index[2] = index2;
//           fp(&kokkos_athread_local_param);
//         }
//       }
//     }
//   }
//   // reduce from all slave cores
//   slave_cores_reduce(&(kokkos_athread_local_param.reduce_result),
//       &(kokkos_athread_local_param.reduce_result), 1, athread_double, OP_add, &buf_reduce, 1);
//   if (athread_tid == 0) {
//     host_param->reduce_result = kokkos_athread_local_param.reduce_result;
//   }
//   return ;
// }
