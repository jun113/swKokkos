CUR_DATE=$(date +'%Y%m%d-%H%M%S')
outName=log/${CUR_DATE}.out
# bsub -p -b -o debug.out -q q_share -J wjl_debug -pr -n 1 -cgsp 64 -ro_size 200 -host_stack 512 -share_size 12288 -cross_size 512 -priv_size 32 -cache_size 128 ./debug_kokkos

bsub -J wjl_swKokkosTest \
	-I -b -cache_size 128 \
	-q q_share \
	-n 1 \
	-cgsp 64 \
        -mpecg 1 \
	-share_size 15360 \
	-host_stack 2048 \
	-o ${outName} ./swKokkosTest | tee ${outName}

	# -o ${outName} ./debug_kokkos 20 4096000 1 0 | tee ${outName}
	# -o ${outName} ./debug_kokkos 15 4048000 1 0 | tee ${outName}
	# -o ${outName} ./debug_kokkos 25 10000000 1 0 | tee ${outName}
	# 4D argc = 9
	# A B RUNTIMES PRINTMATRIX
