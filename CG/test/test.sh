exe="bsub -J wjl-CG -b -cache_size 128 -q q_share -shared -n 1 -cgsp 64 -host_stack 1024 -share_size 15000"
input="cg_matrix_set.csv"
{
  read
  i=1
  while IFS=',' read -r data a b each_nnz
  do
    for data1 in `find "./matrix" -name "$data.mtx"`
    do
      echo -n $data1
       # 执行命令时，重定向 /dev/null 防止其消耗主输入流
        $exe ./CG  $data1  < /dev/null
        # 间隔1秒
      # sleep 1
      # timeout -s 9 4m ./cg_cusparse $data1 
    done
    i=`expr $i + 1`
  done 
} < "$input"

# input="bicg_matrix_set.csv"
# {
#   read
#   i=1
#   while IFS=',' read -r data a b each_nnz
#   do
#     for data1 in `find "./matrix" -name "$data.mtx"`
#     do
#       echo -n $data1
#        # 执行命令时，重定向 /dev/null 防止其消耗主输入流
#        $exe ./bicgstab_solve  $data1   < /dev/null
#       # timeout -s 9 4m ./bicg_cusparse $data1 
#     done
#     i=`expr $i + 1`
#   done 
# } < "$input"
