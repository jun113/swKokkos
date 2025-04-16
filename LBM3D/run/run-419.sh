#!/bin/bash
#SBATCH -J density
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH -p debug
###SBATCH -t 00:10:00		#指定作业最大运行30分钟
#SBATCH --mem=90G		#占用节点全部内存
#SBATCH -N 1			#指定节点数
#SBATCH --ntasks-per-node=1	#指定每个节点的进程数
#SBATCH --gres=dcu:1		#指定每个节点使用4块DCU卡
#SBATCH --exclusive

echo "Start time: `date`"                               #显示开始时间
echo "SLURM_JOB_ID: $SLURM_JOB_ID"                      #显示作业号
echo "SLURM_NNODES: $SLURM_NNODES"                      #显示节点数
echo "SLURM_NTASKS: $SLURM_NTASKS"                      #显示总任务数
echo "SLURM_TASKS_PER_NODE: $SLURM_TASKS_PER_NODE"      #显示每节点的任务数
echo "SLURM_CPUS_PER_TASK: $SLURM_CPUS_PER_TASK"        #显示每个任务使用的CPU数量
echo "SLURM_JOB_PARTITION: $SLURM_JOB_PARTITION"        #显示队列/分区名称
echo "SLURM_SUBMIT_DIR:$SLURM_SUBMIT_DIR"               #显示提交作业目录的路径
echo "SLURM_NODELIST:$SLURM_NODELIST"                   #显示执行节点列表名称

module purge
module load compiler/intel/2017.5.239
module load mathlib/netcdf/4.4.1/intel
module load mpi/hpcx/2.7.4/intel-2017.5.239
module load compiler/rocm/3.9.1
module load compiler/gnu/gcc-9.3.0


##sbatch
./../build/bin/Density.exe

echo "End time: `date`"  #显示结束时间
