#!/bin/bash
#SBATCH -J 11.2_phy	#任务名称
#SBATCH -p free		#队列名
#SBATCH -N 1			#请求的节点数(正常设置为1)
#SBATCH --no-requeue		#若任务失败，不再尝试重新排队
#SBATCH --ntasks-per-node=10	#线程数(不得超过20)
#SBATCH --output=zzz.slurm-%j.out	#日志文件设置(一般无需变动)
#SBATCH --error=zzz.slurm-%j.stderr #错误文件设置(一般无需变动)
#SBATCH --mail-type=END		#通知邮件类型
#SBATCH --mail-user=tianpew@uci.edu    #通知邮件地址


###############

## software
PATH=/share/crsp/lab/bgaut/tianpew/miniforge3/envs/evo/bin:$PATH
#cd $SLURM_SUBMIT_DIR		#进入提交任务时的工作目录
###############

iqtree -s sample97_pt_96.seq.aln.mafft -st DNA -nt 10 -m TVM+R6 -B 1000 -pre sample97_pt_96.seq.aln.mafft.tree

