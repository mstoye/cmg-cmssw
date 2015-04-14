#!/bin/zsh
## make sure the right shell will be used
#$ -S /bin/zsh
## the cpu time for this job
#$ -l h_rt=02:59:00
## the maximum memory usage of this job
#$ -l h_vmem=1900M
## operating system
#$ -l distro=sld6
## architecture
##$ -l arch=amd64
## Name
#$ -N CMGrunner
## stderr and stdout are merged together to stdout
#$ -j y
##(send mail on job's end and abort)
##$ -m a
#$ -l site=hh
## transfer env var from submission host
#$ -V
## set cwd to submission host pwd
#$ -cwd

TaskID=$((SGE_TASK_ID))
echo "SGE_TASK_ID: " $TaskID

JobList=$1
TaskCmd=$(cat $JobList | sed ''$TaskID'q;d')

#eval $(scramv1 runtime -sh);

echo "Going to execute" $TaskCmd
eval $TaskCmd
