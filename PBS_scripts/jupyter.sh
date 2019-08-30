#PBS -q batch
#PBS -l nodes=1:ppn=8
#PBS -l walltime=8:00:00

module load Anaconda/4.2.0
source activate scSplit 

PORT=$(shuf -i10000-11999 -n1)
echo "executing jupyter on http://$(hostname):$PORT"
 
PATH=/home/danaco/.conda/envs/scSplit/bin:$PATH


jupyter-notebook --no-browser --port=$PORT --ip=`hostname -i`
/home/danaco/.conda/envs/scSplit/bin/jupyter notebook list
