#!/bin/bash

# We assume running this from the script directory
# Run as ./gxtb.sh gaussianinputfile.com 
# Where gaussianinputfile.com does not contain any % lines before the route section (#)
# Will generate or attempt to generate both gaussianinputfile.log in the same directory and temporary deletable files
# with the same filename root.


export PATH=$(pwd):/home/laplaza/Software/xtb/_build:$PATH
function is_bin_in_path {
     builtin type -P "$1" &> /dev/null
}
module load gaussian/g16/A.03
is_bin_in_path xtb  && echo "Found xtb." || echo "No xtb found. Exit!"  
is_bin_in_path crest  && echo "Found crest." || echo "No crest found. Exit!" 
is_bin_in_path g16  && echo "Found g16." || echo "No g16 found. Exit!" 
job_directory=$(pwd)
input=${1}
name="${1%%.com}"
inpname="${name}.com"
inpname_xtb="${name}_xtb.com"
outname="${name}.log"
output="${job_directory}/${name}.out"
curdir='$SLURM_SUBMIT_DIR'
gauroute=$(grep "#" ${inpname}|  tr [:upper:] [:lower:] | grep "opt" )
if [[ ${gauroute} =~ "ts" ]]; then
    xtbroute='# IOP(1/18=0) IOP(1/6=500) external="xtb-gaussian -P 8"'
else
    xtbroute='# IOP(1/6=500) external="xtb-gaussian -P 8"'
fi
echo ${xtbroute} | cat - ${inpname} > ${inpname_xtb}
echo "#!/bin/bash
#SBATCH -J ${name}
#SBATCH -o ${output}
#SBATCH --mem=8000
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=3-00:00:00
module load intel intel-mkl
ulimit -s unlimited
export OMP_STACKSIZE=4G
export OMP_MAX_ACTIVE_LEVELS=1
export RAMDIR=/dev/shm/\${SLURM_JOBID}
mkdir -p \${RAMDIR}
cd \${RAMDIR}
export GAUSS_SCRDIR=\${RAMDIR}
export SLURM_SCRDIR=\${RAMDIR}
cp ${curdir}/${inpname_xtb} \${RAMDIR}
#echo 'Running in '\${RAMDIR}' and output to '${outname}' in '${curdir}
g16 < ${inpname_xtb} > ${outname}
cp ${outname} ${curdir}
rm -rf \${RAMDIR}
exit " > ${name}.job
sbatch ${name}.job


