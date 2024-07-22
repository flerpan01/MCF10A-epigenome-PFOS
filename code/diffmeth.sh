#!/bin/bash

#SBATCH -A 
#SBATCH -p node -C mem256GB -N 2
#SBATCH -t 8:00:00

#SBATCH --mail-type=ALL
#SBATCH -J diffmeth

PROJDIR=/home/$USER/proj/PFOS/em_seq/
cd $PROJDIR

ml R_packages/4.1.1

DOSE=(1)
TILE=(NULL 100)

for dose in ${DOSE[@]}; do

	for tile in ${TILE[@]}; do

		if [[ ${tile} == NULL ]]; then
  
		  VAR="cpg"
		  echo running script for invidivual CpG-sites

		elif [[ ${tile} =~ 0 ]]; then
  
	  	VAR="tile"${tile}
	  	echo running script for tiles of ${tile}bp

	  fi

	  R CMD BATCH --no-save --no-restore \
			'--args dose='${dose}' tile='${tile}'' \
			code/diffmeth.R \
			code/logs/diffmeth_${dose}_${VAR}.Rout &

		echo ${dose}_${VAR}

	done
done

wait
