#!/usr/bin/env bash
source activate py36
#for source in setup_*.py
for src in setup_E_FF.py setup_poynting.py setup_Z_CFIE.py setup_nbin.py setup_trans.py setup_Z_obs.py
do
	python $src build_ext -i
done
