#!/bin/sh
#bsub -I -b -q q_test_big -pr -shared -n 2 -cgsp 1 -host_stack 256 -share_size 4096 -cache_size 128 ./vinardo --config config.txt --ligand 1opj_ligand.pdbqt 2>&1 | tee 0123_0.log
#bsub -I -b -swrunarg "-P slave" -q q_test_ss -pr -shared -n 2 -cgsp 64 -host_stack 256 -share_size 4096 -cache_size 128 ./vinardo --config config.txt --ligand 1opj_ligand.pdbqt
#bsub -I -b -q q_test_debug -pr -shared -n 2 -cgsp 64 -host_stack 256 -share_size 4096 -cache_size 128 ./vinardo --config config.txt --ligand 1opj_ligand.pdbqt
#bsub -I -b -perf_float -q q_test_debug -pr -shared -n 2 -cgsp 64 -host_stack 1024 -share_size 4096 -cache_size 128 ./vinardo --config config.txt --ligand 1opj_ligand.pdbq
#bsub -o test_0623.log -I -b -q q_test_ss -pr -n 3 -cgsp 64 -host_stack 4096 -share_size 11000 -cache_size 128 ./vinardo --config vina.in --ligand 1opj_ligand.pdbqt --manageNum 500
#bsub -o test_drugbank_0605.log -I -b -q q_test_ss -pr -n 4001 -cgsp 64 -host_stack 4096 -share_size 11000 -cache_size 128 ./vinardo --config vina.in --ligand 1opj_ligand.pdbqt
bsub -o test1223.log -I -b -q q_sw_expr -pr -n 3 -cgsp 64 -host_stack 4096 -share_size 11000 -cache_size 128 ./vina_reserve --manageNum 1
