# Vina@Qnlm
# =====================================

# About

Vina@Qnlm is a very efficient, accurate and convenient parallel version of molecular docking software.

Vina@qnlm is an improved ultra-massive parallel virtual screening method based on AutoDock Vina on the Sunway TaihuLight domestic supercomputer system. 

We have reconstructed its optimized search algorithm so that it can utilize as many computational resources as possible. 

In terms of performance, it adopts a two-stage parallel structure with powerful molecular docking parallel computing capability, which can achieve billion docking times in a short time. 

In terms of accuracy, the software performs well on the CASF-2016 benchmark test. 

Welcome to use the software and give valuable suggestions, thanks!

# Setup

| Operating system       | Heterogeneous processor        | memory             |
|:----------------------:|:------------------------------:|:------------------:|
|CentOS 6.7 & 6.8        | SW26010 1.45GHz 4MPE+4*64CPEs  | 32GB DDR3          |

# Quick Usage

cd build/linux/release

make 

vim bsub.sh

bash bsub.sh

# Note

Our team focuses on the research work in the direction of drug virtual computing, welcome to communicate and cooperate together.
