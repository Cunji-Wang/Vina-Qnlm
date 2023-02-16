BASE=/home/export/online1/mdt00/shisuan/swyaow/luhao/lh/boost/boost_1_67_0
BOOST_VERSION=1_67_0
#BOOST_INCLUDE =-I $(BASE)include
BOOST_INCLUDE =/home/export/online1/mdt00/shisuan/swyaow/luhao/lh/boost/boost_1_67_0/include
#BOOST_INCLUDE =/home/export/base/qnlm/qnlm01/yumaoxue/boost_1_73_0/boost

#MPI_INCLUDE = /usr/sw/mpi/mpisw9gcc_20201020/include
MPI_INCLUDE = /usr/sw/mpi/mpi_current/include

ATHREAD_INCLUDE = 

LIBS=/home/export/online1/mdt00/shisuan/swyaow/luhao/lh/boost/boost_1_67_0/lib
#BOOST_LIB=-l/home/export/base/qnlm/qnlm01/yumaoxue/boost_1_73_0/libs
CXX=mpic++
HOST_CXX=sw9g++
SLAVE_C=sw9gcc

CC_SLAVE_FLAGS := -mslave -mieee -msimd -mftz -O3
#-fdisable_slave_ioproxy
C_HYBRID_FLAGS := -mhybrid -static
C_PLATFORM=-static


C_OPTIONS= -O3 -DNDEBUG -DATHREAD_MC 
CFLAGS=  -mhost  -mftz -mieee #-OPT:IEEE_arith=2 
BOOST_LIB_VERSION=

include ../../makefile_common


