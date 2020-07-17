##########################################
# Environmen configurations
# Please change the variables, if necessary
##########################################

CUDA_HOME=/usr/local/cuda
CUTIL_HOME=~/NVIDIA_GPU_Computing_SDK

CUDA_INC=$CUDA_HOME/include
CUDA_LIB_DIR=$CUDA_HOME/lib64
CUTIL_INC=$CUTIL_HOME/C/common/inc
CUTIL_LIB_DIR=$CUTIL_HOME/C/lib
CUTIL_LIB_NAME=cutil_x86_64


##########################################
# cleanup
##########################################
echo "Deleting output files ..."
rm *.o
rm test.run


##########################################
# the file gqd_api.cu has included all 
# necessary library files
##########################################
echo "Compiling GQD library (may take several minutes)..." 
nvcc ./test/gqd_api.cu -c -o gqd_lib.o -I./src -I$CUTIL_INC -arch sm_13 --opencc-options -OPT:Olimit=0


##########################################
# the test case needs QD library
##########################################

echo "Compiling test cases (need the QD library and OpenMP support)..."
g++ -o test.run test/gqd_util.cpp test/test_main.cpp gqd_lib.o -lcudart -lcuda -l$CUTIL_LIB_NAME -L$CUDA_LIB_DIR -L$CUTIL_LIB_DIR -I./src -I$CUTIL_INC -I$CUDA_INC -lqd -O3 -fopenmp


##########################################
# run the test cases
##########################################
echo "Testing GQD library..."
./test.run


