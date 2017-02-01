# Base Config
include ConfigBase.make

# Here add -lpthread if threading is needed and also 
# set -DUSE_PTHREADS below
LDFLAGS += -lpthread

# Enable pthreads
CPPFLAGS += -DUSE_PTHREADS

# This enables partial debugging (make sure to comment out previous line)
CPPFLAGS +=   -g3

# Here add your lapack and blas libraries or say NO_LAPACK
# CPPFLAGS += -DNO_LAPACK
# If on MacOs please say LDFLAGS += -framework Accelerate
LDFLAGS += /usr/lib64/liblapack.so.3 /usr/lib64/libblas.so.3
