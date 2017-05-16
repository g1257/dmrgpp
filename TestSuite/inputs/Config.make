# Base Config
include ConfigBase.make

# Here add your lapack and blas libraries or say NO_LAPACK
# CPPFLAGS += -DNO_LAPACK
# If on MacOs please say LDFLAGS += -framework Accelerate
LDFLAGS += -llapack -lblas

# Here add -lpthread if threading is needed and also 
# set -DUSE_PTHREADS below
LDFLAGS += -lpthread

# Enable pthreads
CPPFLAGS += -DUSE_PTHREADS

# This disables debugging
CPPFLAGS += -DNDEBUG

# Optimization level here
CPPFLAGS += -O3

# Specify the strip command to use (or use true to disable) 
STRIP_COMMAND = strip 
