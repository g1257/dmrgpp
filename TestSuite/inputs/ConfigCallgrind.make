# Base Config
include ConfigBase.make

# This disables debugging
CPPFLAGS += -DNDEBUG

# Optimization level here
CPPFLAGS += -O3

# This enables partial debugging (make sure to comment out previous line)
CPPFLAGS +=   -g3

# Here add your lapack and blas libraries or say NO_LAPACK
# CPPFLAGS += -DNO_LAPACK
# If on MacOs please say LDFLAGS += -framework Accelerate
LDFLAGS += -lblas -llapack

