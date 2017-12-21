# PsimagLite support is needed by DMRG++
LDFLAGS = -L../../PsimagLite/lib -L../../../PsimagLite/lib -lpsimaglite

# Compiler to use. If using MPI then say mpicxx here (or mpic++)
# and also say -DUSE_MPI below
CXX = g++

# We're using ansi C++
CPPFLAGS += -pedantic -std=c++98

# Enable MPI (you must set the proper
# compiler wrapper under CXX above)
# CPPFLAGS += -DUSE_MPI

# Enable warnings and treat warnings as errors
CPPFLAGS += -Wall -Werror -Wendif-labels

# This enables additional debugging
#CPPFLAGS += -D_GLIBCXX_DEBUG -D_GLIBCXX_PROFILE

# This makes the code use long instead of short integers
#CPPFLAGS +=-DUSE_LONG

# This makes the code use float instead of double
#CPPFLAGS += -DUSE_FLOAT

# This enables signals
#CPPFLAGS +=-DUSE_SIGNALS

# This enables gsl support
#CPPFLAGS +=-DUSE_GSL
#LDFLAGS += -lgsl -lgslcblas

# This enables the custom allocator (use only for debugging)
#CPPFLAGS += -DUSE_CUSTOM_ALLOCATOR

#Change basis even for un-needed operators
#CPPFLAGS += -DOPERATORS_CHANGE_ALL

# Disable KronUtil
#CPPFLAGS += -DDO_NOT_USE_KRON_UTIL

#Add directory to linker where libkronutil.a resides
LDFLAGS += -LKronUtil

# Specify the strip command to use (or use true to disable)
STRIP_COMMAND = true

#This enables boost support that is needed for Ainur
#CPPFLAGS += -DUSE_BOOST

#This enables the PLUGIN_SC repository
#CPPFLAGS += -DPLUGIN_SC

#This adds linkage for the PLUGIN_SC libraries
#LDFLAGS += -L../../dmrgppPluginSc/src -L../../../dmrgppPluginSc/src
#LDFLAGS +=  -ldmrgppPluginSc
#LDFLAGS +=  -lgomp

