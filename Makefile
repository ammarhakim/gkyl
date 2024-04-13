# -*- makefile-gmake -*-

# Sample Makefile to use installed gkylzero library: copy and modify
# for your needs

ARCH_FLAGS ?= -march=native
CUDA_ARCH ?= 70
# Warning flags: -Wall -Wno-unused-variable -Wno-unused-function -Wno-missing-braces
LDFLAGS = 
PREFIX ?= ${HOME}/gkylsoft
INSTALL_PREFIX ?= ${PREFIX}

# determine OS we are running on
UNAME = $(shell uname)

# Default lapack include and libraries: we prefer linking to static library
LAPACK_INC = $(PREFIX)/OpenBLAS/include
LAPACK_LIB_DIR = $(PREFIX)/OpenBLAS/lib
LAPACK_LIB = -lopenblas

# SuperLU includes and librararies
SUPERLU_INC = $(PREFIX)/superlu/include
ifeq ($(UNAME_S),Linux)
	SUPERLU_LIB_DIR = $(PREFIX)/superlu/lib64
	SUPERLU_LIB = $(PREFIX)/superlu/lib64/libsuperlu.a
else
	SUPERLU_LIB_DIR = $(PREFIX)/superlu/lib
	SUPERLU_LIB = $(PREFIX)/superlu/lib/libsuperlu.a
endif

BUILD_DIR ?= build

# Include config.mak file (if it exists)
-include $(PREFIX)/gkylzero/share/config.mak

CFLAGS = -O3 -g -I.

G0_INC_DIR = ${PREFIX}/gkylzero/include
G0_LIB_DIR = ${PREFIX}/gkylzero/lib
G0_LIB = -lgkylzero

USING_NVCC =
NVCC_FLAGS =
CUDA_LIBS =
ifeq ($(CC), nvcc)
       USING_NVCC = yes
       CFLAGS = -O3 -g --forward-unknown-to-host-compiler --use_fast_math -ffast-math -MMD -MP -fPIC
       NVCC_FLAGS = -x cu -dc -arch=sm_${CUDA_ARCH} --compiler-options="-fPIC"
       LDFLAGS += -arch=sm_${CUDA_ARCH}
       ifdef CUDAMATH_LIBDIR
              CUDA_LIBS = -L${CUDAMATH_LIBDIR}
       else
              CUDA_LIBS =
       endif
       CUDA_LIBS += -lcublas -lcusparse -lcusolver
endif

G0_LIBS = ${G0_LIB} ${CUDA_LIBS} -lm -lpthread
G0_RPATH = -Wl,-rpath,${G0_LIB_DIR}

# determine OS we are running on
UNAME = $(shell uname)

# On OSX we should use Accelerate framework
ifeq ($(UNAME), Darwin)
	LAPACK_INC = . # dummy
	LAPACK_LIB = -framework Accelerate
	CFLAGS += -DGKYL_USING_FRAMEWORK_ACCELERATE
endif

# Read MPI paths and flags if needed 
USING_MPI =
MPI_INC_DIR = zero # dummy
MPI_LIB_DIR = .
ifeq (${USE_MPI}, 1)
	USING_MPI = yes
	MPI_INC_DIR = ${CONF_MPI_INC_DIR}
	MPI_LIB_DIR = ${CONF_MPI_LIB_DIR}
	MPI_LIBS = -lmpi
	CFLAGS += -DGKYL_HAVE_MPI
endif

# Read NCCL paths and flags if needed (needs MPI and NVCC)
USING_NCCL =
NCCL_INC_DIR = zero # dummy
NCCL_LIB_DIR = .
ifeq (${USE_NCCL}, 1)
ifdef USING_MPI
ifdef USING_NVCC
	USING_NCCL = yes
	NCCL_INC_DIR = ${CONF_NCCL_INC_DIR}
	NCCL_LIB_DIR = ${CONF_NCCL_LIB_DIR}
	NCCL_LIBS = -lnccl
	CFLAGS += -DGKYL_HAVE_NCCL
endif
endif
endif

# Read LUA paths and flags if needed 
USING_LUA =
LUA_INC_DIR = zero # dummy
LUA_LIB_DIR = .
ifeq (${USE_LUA}, 1)
	USING_LUA = yes
	LUA_INC_DIR = ${CONF_LUA_INC_DIR}
	LUA_LIB_DIR = ${CONF_LUA_LIB_DIR}
	LUA_LIBS = -l${CONF_LUA_LIB}
	CFLAGS += -DGKYL_HAVE_LUA
endif

G2_INCLUDES = -IComm -ILib
INCLUDES = -I${G0_INC_DIR} -I${LAPACK_INC} -I${SUPERLU_INC} -I${MPI_INC_DIR} -I${LUA_INC_DIR} -I${NCCL_INC_DIR} ${G2_INCLUDES}
LIB_DIRS = -L${LAPACK_LIB_DIR} -L${SUPERLU_LIB_DIR} -L${MPI_LIB_DIR} -L${LUA_LIB_DIR} -L${NCCL_LIB_DIR}
EXT_LIBS = ${LAPACK_LIB} ${SUPERLU_LIB} ${MPI_LIBS} ${LUA_LIBS} ${NCCL_LIBS} -lm -lpthread -ldl

# Directories containing source code
SRC_DIRS := Comm Lib

# List of source files
SRCS := $(shell find $(SRC_DIRS) -name '*.c')

OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

# main build target
all: gkyl ## Build Gkeyll executable (gkyl)

# Build commands for C source
$(BUILD_DIR)/%.c.o: %.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CFLAGS) $(ARCH_FLAGS) $(INCLUDES) -c $< -o $@

# Build commands for C source
$(BUILD_DIR)/sqlite3/sqlite3.c.o: sqlite3/sqlite3.c
	$(MKDIR_P) $(dir $@)
	$(CC)  -Wno-implicit-int-float-conversion -g -O -c $< -o $@

gkyl: ${OBJS} gkyl.c ## Build main Gkeyll executable (gkyl)
	${CC} ${CFLAGS} ${INCLUDES} gkyl.c $(OBJS) -o gkyl -L${G0_LIB_DIR} ${G0_RPATH} ${G0_LIBS} ${LIB_DIRS} ${EXT_LIBS}

install: all ## Install library and headers
# Construct install directories
	$(MKDIR_P) ${INSTALL_PREFIX}/gkyl/bin
	${MKDIR_P} ${INSTALL_PREFIX}/gkyl/lib
# Copy main executable
	cp -f gkyl ${INSTALL_PREFIX}/gkyl/bin/gkyl
# Copy Lua code from various directories
	${MKDIR_P} ${INSTALL_PREFIX}/gkyl/bin/xsys
	find xsys -name '*.lua' | xargs cp --parents -fv -t ${INSTALL_PREFIX}/gkyl/bin
#
	${MKDIR_P} ${INSTALL_PREFIX}/gkyl/bin/sci
	find sci -name '*.lua' | xargs cp --parents -fv -t ${INSTALL_PREFIX}/gkyl/bin
#
	${MKDIR_P} ${INSTALL_PREFIX}/gkyl/bin/Tool
	find Tool -name '*.lua' | xargs cp --parents -fv -t ${INSTALL_PREFIX}/gkyl/bin
#
	${MKDIR_P} ${INSTALL_PREFIX}/gkyl/bin/sqlite3
	find sqlite3 -name '*.lua' | xargs cp --parents -fv -t ${INSTALL_PREFIX}/gkyl/bin
#
	${MKDIR_P} ${INSTALL_PREFIX}/gkyl/bin/Lib
	find Lib -name '*.lua' | xargs cp --parents -fv -t ${INSTALL_PREFIX}/gkyl/bin
#
	${MKDIR_P} ${INSTALL_PREFIX}/gkyl/bin/Grid
	find Grid -name '*.lua' | xargs cp --parents -fv -t ${INSTALL_PREFIX}/gkyl/bin
#
	${MKDIR_P} ${INSTALL_PREFIX}/gkyl/bin/DataStruct
	find DataStruct -name '*.lua' | xargs cp --parents -fv -t ${INSTALL_PREFIX}/gkyl/bin
#
	${MKDIR_P} ${INSTALL_PREFIX}/gkyl/bin/Eq
	find Eq -name '*.lua' | xargs cp --parents -fv -t ${INSTALL_PREFIX}/gkyl/bin
#
	${MKDIR_P} ${INSTALL_PREFIX}/gkyl/bin/Updater
	find Updater -name '*.lua' | xargs cp --parents -fv -t ${INSTALL_PREFIX}/gkyl/bin
#
	${MKDIR_P} ${INSTALL_PREFIX}/gkyl/bin/App
	find App -name '*.lua' | xargs cp --parents -fv -t ${INSTALL_PREFIX}/gkyl/bin
#
	${MKDIR_P} ${INSTALL_PREFIX}/gkyl/bin/Cuda
	find Cuda -name '*.lua' | xargs cp --parents -fv -t ${INSTALL_PREFIX}/gkyl/bin
#
	${MKDIR_P} ${INSTALL_PREFIX}/gkyl/bin/Comm
	find Comm -name '*.lua' | xargs cp --parents -fv -t ${INSTALL_PREFIX}/gkyl/bin
#
	${MKDIR_P} ${INSTALL_PREFIX}/gkyl/bin/Io
	find Io -name '*.lua' | xargs cp --parents -fv -t ${INSTALL_PREFIX}/gkyl/bin
#
	${MKDIR_P} ${INSTALL_PREFIX}/gkyl/bin/Basis
	find Basis -name '*.lua' | xargs cp --parents -fv -t ${INSTALL_PREFIX}/gkyl/bin

clean:
	rm -rf gkyl gkyl.dSYM
	rm -rf ${BUILD_DIR}

# command to make dir
MKDIR_P ?= mkdir -p

# From: https://www.client9.com/self-documenting-makefiles/
.PHONY: help
help: ## Show help
	@echo "Makefile help. You can set parameters on the command line:"
	@echo ""
	@awk -F ':|##' '/^[^\t].+?:.*?##/ {\
        printf "\033[36m%-30s\033[0m %s\n", $$1, $$NF \
        }' $(MAKEFILE_LIST)
