SHELL=/bin/bash

CC=g++-4.7
CPPFLAGS += -std=c++11 -I include -W -Wall -mavx
CPPFLAGS += -O3

LDFLAGS =
LDLIBS = -lm -ltbb

TBB_DIR = /Users/richardworrall/tbb

TBB_INC_DIR = $(TBB_DIR)/include

TBB_LIB_DIR = $(TBB_DIR)/lib

CPPFLAGS += -I $(TBB_INC_DIR)
LDFLAGS += -L $(TBB_LIB_DIR)

# For your makefile, add TBB and OpenCL as appropriate
all: src/bitecoin_client