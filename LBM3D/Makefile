#!/bin/bash

FOLDER=.
KOKKOS_PATH= /public/home/wangwu/Kokkos/kokkos-master
KOKKOS_DEVICES = HIP
#KOKKOS_USE_TPLS = hwloc
KOKKOS_OPTIONS="disable_deprecated_code"
#HWLOC_PATH=/usr/include/hwloc

SRC_PATH = $(FOLDER)/source
LBM_PATH = $(SRC_PATH)/include
RELEASE = $(FOLDER)/release
EXE_NAME = $(FOLDER)/release/lbm

SRC = $(wildcard $(SRC_PATH)/*.cpp)
HDR = $(wildcard *.h)

default: build
#	echo "Start Build"

CXX = hipcc
EXE = ${EXE_NAME}.hip
KOKKOS_ARCH = "ZEN,VEGA906"
CXXFLAGS = -Werror -O3 -march=native -mtune=native -DNDEBUG                

INCLUDE = -I$(LBM_PATH)

override CXXFLAGS += -I./ 
LINK = ${CXX}
LINKFLAGS =

DEPFLAGS = -M

vpath %.cpp $(sort $(dir $(SRC)))

OBJ = $(addprefix $(RELEASE)/,$(notdir $(SRC:.cpp=.o)))
LIB =

include $(KOKKOS_PATH)/Makefile.kokkos

$(warning $(OBJ) $(EXE) $(sort $(dir $(SRC))))

build: $(EXE)

$(EXE): $(OBJ) $(KOKKOS_LINK_DEPENDS) 
	$(LINK) $(KOKKOS_LDFLAGS) $(LINKFLAGS) $(EXTRA_PATH) $(OBJ) $(KOKKOS_LIBS) $(LIB) -o $(EXE)

clean: kokkos-clean
	rm -f $(RELEASE)/*.o $(RELEASE)/*.cuda $(RELEASE)/*.host $(RELEASE)/*.a $(RELEASE)/*.h $(RELEASE)/*.tmp $(RELEASE)/*.x $(EXE)

# Compilation rules
$(RELEASE)/%.o:%.cpp $(KOKKOS_CPP_DEPENDS) $(HDR)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CXXFLAGS) $(INCLUDE) -c $< -o $@
