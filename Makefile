# 3D-FDTD: make, make tests, make check-tests, make regression, make check-layers
PROJECT_NAME = 3D-FDTD
PROJECT_BRIEF = "3D parallel FDTD solver"
OUTPUT_DIRECTORY = docs/output
DOXYGEN_INPUT = README.md
DOXYGEN_MAINPAGE = README.md
DOXYGEN_IMAGE_PATH = docs/dev/Images docs/geometry/shapes

CXX = mpicxx
TARGETS = fdtd
SRC = main.cc $(patsubst src/%,%,$(wildcard src/*/*.cc src/*/*/*.cc))

TEST_TARGETS = test_config test_machine test_carray test_h5io test_geometry
TEST_SRC = $(TEST_TARGETS:%=%.cc)
# the tests link the modules without solver state
TEST_DEPS = $(addprefix src/,config/ini.o config/sim_config.o config/chi_tensors.o io/h5io.o geometry/csg.o geometry/anapole.o)

# parallel HDF5 through pkg-config; another name: make HDF5_PKG=hdf5-openmpi
HDF5_PKG ?= hdf5
CPPFLAGS += $(shell pkg-config --cflags $(HDF5_PKG))
LDLIBS += $(shell pkg-config --libs $(HDF5_PKG))

CPPFLAGS += -DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX
CPPFLAGS += -DFDTD_VERSION=\"$(shell git describe --always --dirty 2>/dev/null || echo dev)\"

include rules.mk

CXXFLAGS += -Wall -Wextra

# optional DAMM kernels: make DAMM=/path/to/damm, see docs/dev/Roadmap.md
ifdef DAMM
CPPFLAGS += -DFDTD_DAMM -I$(DAMM)/inc
CXXFLAGS += -fopenmp
LDFLAGS += -fopenmp
endif

MPIEXEC ?= mpirun -np 4
check-tests: tests
	set -e; for t in $(TEST_TARGETS); do ./$$t; done
	$(MPIEXEC) ./test_h5io

regression: $(TARGETS)
	cd test/regression && python3 run_tests.py --program ../../fdtd $(REGRESSION_ARGS)

check-layers:
	python3 tools/check_layers.py --strict core io grid fields materials boundary update geometry sources monitors

.PHONY: check-tests regression check-layers
