# A sample Makefile for building Google Test and using it in user
# tests.  Please tweak it to suit your environment and project.  You
# may want to move it to your project's root directory.
#
# SYNOPSIS:
#
#   make [all]  - makes everything.
#   make TARGET - makes the given target.
#   make clean  - removes all files generated by make.

# Please tweak the following variable definitions as needed by your
# project, except GTEST_HEADERS, which you can use in your own targets
# but shouldn't modify.

# Points to the root of Google Test, relative to where this file is.
# Remember to tweak this if you move this file.
GTEST_DIR = .

# Where to find user code.
USER_DIR = ..

BIN_DIR = ./bin

LIBDIR = ../../lib/
OS = MAC
RLIB = -lR
RLIBS = $(LIBDIR)libR.a #libR.dylib $(LIBDIR)libRblas.dylib
LIB= -lmethods
MDRPDTDIR = mdrpdt

ifeq ($(OS),UNIX)
	CXXFLAGS += -DUNIX
endif
ifeq ($(OS),WIN)
	CXXFLAGS += -DWIN
endif
ifeq ($(OS),MAC)
	CXXFLAGS += -DMAC
endif

# Flags passed to the preprocessor.
CPPFLAGS += -I$(GTEST_DIR) -I$(GTEST_DIR)/include -I..

# Flags passed to the C++ compiler.
CXXFLAGS += -g -Wall -Wextra

# All tests produced by this Makefile.  Remember to add new tests you
# created to the list.
TESTS = data_set_test helper_test marker_test

# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h

# House-keeping build targets.

all : $(TESTS)

clean :
	rm -f $(BIN_DIR)/* gtest.a gtest_main.a *.o

# Builds gtest.a and gtest_main.a.

# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

# Builds a sample test.  A test should link with either gtest.a or
# gtest_main.a, depending on whether it defines its own main()
# function.

data_set_test.o : $(GTEST_DIR)/data_set_test.cc $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -L$(LIBDIR) -c $(GTEST_DIR)/data_set_test.cc

data_set_test : data_set_test.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -L$(LIBDIR) $(LIB) $^ -o $(BIN_DIR)/$@

helper_test.o : $(GTEST_DIR)/helper_test.cc $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -L$(LIBDIR) -c $(GTEST_DIR)/helper_test.cc

helper_test : helper_test.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -L$(LIBDIR) $(LIB) $^ -o $(BIN_DIR)/$@

marker_test.o : $(GTEST_DIR)/marker_test.cc $(GTEST_HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -L$(LIBDIR) -c $(GTEST_DIR)/marker_test.cc

marker_test : marker_test.o gtest_main.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -L$(LIBDIR) $(LIB) $^ -o $(BIN_DIR)/$@
