target := plato-viewer

CXX := g++
PKGCFLAG := `/opt/local/bin/pkg-config --cflags gtkmm-2.4 libglademm-2.4 glibmm-2.4 gthread-2.0 gtkglextmm-1.2`
PKGLIB := `/opt/local/bin/pkg-config --libs gtkmm-2.4 glibmm-2.4 libglademm-2.4 gthread-2.0 gtkglextmm-1.2`
LIBS := -lmethods -lm -lyaml
METHODDIR := /Users/gilesjt/development/libraries/method_lib
OTHERDIR := /opt/local/include
LIBDIR := -L/Users/gilesjt/development/libraries/lib
INCDIR := -I. -I$(METHODDIR) -I$(OTHERDIR)
CXXFLAGS := -O3 -DMAC -g $(INCDIR) $(PKGCFLAG)

sources := $(wildcard *.cpp)
objects := $(sources:.cpp=.o)

.PHONY : all clean
all : $(target)
	@echo done.
	@echo

$(target) : $(objects)
	@echo linking...
	@$(CXX) $(LIBDIR) $^ $(LIBS) $(PKGLIB) -o $@

.cpp.o :
	@echo compiling $<
	@$(CXX) -c $(CXXFLAGS) $< -o $@

clean :
	rm -f $(target) $(objects)
