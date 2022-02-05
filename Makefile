QUIET ?= @

CXX = g++ 

# folder names 
EXTDIR = ext
SRCDIR = src
OBJDIR = obj
INCDIR = inc
LIBDIR = lib
BINDIR = bin
TESTDIR = test


# folders for external libraries
WIGDIR = $(EXTDIR)/wigxjpf
FASTWIGDIR = $(EXTDIR)/fastwigxj

###############################################################################################

# call 'make DEBUG=1' to compile in debug mode
ifeq ($(DEBUG), 1)
CXXFLAGS =  -std=c++11 -g -O0  -fopenmp -Wall -Og -fPIC -I$(WIGDIR)/inc -I$(FASTWIGDIR)/inc -Isrc -I$(INCDIR)/  # debug
else
CXXFLAGS =  -std=c++11 -fopenmp -Wall -fPIC -I$(WIGDIR)/inc -I$(FASTWIGDIR)/inc -Isrc -I$(INCDIR)/ # optimized
endif

LDFLAGS = $(CXXFLAGS) -L$(WIGDIR)/lib/ -L$(FASTWIGDIR)/lib/ -Wl,--no-as-needed -Wl,--allow-shlib-undefined

LDLIBS = -ldl -lgsl -lgslcblas -lwigxjpf -lfastwigxj -lm -lwigxjpf_quadmath -lquadmath

###############################################################################################

INCS = inc/jsymbols.h inc/mcmc.h inc/ampl.h inc/utilities.h inc/error.h inc/common.h inc/config.h inc/setup.h

_OBJS = libshared.o jsymbols.o mcmc.o first_file.o ampl.o setup.o

_TESTS = test_jsymbols test_sampler

OBJS = $(patsubst %,$(OBJDIR)/%,$(_OBJS))
TESTS = $(patsubst %,$(BINDIR)/%,$(_TESTS))
TOOLS = $(patsubst %,$(BINDIR)/%,$(_TOOLS))


# library/src object files     REMOVED $(INCS) AT END LINE
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp 
	@echo "   CXX    $@"
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c -o $@ $< 

# shared library
$(LIBDIR)/libshared.so: $(OBJS)
	@echo "   CXX    $@"
	mkdir -p $(dir $@)
	$(CXX) -shared $(OBJS) -o $@ $(LDFLAGS) $(LDLIBS) 
	
# compile test programs
$(BINDIR)/%: $(TESTDIR)/%.cpp $(OBJS)
	@echo "   CXX    $@"
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -o $@ -Iinc/ $< -Llib/ $(LDFLAGS) -lshared $(LDLIBS)	

###############################################################################################

.DEFAULT_GOAL := default

default: all

all: lib tests 

lib: $(LIBDIR)/libshared.so

tests: lib $(TESTS)

clean:
	rm -rf $(OBJDIR)
	rm -rf $(LIBDIR)
	rm -rf $(BINDIR)

