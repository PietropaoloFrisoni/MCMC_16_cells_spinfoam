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
TOOLSDIR = tools

# folders for external libraries
WIGDIR = $(EXTDIR)/wigxjpf
FASTWIGDIR = $(EXTDIR)/fastwigxj

# folders with header files for parallel hash tables
PARALLELHASHMAPDIR = $(EXTDIR)/parallel_hashmap

###############################################################################################

# call 'make DEBUG=1' to compile in debug mode 
ifeq ($(DEBUG), 1)
CXXFLAGS =  -std=c++17 -fopenmp -Wall -g -Og -fPIC  \
            -I$(WIGDIR)/inc -I$(FASTWIGDIR)/inc -Isrc -I$(INCDIR)/ -I$(PARALLELHASHMAPDIR) -DDEBUG_ON=1 # debug
else
CXXFLAGS =  -std=c++17 -fopenmp -O3 -fPIC -march=native -fno-math-errno  \
            -I$(WIGDIR)/inc -I$(FASTWIGDIR)/inc -Isrc -I$(INCDIR)/ -I$(PARALLELHASHMAPDIR) -DDEBUG_OFF=1 # optimized
endif

LDFLAGS = $(CXXFLAGS) -L$(WIGDIR)/lib/ -L$(FASTWIGDIR)/lib/ -Wl,--no-as-needed -Wl,--allow-shlib-undefined

LDLIBS = -ldl -lgsl -lgslcblas -lwigxjpf -lfastwigxj -lm -lwigxjpf_quadmath -lquadmath

###############################################################################################

INCS = inc/mcmc.h inc/error.h inc/common.h inc/hash_21j_symbols.h

_OBJS =  mcmc.o  hash_21j_symbols.o python_mirror.o

_TOOLS = test_jsymbols Metropolis_Hastings_parallel_run Hashing_21j 

OBJS = $(patsubst %,$(OBJDIR)/%,$(_OBJS))
TESTS = $(patsubst %,$(BINDIR)/%,$(_TESTS))
TOOLS = $(patsubst %,$(BINDIR)/%,$(_TOOLS))


# library/src object files 
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(INCS)
	@echo "   CXX    $@"
	$(QUIET)mkdir -p $(dir $@)
	$(QUIET)$(CXX) $(CXXFLAGS) -c -o $@ $< 


# shared library
$(LIBDIR)/libspinfoam.so: $(OBJS)
	@echo "   CXX    $@"
	$(QUIET)mkdir -p $(dir $@)
	$(QUIET)$(CXX) -shared $(OBJS) -o $@ $(LDFLAGS) $(LDLIBS) 

	
# compile tools programs
$(BINDIR)/%: $(TOOLSDIR)/%.cpp $(OBJS)
	@echo "   CXX    $@"
	$(QUIET)mkdir -p $(dir $@)
	$(QUIET)$(CXX) $(CXXFLAGS) -o $@ -Iinc/ $< -Llib/ $(LDFLAGS) -lspinfoam $(LDLIBS)	

###############################################################################################

.DEFAULT_GOAL := default

default: all

all: lib generate_executable 

lib: $(LIBDIR)/libspinfoam.so

generate_executable: lib $(TOOLS)

clean:
	rm -rf $(OBJDIR)
	rm -rf $(LIBDIR)
	rm -rf $(BINDIR)

