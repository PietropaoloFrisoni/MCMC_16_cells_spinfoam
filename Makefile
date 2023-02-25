QUIET ?= @

CXX = g++ 

# folder names 
EXTDIR = ext
SRCDIR = src
OBJDIR = obj
INCDIR = inc
LIBDIR = lib
BINDIR = bin
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

INCS = $(INCDIR)/chain_class.h $(INCDIR)/mcmc.h $(INCDIR)/error.h $(INCDIR)/common.h $(INCDIR)/hash_21j_symbols.h

_OBJS = chain_class.o mcmc.o hash_21j_symbols.o python_mirror.o

_TOOLS = Metropolis_Hastings_parallel_run Hashing_21j 

OBJS = $(patsubst %,$(OBJDIR)/%,$(_OBJS))
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

all: precompile_wigxjpf_fastwigxj lib generate_executable 

precompile_wigxjpf_fastwigxj: 
	@echo "Compiling wigxjpf..."
	$(QUIET) cd $(EXTDIR)/wigxjpf && make
	@echo "Compiling fastwigxj..."
	$(QUIET) cd $(EXTDIR)/fastwigxj && make

lib: $(LIBDIR)/libspinfoam.so

generate_executable: lib $(TOOLS)

clean:
	rm -rf $(OBJDIR)
	rm -rf $(LIBDIR)
	rm -rf $(BINDIR)
	@echo "cleaning wigxjpf..."
	$(QUIET) cd $(EXTDIR)/wigxjpf && make clean
	@echo "cleaning fastwigxj..."
	$(QUIET) cd $(EXTDIR)/fastwigxj && make clean

