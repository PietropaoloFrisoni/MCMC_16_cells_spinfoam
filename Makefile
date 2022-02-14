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

# call 'make DEBUG=1' to compile in debug mode (FAKE)
ifeq ($(DEBUG), 1)
CXXFLAGS =  -std=c++11 -g -O0  -fopenmp -Wall -Og -fPIC -I$(WIGDIR)/inc -I$(FASTWIGDIR)/inc -Isrc -I$(INCDIR)/ -Iparallel_hashmap # debug
else
CXXFLAGS =  -std=c++11 -g -fopenmp -Wall -fPIC -I$(WIGDIR)/inc -I$(FASTWIGDIR)/inc -Isrc -I$(INCDIR)/ -I$(WIGDIR)/inc -Iparallel_hashmap # optimized
endif

LDFLAGS = $(CXXFLAGS) -L$(WIGDIR)/lib/ -L$(FASTWIGDIR)/lib/ -Wl,--no-as-needed -Wl,--allow-shlib-undefined

LDLIBS = -ldl -lgsl -lgslcblas -lwigxjpf -lfastwigxj -lm -lwigxjpf_quadmath -lquadmath

###############################################################################################

INCS = inc/mcmc.h inc/error.h inc/common.h inc/hash_21j_symbols.h

_OBJS =  mcmc.o  hash_21j_symbols.o

_TESTS = test_jsymbols test_sampler Hashing_21j

OBJS = $(patsubst %,$(OBJDIR)/%,$(_OBJS))
TESTS = $(patsubst %,$(BINDIR)/%,$(_TESTS))
TOOLS = $(patsubst %,$(BINDIR)/%,$(_TOOLS))


# library/src object files 
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(INCS)
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

