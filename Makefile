CPPFLAGS := -g -I. $(shell root-config --cflags)
CXXFLAGS := -Wall -Wextra -ggdb3 -O9 -march=native -mfpmath=sse -msse3 -ftree-vectorize -pipe -DNDEBUG
LDFLAGS := $(CXXFLAGS) `root-config --glibs`

SUFFIXES := .o .cc .cpp

EXES := $(basename $(wildcard *.cpp))
SRCS := $(wildcard *.cc) $(wildcard *.cpp)
OBJS := $(patsubst %.cc, %.o, $(patsubst %.cpp, %.o, $(SRCS)))
DEPS := $(OBJS:.o=.P)

define cxx_compile_with_dependency_creation
	$(COMPILE.cc) -MD -o $@ $<
	@sed -e 's|.*:|$*.o:|' <$*.d >$*.P
	@sed -e 's/.*://' -e 's/\\$$//' <$*.d | fmt -1 | \
	  sed -e 's/^ *//' -e 's/$$/:/' >>$*.P
	@rm -f $*.d
endef

define cxx_link_rule
	$(CXX) $^ $(LDFLAGS) $(LOADLIBES) $(LDLIBS) -o $@
endef

%.o: %.cc
	$(call cxx_compile_with_dependency_creation)
        
%.o: %.cpp
	$(call cxx_compile_with_dependency_creation)

%: %.o
	$(call cxx_link_rule)

%.Dict.cc %.Dict.h: %.h %.LinkDef.h
	$(ROOTSYS)/bin/rootcint -f $@ -c -I. $^



#all: $(EXES)

wzAnalysis: wzAnalysis.o WZEvent.o EventTree.o Leptons.o

wzSelectionYields: wzSelectionYields.o WZEvent.o EventTree.o Leptons.o


clean:
	- $(RM) *.o *.Dict.cc *.Dict.h $(addsuffix .o, $(EXES)) Dependencies.make $(EXES) *.P

-include $(DEPS)
