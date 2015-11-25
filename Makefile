CXXFLAGS += -O2 -std=c++11 -g -Wall -Wextra
LFLAGS =
LIBS =
SRCFILES := $(wildcard *.cpp) $(wildcard include/*.cpp)
OBJFILES := $(patsubst %.cpp,%.o,$(wildcard *.cpp)) $(patsubst %.cpp,%.o,$(wildcard include/*.cpp))
EXEC = run

all: $(EXEC)
.PHONY: all clean

# Automatically determine dependencies of source files
depend: .depend
.depend: $(SRCFILES)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^ > ./.depend;
include .depend

$(EXEC): $(OBJFILES)
	$(CXX) $(LFLAGS) -o $@ $^ $(LIBS)

clean:
	-rm -f $(OBJFILES) $(EXEC) .depend
