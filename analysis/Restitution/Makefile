CXXFLAGS += -O2 -std=c++11 -g
OUTPUT_OPTION = -MMD -MP -o $@
LFLAGS =
LIBS =

EXE = Restitution_Portrait
EXE_SRC = exe/model_sim.cpp

.PHONY: all clean
all: $(EXE)

MODEL_SRC = $(wildcard ../../src/*.cpp)

# Dependencies are automatically generated
SRC := $(MODEL_SRC) $(EXE_SRC)
OBJ := $(SRC:.cpp=.o)
DEP := $(SRC:.cpp=.d)
-include $(DEP)

$(EXE): $(OBJ)
	$(CXX) $(CXXFLAGS) $(LFLAGS) -o $@ $^ $(LIBS)

clean:
	-rm -f $(OBJ) $(EXE) $(DEP)
