# Build the Catch2 unit-test binary. Run `make` in bin/ first so the core
# object files exist. This Makefile compiles:
#   - Catch2 amalgamated (vendored)
#   - Unit test .cpp files in this directory
# Then links against bin/*.o to reuse the compiled core.

CXX      := g++
# Catch2 v3 requires C++17; the core VikNGS code compiles fine at that level
# since the .o files are already built (we just link them).
CXXFLAGS := -std=c++17 -O2 -D_GNU_SOURCE
INCLUDES := -I. -I../../src -I../../src/Parser

BUILD    := build
BIN_OBJS := $(wildcard ../../bin/*.o)
# Exclude the CLI main object files — they pull in their own main()
BIN_OBJS := $(filter-out ../../bin/vikNGScmd.o ../../bin/simcmd.o,$(BIN_OBJS))

TEST_SRCS := $(wildcard test_*.cpp)
TEST_OBJS := $(TEST_SRCS:%.cpp=$(BUILD)/%.o)
CATCH_OBJ := $(BUILD)/catch_amalgamated.o

TARGET := ../../bin/test-unit

$(TARGET): $(CATCH_OBJ) $(TEST_OBJS) $(BIN_OBJS) | $(BUILD)
	$(CXX) $(CXXFLAGS) $(CATCH_OBJ) $(TEST_OBJS) $(BIN_OBJS) -pthread -o $@

$(CATCH_OBJ): catch_amalgamated.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) -c catch_amalgamated.cpp -o $@

$(BUILD)/%.o: %.cpp | $(BUILD)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(BUILD):
	mkdir -p $(BUILD)

clean:
	rm -rf $(BUILD) $(TARGET)

.PHONY: clean
