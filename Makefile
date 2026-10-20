CXX ?= g++

TARGET := build/optymalizacja
SOURCES := main.cpp matrix.cpp ode_solver.cpp opt_alg.cpp solution.cpp user_funs.cpp
OBJECTS := $(SOURCES:%.cpp=build/%.o)
DEPS := $(OBJECTS:.o=.d)

CXXFLAGS ?= -std=c++17 -Wall -Wextra -pedantic
CPPFLAGS ?= -MMD -MP
LDFLAGS ?=
LDLIBS ?=

.PHONY: all run debug clean help

all: $(TARGET)

$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

build/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

debug: CXXFLAGS += -O0 -g
debug: clean all

clean:
	rm -rf build

help:
	@printf '%s\n' \
		'make        Build the project' \
		'make run    Build and run the program' \
		'make debug  Clean and build with debug symbols' \
		'make clean  Remove build artifacts'

-include $(DEPS)
