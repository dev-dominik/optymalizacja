CXX ?= c++

TARGET := optymalizacja
BUILD_DIR := build

CPPFLAGS ?= -MMD -MP
CXXFLAGS ?= -std=c++17 -Wall -Wextra -pedantic

SOURCES := main.cpp matrix.cpp ode_solver.cpp opt_alg.cpp solution.cpp user_funs.cpp
OBJECTS := $(SOURCES:%.cpp=$(BUILD_DIR)/%.o)
DEPENDS := $(OBJECTS:.o=.d)

.PHONY: all build run clean rebuild

all: build

build: $(BUILD_DIR)/$(TARGET)

run: build
	./$(BUILD_DIR)/$(TARGET) $(ARGS)

$(BUILD_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR)

rebuild: clean build

-include $(DEPENDS)
