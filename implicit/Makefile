CXX = g++
CXXFLAGS = -g -std=c++11 -Wall -O0

SRC_DIR = src
BUILD_DIR = build
BIN_DIR = bin

HDRS = $(SRC_DIR)/species.h $(SRC_DIR)/plasma.h
SRCS = $(SRC_DIR)/main.cpp $(SRC_DIR)/species.cpp $(SRC_DIR)/plasma.cpp 
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

LDLIBS = -lblas -llapack 

all: $(BIN_DIR)/plas

$(BIN_DIR)/plas: $(OBJS) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDLIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(HDRS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR) $(BIN_DIR):
	mkdir -p $@

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)