CXX = g++
CXXFLAGS = -std=c++11 -O3 -fopenmp -I./include -w
LDFLAGS = -L./lib -lz -lspoa
LIBDIR := -L.

SRCS = ./src/main.cpp ./src/lightWightMSA.cpp ./src/sequenceUtils.cpp ./src/edlib.cpp ./src/monomer.cpp ./src/sampleDBSCAN.cpp ./src/monoRefine.cpp ./src/hor.cpp ./src/genome.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = centroAnno

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET) $(LDFLAGS) 

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
