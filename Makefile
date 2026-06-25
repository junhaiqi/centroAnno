CXX = g++
CXXFLAGS = -std=c++11 -O3 -fopenmp -I./include -w
LDFLAGS = -L./lib -lz -lspoa

SRCS = ./src/main.cpp ./src/lightWightMSA.cpp ./src/sequenceUtils.cpp ./src/edlib.cpp ./src/monomer.cpp ./src/sampleDBSCAN.cpp ./src/monoRefine.cpp ./src/hor.cpp ./src/genome.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = centroAnno

# spoa library (built from source for local compatibility)
SPOA_VERSION = 4.1.5
SPOA_SRCS = lib/spoa/src/alignment_engine.cpp \
            lib/spoa/src/graph.cpp \
            lib/spoa/src/simd_alignment_engine_dispatcher.cpp \
            lib/spoa/src/sisd_alignment_engine.cpp \
            lib/spoa/src/version.cpp
SPOA_OBJS = $(SPOA_SRCS:.cpp=.o)
SPOA_LIB = lib/libspoa.a
SPOA_INC = -I./lib/spoa/include -I./lib/spoa/src

all: $(TARGET)

$(TARGET): $(OBJS) $(SPOA_LIB)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET) $(LDFLAGS)

# Pattern rule for centroAnno sources
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Pattern rule for spoa sources (needs spoa headers)
lib/spoa/src/%.o: lib/spoa/src/%.cpp lib/spoa/src/spoa_config.h
	$(CXX) $(CXXFLAGS) $(SPOA_INC) -c $< -o $@

# Build static library from spoa objects
$(SPOA_LIB): $(SPOA_OBJS)
	ar rcs $@ $(SPOA_OBJS)

# Generate spoa_config.h from template
lib/spoa/src/spoa_config.h: lib/spoa/src/spoa_config.h.in
	sed 's/@SPOA_VERSION@/$(SPOA_VERSION)/' $< > $@

clean:
	rm -f $(OBJS) $(TARGET) $(SPOA_OBJS) $(SPOA_LIB) lib/spoa/src/spoa_config.h
