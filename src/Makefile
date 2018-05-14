OBJS=iterative_fft.o
CXXFLAGS=-O2 -std=c++17 -Wall --pedantic
LDFLAGS=-lm
CXX=clang++
all: fft

fft: $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	$(RM) $(OBJS) fft
