OBJS=fft.o
CXXFLAGS=-g -std=c++17 -Wall --pedantic
LDFLAGS=-lm
CXX=clang++
all: fft

fft: $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	$(RM) $(OBJS) fft
