OBJS=fft.o
CXXFLAGS=-O2 -Wall --pedantic
LDFLAGS=-lm

all: fft

fft: $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

clean:
	$(RM) $(OBJS) fft
