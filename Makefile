CXX = c++
CXXFLAGS = -Wall -O3 -std=c++17
sources = $(wildcard *.cpp)

$(info $(sources))

all: main
main: $(sources)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -rf *.o main
