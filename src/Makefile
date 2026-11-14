CC = mkoctfile
FLAGS = -Wall -Wextra -Ofast 

srcs = $(wildcard *.cpp)
octs = $(srcs:.cpp=.oct)

all: $(octs)

%.oct: %.cpp
	$(CC) $(FLAGS) $< -o $@

clean:
	rm -f *.oct
	touch *.cpp
