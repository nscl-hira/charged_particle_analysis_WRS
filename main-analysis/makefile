COMPILER = g++
INCLUDE := `root-config --libs --cflags`
INCLUDE := ${INCLUDE} -I./include

.PHONY: all clean

all: anal

% : %.cpp src/hira.cpp src/reader.cpp
	${COMPILER} $^ -o $@.exe ${INCLUDE} 

clean:
	rm *.exe
