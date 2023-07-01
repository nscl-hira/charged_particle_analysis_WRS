COMPILER := g++
INCLUDE := `root-config --libs --cflags`
FLAGS := -O2 -fPIC

SRC_DIR := ${PROJECT_DIR}/src
SUB_DIR := ${SRC_DIR}/e15190 ${SRC_DIR}/e15190/HiRA

VPATH := ${SRC_DIR} ${SUB_DIR}
INCLUDE += -I${SRC_DIR} -I${SRC_DIR}/e15190 -I${SRC_DIR}/e15190/HiRA

SRC := ${shell find ${SRC_DIR} -name "*.cpp"}

.PHONY: all clean

all: anal

anal : anal.cpp anal.hh HistogramManager.cpp HistogramManager.hh ${SRC} 
	${COMPILER} ${PROJECT_DIR}/analysis/$^ -o $@.exe ${INCLUDE}

clean:
	rm -f *.exe
