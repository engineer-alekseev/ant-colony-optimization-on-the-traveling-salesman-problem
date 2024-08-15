#all, clean, ant, gauss, winograd
CC = g++ -fopenmp
CC_FLAGS = -Wall -Wextra -Werror -std=c++17 #-pedantic-errors
MFLAGS=-MMD -MP
SANITIZER_FLAGS = -fsanitize=address
SRCS=$(wildcard *.cc,*/*.cc)
OBJ_FILES=$(patsubst %.cc,%.o,$(SRCS))
REPORT_DIR = report


all: ant

ant: clean
	$(CC) $(CC_FLAGS) interface/interface.cc -o app
	./app ant


clean:
	rm -rf *.o *.a run *.dSYM *.d app *Result.txt *Generated.txt

.PHONY:
	all, clean, ant

clang:
	clang-format -i --style=Google **/*.cc **/*.h
	clang-format -n --style=Google **/*.cc **/*.h

cppcheck:
	cppcheck --enable=all --language=c++ --std=c++17 --suppress=missingIncludeSystem console/console.h aco/* matrix/*
