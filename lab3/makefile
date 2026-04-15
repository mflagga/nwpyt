ifeq ($(OS),Windows_NT)
	PC = python
	RM = del /q
else
	PC = python3
	RM = rm -f
endif

.SUFFIXES:

CXX = nvcc
EXEC = main.out
CFLAGS = -arch=native -Xcompiler -Wall,-Wextra -Werror all-warnings -O3
LIBS = -lm
SRCS = main.cu
HDRS = header.cuh
SCRIPT = plot.py

all:

$(EXEC): $(SRCS) $(HDRS)
	$(CXX) $(SRCS) $(CFLAGS) $(LIBS) -o $@

clean:
	$(RM) $(EXEC)

.PHONY: all clean
