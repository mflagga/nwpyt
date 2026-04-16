ifeq ($(OS),Windows_NT)
	PC = python
	RM = del /q
	CFLAGS = -arch=native -O3
	LIBS = 
else
	PC = python3
	RM = rm -f
	CFLAGS = -arch=native -Xcompiler -Wall,-Wextra -Werror all-warnings -O3
	LIBS = -lm
endif

.SUFFIXES:

CXX = nvcc
EXEC = main.out
SRCS = main.cu
HDRS = header.cuh
SCRIPT = plot.py
DATA = S1a.csv
MEDIA = img/S1mask.png

all: $(MEDIA)

$(EXEC): $(SRCS) $(HDRS)
	$(CXX) $(SRCS) $(CFLAGS) $(LIBS) -o $@

$(DATA): $(EXEC)
	./$<

$(MEDIA): $(SCRIPT) $(DATA)
	$(PC) $(SCRIPT)

clean:
	$(RM) $(EXEC)
	$(RM) $(DATA)
	$(RM) $(MEDIA)

.PHONY: all clean
