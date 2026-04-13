ifeq ($(OS),Windows_NT)
	PC = python
	RM = del /q
	VR = echo "nie ma valgrinda na windowsie"
else
	PC = python3
	RM = rm -f
	VR = valgrind ./$(EXEC)
endif

.SUFFIXES:

CXX = gcc
LIBS = -lm
CFLAGS = -Wall -Wextra -Werror -Wpedantic -O3
SRCS = main.c
HDRS = header.h
EXEC = main.out
SCRIPT = wykres.py
DATA = Aa.csv Amisc.csv Ab.csv Ba.csv BEL.csv Bmisc.csv Ca.csv Cmisc.csv
MEDIA = Aa.png Ab.png Ba.png Ca.png

all: $(MEDIA)

$(EXEC): $(SRCS) $(HDRS)
	$(CXX) $(SRCS) $(CFLAGS) $(LIBS) -o $@

$(DATA): $(EXEC)
	./$<

$(MEDIA): $(SCRIPT) $(DATA)
	$(PC) $<

clean:
	$(RM) $(EXEC)
	$(RM) $(DATA)
	$(RM) $(MEDIA)

safe: $(EXEC)
	$(VR)

.PHONY: all clean safe
