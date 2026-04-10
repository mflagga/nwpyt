ifeq ($(OS),Windows_NT)
	PC = python
	RM = del /q
else
	PC = python3
	RM = rm -f
endif

.SUFFIXES:

CXX = gcc
LIBS = -lm
CFLAGS = -Wall -Wextra -Werror -Wpedantic -O3
SRCS = main.c
HDRS = header.h
EXEC = main.out
SCRIPT = wykres.py
DATA = Aa.csv Amisc.csv Ab.csv Ba.csv BEL.csv Bmisc.csv
MEDIA = Aa.png Ab.png Ba.png

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

.PHONY: all clean