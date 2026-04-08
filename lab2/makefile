.SUFFIXES:

CXX = gcc
LIBS = -lm
CFLAGS = -Wall -Wextra -Werror -Wpedantic -O3
SRCS = main.c
HDRS = header.h
EXEC = main.out
PC = python3
SCRIPT = wykres.py
DATA = Aa.csv Amisc.csv
MEDIA = Aa.png

all: $(MEDIA)

$(EXEC): $(SRCS) $(HDRS)
	$(CXX) $(SRCS) $(CFLAGS) $(LIBS) -o $@

$(DATA): $(EXEC)
	./$<

$(MEDIA): $(SCRIPT) $(DATA)
	$(PC) $<

clean:
	rm -f $(EXEC)
	rm -f $(DATA)
	rm -f $(MEDIA)

.PHONY: all clean