ifeq ($(OS),Windows_NT)
	PC = python
	RM = del /q
	VR = echo "nie ma valgrinda na windowsie"
else
	ifdef VIRTUAL_ENV
		PC = python
	else
		PC = python3
	endif
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
SCRIPT = plot.py
DATA = 
MEDIA = 

all: gotowe

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

gotowe: $(EXEC)
	./$(EXEC)

.PHONY: all clean safe gotowe
