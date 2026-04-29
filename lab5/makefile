ifeq ($(OS),Windows_NT)
	PC = python
	RM = del /q
	VR = @echo nie ma valgrinda na windowsie
	AVE = ../nwenv/Scripts/activate
else
	ifdef VIRTUAL_ENV
		PC = python
	else
		PC = python3
	endif
	RM = rm -f
	VR = valgrind ./$(EXEC)
	AVE = source ../nwenv/bin/activate
endif

.SUFFIXES:

CXX = gcc
LIBS = -lm
CFLAGS = -Wall -Wextra -Werror -Wpedantic -O3
SRCS = main.c
HDRS = header.h
EXEC = main.out
SCRIPT = plot.py
DATA = psifile.csv
MEDIA = psi.png

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

safe: $(EXEC)
	$(VR)

ve: 
	@echo Run: $(AVE)

gotowe: $(EXEC)
	./$(EXEC)

.PHONY: all clean safe gotowe ve
