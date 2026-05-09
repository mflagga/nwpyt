ifeq ($(OS),Windows_NT)
	PC = python
	RM = del /q
	RMDIR = rmdir /s /q
	MV = move
	DEVNULL = NUL
	VR = @echo nie ma valgrinda na windowsie
	AVE = ../nwenv/Scripts/activate
else
	ifdef VIRTUAL_ENV
		PC = python
	else
		PC = python3
	endif
	RM = rm -f
	RMDIR = rm -rf
	MV = mv
	DEVNULL = /dev/null
	VR = valgrind ./$(EXEC)
	AVE = source ../nwenv/bin/activate
endif

.SUFFIXES:

CXX = gcc
LIBS = -lm -lgsl -lgslcblas
CFLAGS = -Wall -Wextra -Werror -Wpedantic -O3
SRCS = main.c
HDRS = header.h
EXEC = main.out
SCRIPT = plot.py
DATA = Vmap.csv misc.csv Vnmap.csv Gfile.csv
MEDIA = Vmap.png Vnmap.png Gfile.png
TEX = sprawozdanie.tex
REPORT = sprawozdanie.pdf

all: $(REPORT)

$(EXEC): $(SRCS) $(HDRS)
	$(CXX) $(SRCS) $(CFLAGS) $(LIBS) -o $@

$(DATA): $(EXEC)
	./$<

$(MEDIA): $(SCRIPT) $(DATA)
	$(PC) $(SCRIPT)

$(REPORT): $(TEX) $(MEDIA)
	latexmk -pdf -silent -output-directory=build $< > $(DEVNULL)
	$(MV) ./build/$@ ./$@

clean:
	$(RM) $(EXEC)
	$(RM) $(DATA)
	$(RM) $(MEDIA)
	$(RMDIR) build

safe: $(EXEC)
	$(VR)

ve:
	@echo $(AVE)

gotowe: $(EXEC)
	./$(EXEC)

.PHONY: all clean safe gotowe ve
