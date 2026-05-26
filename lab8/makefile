ifeq ($(OS),Windows_NT)
	PC = python
	RM = del /q
	RMDIR = rmdir /s /q
	MV = move
	DEVNULL = NUL 2>&1
	VR = @echo nie ma valgrinda na windowsie
	SLASH = $(strip \)
	AVE = ..$(SLASH)nwenv$(SLASH)Scripts$(SLASH)activate
else
	ifdef VIRTUAL_ENV
		PC = python
	else
		PC = python3
	endif
	RM = rm -f
	RMDIR = rm -rf
	MV = mv
	SLASH = /
	DEVNULL = /dev/null
	VR = valgrind ./$(EXEC)
	AVE = source ../nwenv/bin/activate
endif

.SUFFIXES:

CXX = gcc
LIBS = -lm
CFLAGS = -Wall -Wextra -Werror -Wpedantic -O3
SRCS = main.c def.c
OBJS = main.o def.o
HDRS = dek.h
SCRIPT = plot.py
TEX = spr.tex
EXEC = main.out
REPORT = sprawozdanie.pdf
DATA = 
MEDIA = 

all: $(REPORT)

%.o: %.c $(HDRS)
	$(CXX) $< $(CFLAGS) -o $@ -c

$(EXEC): $(OBJS)
	$(CXX) $(OBJS) $(LIBS) -o $@

$(DATA): $(EXEC)
	./$<

$(MEDIA): $(SCRIPT) $(DATA)
	$(PC) $<

$(REPORT): $(TEX) $(MEDIA)
	latexmk -pdf -interaction=nonstopmode -output-directory=build $< > $(DEVNULL)
	$(MV) build$(SLASH)spr.pdf $@ > $(DEVNULL)

clean:
	$(RM) $(EXEC)
	$(RM) $(OBJS)
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