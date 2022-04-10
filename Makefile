CC        = gcc
EXEC      = nanomd
CFLAGS    = -Wall -Wextra -Wunused
CDEBUG    = -g -ggdb -dH
COPTIMIZE = -Wuninitialized -O3
CLIBS     = -lm

CSRCS     = $(wildcard *.c)
CHDRS     = $(wildcard *.h)
TXTS      = $(wildcard *.txt README* LICENSE*)
SCRIPTS   = $(wildcard Makefile* *.sh)

NAME	:= "NanoMD"
VERSION	:= $(shell sed -n 's/.*VERSION \"\(.*\)\".*/\1/p' MDTest.c)
CPUARCH	:= $(shell uname -m)

ifeq ($(MAKECMDGOALS),debug)
	CFLAGS += $(CDEBUG)
	EXEC   := $(addsuffix -debug, $(EXEC))
else ifeq ($(MAKECMDGOALS),pack)
	EXEC   := $(addsuffix -v$(VERSION), $(EXEC))
else
	CFLAGS += $(COPTIMIZE)
endif

.PHONY: all clean pack

all: clean bin

debug: bin

bin:
	@echo :: Compiling \"$(NAME) v$(VERSION)\" \($(CPUARCH)\) ...
	$(CC) $(CFLAGS) $(CSRCS) -o $(EXEC) $(CLIBS)
	@echo :: Done

clean:
	@echo :: Cleaning up ...
	@rm -f $(EXEC) $(EXEC)-gdb $(EXEC)-debug $(EXEC)-v$(VERSION).tar.gz

pack:
	@echo :: Packing files ...
	tar -cvzhf $(EXEC).tar.gz $(CSRCS) $(CHDRS) $(TXTS) $(SCRIPTS)
	@echo :: Done
