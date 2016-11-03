# Makefile for calSOV.cpp
CC         = g++
CFLAGS     = -Wall -O3
CFLAGS_DEBUG = -Wall -g3 -DDEBUG
LIBS       = -lm
LIB_PATH   = lib
INCLUDE    = ./
SRC        = calSOV.cpp
OBJ        = $(SRC:.cpp=.o) 
OBJ_DEBUG  = $(SRC:.cpp=_debug.o) 
EXE        = calSOV
EXE_DEBUG  = calSOV_debug
RM         = /bin/rm -f
CP         = /bin/cp -f
BINPATH    = /usr/local/bin

all: $(OBJ)
	$(CC) $(CFLAGS)  $(OBJ) -L$(LIB_PATH) -o $(EXE)  $(LIBS)

# -c flag tells the compiler to create only OBJ files
$(OBJ): $(SRC)
	$(CC) $(CFLAGS) -c -I$(INCLUDE) $(SRC) 




debug: $(OBJ_DEBUG)
	$(CC) $(CFLAGS_DEBUG) $(OBJ_DEBUG) -o $(EXE_DEBUG) $(LIBS)

$(OBJ_DEBUG): $(SRC)
	$(CC) $(CFLAGS_DEBUG) -c -I$(INCLUDE) $<  -o $@

install:
	$(CP)  $(EXE)  $(BINPATH)

clean:
	$(RM)  $(EXE) $(OBJ)
	$(RM)  $(EXE_DEBUG) $(OBJ_DEBUG)
