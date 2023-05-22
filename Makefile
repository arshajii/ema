TARGET = ema
LIBS = -L$(BWADIR) -lbwa -lm -lz -lpthread
CC = gcc
WARNINGS = -Wall -Wextra
CFLAGS = -std=gnu99 -march=x86-64 -O3 -fopenmp -flto -fstrict-aliasing $(WARNINGS)
LFLAGS = -lstdc++ -march=x86-64 -O3 -flto -fopenmp -lpthread

CXX = g++
CPPFLAGS = -c -std=c++11 -O3 -march=x86-64 -pthread
LDFLAGS = -pthread

#CFLAGS = -std=gnu99 -fstrict-aliasing -fopenmp -ggdb -O0 $(WARNINGS)
#LFLAGS = -fopenmp -lpthread

SRCDIR = src
OBJDIR = obj
INCDIR = include
CPPDIR = cpp
BWADIR = bwa

.PHONY: default all preproc bwa clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(wildcard $(SRCDIR)/*.c)) $(patsubst $(CPPDIR)/%.cc, $(CPPDIR)/%.o, $(wildcard $(CPPDIR)/*.cc))
HEADERS = $(wildcard $(INCDIR)/*.h $(BWADIR)/*.h)

preproc:
	$(MAKE) -C $(CPPDIR)

bwa:
	$(MAKE) -C $(BWADIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.c $(HEADERS)
	$(CC) $(CFLAGS) -I. -I$(INCDIR) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS) preproc bwa
	$(CXX) $(LDFLAGS) $(LFLAGS) $(OBJECTS) $(LIBS) -o $@

clean:
	-rm -f $(OBJDIR)/*.o
	$(MAKE) clean -C $(CPPDIR)
	$(MAKE) clean -C $(BWADIR)

