TARGET = ema
LIBS = -L$(BWADIR)/bwa -lbwa -lm -lz -lpthread
CC = gcc
WARNINGS = -Wall -Wextra -Werror
CFLAGS = -std=gnu99 -march=native -O3 -flto -fstrict-aliasing $(WARNINGS)
LFLAGS = -march=native -O3 -flto
#CFLAGS = -std=gnu99 -fstrict-aliasing -ggdb -O0 $(WARNINGS)

SRCDIR = src
OBJDIR = obj
INCDIR = include
BWADIR = .

.PHONY: default all bwa clean

default: $(TARGET)
all: default

OBJECTS = $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(wildcard $(SRCDIR)/*.c))
HEADERS = $(wildcard $(INCDIR)/*.h $(BWADIR)/bwa/*.h)

bwa:
	$(MAKE) -C $(BWADIR)/bwa

$(OBJDIR)/%.o: $(SRCDIR)/%.c $(HEADERS)
	$(CC) $(CFLAGS) -I$(INCDIR) -I$(BWADIR) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS) bwa
	$(CC) $(LFLAGS) $(OBJECTS) $(LIBS) -o $@

clean:
	-rm -f $(OBJDIR)/*.o
	$(MAKE) clean -C $(BWADIR)/bwa

