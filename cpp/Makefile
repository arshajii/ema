GIT_VERSION:=$(shell git describe --dirty --always --tags)
SOURCES:=$(wildcard ./*.cc)
OBJECTS:=$(SOURCES:.cc=.o)
DEP := $(OBJECTS:.o=.d)

CPPFLAGS = -c -MMD -MP -I. -std=c++11

.PHONY: all clean

all: CPPFLAGS+=-g -O3
all: $(SOURCES) $(OBJECTS)

.cc.o:
	$(CXX) $(CPPFLAGS) -DGITVER=\"$(GIT_VERSION)\" $< -o $@

-include $(DEP)

clean:
	rm -rf $(OBJECTS) $(DEP)
