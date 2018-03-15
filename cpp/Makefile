GIT_VERSION:=$(shell git describe --dirty --always --tags)
SOURCES:=$(wildcard ./*.cc)
OBJECTS:=$(SOURCES:.cc=.o)
DEP := $(OBJECTS:.o=.d)

CPPFLAGS += -MMD -MP -I.

.PHONY: all clean

EXECUTABLE=ema-preprocess
EXE=ema-preprocess

all: CPPFLAGS+=-g -O2
all: $(SOURCES) $(OBJECTS)

sanitize: CXX=g++
sanitize: CPPFLAGS+= -g -O1 -fno-omit-frame-pointer -fsanitize=address
sanitize: LDFLAGS+=-fno-omit-frame-pointer -fsanitize=address
sanitize: $(SOURCES) $(EXECUTABLE)

release: CPPFLAGS+=-g -O3 -DNDEBUG
release: $(SOURCES) $(EXECUTABLE)

debug: CPPFLAGS+=-g
debug: $(SOURCES) $(EXECUTABLE)

superdebug: CPPFLAGS+=-g -O0 -fno-inline
superdebug: $(SOURCES) $(EXECUTABLE)

profile: CPPFLAGS+=-g -pg -O2
profile: LDFLAGS+=-pg
profile: $(SOURCES) $(EXECUTABLE)

gprofile: CXX=g++
gprofile: LDFLAGS=-Wl,--no-as-needed,-lprofiler,--as-needed -ltcmalloc -lrt -lz -fopenmp
gprofile: CPPFLAGS+=-g -O2
gprofile: $(SOURCES) $(EXECUTABLE)

#$(EXECUTABLE): $(OBJECTS)
#	$(CXX) $(OBJECTS) $(LDFLAGS) -o $(EXE)

.cc.o:
	$(CXX) $(CPPFLAGS) -DGITVER=\"$(GIT_VERSION)\" $< -o $@

-include $(DEP)

clean:
	rm -rf $(OBJECTS) $(DEP)
