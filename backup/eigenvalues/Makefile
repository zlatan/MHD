CC=gcc
CFLAGS=-c -Wall -Ic:/gsl
LDFLAGS= -Lc:/gsl
LIBS= -lgsl
SOURCES=eigenvalues.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=hello

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
    $(CC) $(LDFLAGS) $(OBJECTS) $(LIBS) -o $@

.c.o:
    $(CC) $(CFLAGS) $< -o $@

