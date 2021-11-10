CC=gcc
CFLAGS= -Wall -Werror -pedantic -fPIC -O3
LIBS= -lgsl -lgslcblas -lm
LDFLAGS= -shared
RM= rm -f

TARGET= cloglik.so
SOURCES= $(wildcard *.c)
HEADERS= $(wildcard *.h)

.PHONY: all
all: ${TARGET}

$(TARGET): $(SOURCES)
	$(CC) $(CFLAGS) ${LDFLAGS} -o $@  $^  $(LIBS) 

.PHONY: clean
clean:
	-${RM} ${TARGET}
