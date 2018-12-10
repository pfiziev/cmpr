CC = gcc
OS := $(shell uname -s)
FLAGS = -Wall -O3 -std=gnu99 -g

ifeq ($(OS),Linux)
	LDFLAGS_MULTIPLE_EC = -Wl,-soname,csdelta.so
endif
ifeq ($(OS),Darwin)
	LDFLAGS_MULTIPLE_EC = -Wl,-install_name,csdelta.so
endif

all: csdelta.so

csdelta.so : csdelta.o
	$(CC) -shared $(LDFLAGS_MULTIPLE_EC) $(FLAGS) -o csdelta.so csdelta.o

csdelta.o : csdelta.c
	$(CC) -c -fPIC $(FLAGS) csdelta.c -o csdelta.o


clean:
	rm -f csdelta.so csdelta.o

