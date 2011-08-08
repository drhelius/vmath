TESTS = test3
CPPFLAGS = -I../src -g2 -Wall
LDFLAGS = -ldl

# if this flags doesn't work, run pkg-config cppunit --cflags or pkg-config cppunit --libs 
# respectively and past outpt to CPPFLAGS or LDFLAGS respectively
CPPFLAGS += -I/opt/local/include 
LDFLAGS += -L/opt/local/lib -lcppunit 


all: $(TESTS)

test2: test2.o vmath.o
	g++ $(LDFLAGS) -o test2 test2.o vmath.o

test3: test3.o vmath.o
	g++ $(LDFLAGS) -o test3 test3.o vmath.o

vmath.o: ../src/vmath.cpp ../src/vmath.h
	g++ $(CPPFLAGS) -c $< -o $@

%.o: %.cpp
	g++ $(CPPFLAGS) -c $< -o $@

clean:
	rm -f test2 test3 test2.o test3.o vmath.o
