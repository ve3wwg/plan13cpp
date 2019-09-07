.PHONY:	all clean clobber install

all:	main

install: all

OBJS	= plan13.o main.o

LDFLAGS = -L$(PREFIX)/lib 

main: $(OBJS)
	$(CXX) $(OBJS) -o main

clean:	
	rm -f *.o *.x1o a.out core core.*

clobber: clean
	@rm -f .errs.t main

include Makefile.incl
PREFIX	= $(HOME)/local

plan13.o: plan13.hpp

# End Makefile
