# programs
CC				:= gcc
CXX				:= g++
MV				:= mv -f
RM				:= rm -f

# flags
INCLUDE_CFLAGS	= -I./ -I/target/include -I./include
CFLAGS			= 	$(INCLUDE_CFLAGS) -g -O2 -Wall -fPIC
LDFLAGS			= 	-L/target/lib
LDFLAGS			+= 	-luhd -lliquid -lm -lc -lboost_system-mt -lboost_thread-mt
LDFLAGS			+= 	-lboost_program_options-mt -lpthread -lvolk -lfftw3f
env				 	= 	LD_LIBRARY_PATH="/target/lib"

sup_src			:=									\
	framing.cc										\

sup_obj					= $(patsubst %.cc, %.o, $(sup_src))

all							: main.exe

main.exe				: $(sup_obj) main.o
	$(env) $(CXX) $^ -o $@ $(LDFLAGS)

main.o					: main.cc
	$(CXX) $(CFLAGS) -c $< -o $@

$(sup_obj)			: %.o :	%.cc
	$(CXX) $(CFLAGS) -c $< -o $@

clean						:
	$(RM) main.exe main.o
	$(RM) $(sup_obj)
