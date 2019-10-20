
# CXX=g++
CXX=clang++
AR=ar

FLAGS=-O2 -Wshadow -Wunused -std=c++14 -fpic

GD_INCLUDES_L=-lgd -lpng -lz -lfreetype -lm
#BOOST_INCLUDES_L=-lboost_regex -lboost_program_options
PYTHON_INCLUDES_C=-I/usr/include/python3.7m

HEADERS=auxiliary.hh canvas.hh

OBJECTS=geometry3.o canvas.o
OBJECTS_STANDALONE=trees.o
OBJECTS_LIB=interface.o

OBJECTS_P=$(patsubst %.o, build/%.o, $(OBJECTS))
OBJECTS_STANDALONE_P=$(patsubst %.o, build/%.o, $(OBJECTS_STANDALONE))
OBJECTS_LIB_P=$(patsubst %.o, build/%.o, $(OBJECTS_LIB))

build/%.o: %.cc $(HEADERS) Makefile
	#$(CXX) -DGRAPHICS_DEBUG $(FLAGS) $(PYTHON_INCLUDES_C) -c $< -o $@
	$(CXX) $(FLAGS) $(PYTHON_INCLUDES_C) -c $< -o $@

trees: $(OBJECTS_P) $(OBJECTS_STANDALONE_P) Makefile
	$(CXX) $(OBJECTS_P) $(OBJECTS_STANDALONE_P) $(GD_INCLUDES_L) -o $@

libarttrees.so: $(OBJECTS_P) $(OBJECTS_LIB_P) Makefile
	$(CXX) $(FLAGS) -shared $(OBJECTS_P) $(OBJECTS_LIB_P) $(GD_INCLUDES_L) -o $@


.PHONY: test
test: Makefile test.cc scene.hh auxiliary.hh
	$(CXX) $(FLAGS) test.cc -o test
	./test

# .PHONY: test-circ
# test-circ: Makefile test-circ.cc circ_360.hh
# 	$(CXX) -g -O2 -std=c++14 test-circ.cc -o test-circ
# 	./test-circ

.PHONY: clean distclean
clean:
	rm -f out.png debug*

distclean: clean
	rm -f pano pano-debug test a.out libartpano.so build/*o



