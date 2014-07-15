LFLAGS = -lm -lfftw3  

main: main.cpp 
	g++ $(LFLAGS) -o $@ $(MOBLIB) $(SUBFILES) main.cpp $(LIBS)
