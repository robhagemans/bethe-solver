
libs:		exception.o basics.o bethe-takahashi.o strip.o square.o iso-solver.o configuration.o

clean:
		rm *.o 
		
basics.o:	basics.cc  basics.h 
		g++  -c   basics.cc  -O3 -o basics.o -D_GLIBCXX_DEBUG


exception.o: 	exception.cc exception.h
		g++  -c   exception.cc  -O3  -o exception.o -D_GLIBCXX_DEBUG
		
square.o:       square.cc exception.h basics.h
		g++  -c   square.cc -O3 -o square.o -D_GLIBCXX_DEBUG


configuration.o: configuration.cc configuration.h
		g++ -c configuration.cc -O3 -o configuration.o -D_GLIBCXX_DEBUG


strip.o:	strip.cc strip.h  basics.h configuration.h
		g++  -c   strip.cc  -O3 -o strip.o -D_GLIBCXX_DEBUG

bethe-takahashi.o:   bethe-takahashi.cc bethe-takahashi.h square.h exception.h strip.h basics.h 
		g++  -c   bethe-takahashi.cc  -O3 -o bethe-takahashi.o -D_GLIBCXX_DEBUG

		

iso-solver.o:	iso-solver.cc iso-solver.h bethe-takahashi.h strip.h exception.h basics.h
		g++ -c iso-solver.cc -O3 -o iso-solver.o -D_GLIBCXX_DEBUG



