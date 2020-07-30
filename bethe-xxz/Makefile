libs:		exception.o bethe.o square.o det.o strip.o young.o chain.o base.o state.o xxz-chain.o xxz-base.o xxz-state.o gap-chain.o gap-base.o gap-state.o iso-chain.o iso-base.o iso-state.o scan.o process.o iso-deviated.o main.o quantity.o generic.o recursive.o

parallel:	parallel.o

clean:
		rm *.o 

exception.o: 	exception.cc exception.h
		g++  -c   exception.cc  -O3  -o exception.o
		
main.o: 	main.cc exception.h		
		g++  -c   main.cc  -O3 -o main.o
		
bethe.o:	bethe.cc  bethe.h exception.h
		g++  -c   bethe.cc  -O3 -o bethe.o

square.o:	square.cc square.def square.h exception.h bethe.h
		g++  -c   square.cc  -O3  -o square.o

matrix.o:	matrix.cc matrix.h exception.h
		g++  -c   matrix.cc  -O3  -o matrix.o

det.o:		det.cc det.def det.h square.h square.def bethe.h
		g++  -c   det.cc  -O3 -o det.o

young.o:	young.cc young.h exception.h
		g++  -c   young.cc  -O3 -o young.o

strip.o:	strip.cc strip.h  bethe.h base.h
		g++  -c   strip.cc  -O3 -o strip.o
		
chain.o:	chain.cc chain.h  bethe.h 
		g++  -c   chain.cc  -O3 -o chain.o
		
base.o: 	base.cc base.h bethe.h chain.h young.h
		g++  -c   base.cc  -O3 -o base.o
		
state.o:	state.cc state.h  bethe.h strip.h strip.cc square.h det.h chain.h base.h
		g++  -c   state.cc  -O3 -o state.o

xxz-chain.o:	xxz-chain.cc xxz-chain.h exception.h  bethe.h chain.h
		g++  -c   xxz-chain.cc  -O3 -o xxz-chain.o
		
xxz-base.o: 	xxz-base.cc xxz-base.h base.h
		g++  -c   xxz-base.cc  -O3 -o xxz-base.o
		
xxz-state.o:	xxz-state.cc xxz-state.h  bethe.h state.h det.h strip.h strip.cc
		g++  -c   xxz-state.cc  -O3 -o xxz-state.o

iso-chain.o:	iso-chain.cc iso-chain.h exception.h  bethe.h chain.h
		g++  -c   iso-chain.cc  -O3 -o iso-chain.o
		
iso-base.o: 	iso-base.cc iso-base.h base.h
		g++  -c   iso-base.cc  -O3 -o iso-base.o
		
iso-state.o:	iso-state.cc iso-state.h   bethe.h state.h det.h strip.h strip.cc
		g++  -c   iso-state.cc  -O3 -o iso-state.o

iso-deviated.o:	iso-deviated.cc iso-deviated.h iso-state.h strip.h strip.cc
		g++  -c   iso-deviated.cc  -O3 -o iso-deviated.o


gap-chain.o:	gap-chain.cc gap-chain.h exception.h  bethe.h chain.h
		g++  -c   gap-chain.cc  -O3 -o gap-chain.o
		
gap-base.o: 	gap-base.cc gap-base.h base.h
		g++  -c   gap-base.cc  -O3 -o gap-base.o
		
gap-state.o:	gap-state.cc gap-state.h  bethe.h state.h det.h strip.h strip.cc
		g++  -c   gap-state.cc  -O3 -o gap-state.o

generic.o: 	generic.cc generic.h gap-state.h gap-base.h gap-chain.h iso-state.h iso-base.h iso-chain.h xxz-state.h xxz-base.h xxz-chain.h bethe.h matrix.h
		g++  -c   generic.cc  -O3 -o generic.o

quantity.o: quantity.cc quantity.h generic.h gap-state.h gap-base.h gap-chain.h iso-state.h iso-base.h iso-chain.h xxz-state.h xxz-base.h xxz-chain.h bethe.h matrix.h
		g++  -c   quantity.cc  -O3 -o quantity.o

scan.o: 	scan.cc scan.h quantity.h generic.h gap-state.h gap-base.h gap-chain.h iso-state.h iso-base.h iso-chain.h xxz-state.h xxz-base.h xxz-chain.h bethe.h matrix.h
		g++  -c   scan.cc  -O3 -o scan.o


recursive.o: 	recursive.cc recursive.h scan.h quantity.h generic.h gap-state.h gap-base.h gap-chain.h iso-state.h iso-base.h iso-chain.h xxz-state.h xxz-base.h xxz-chain.h bethe.h matrix.h
		g++  -c   recursive.cc  -O3 -o recursive.o


process.o:	process.cc process.h scan.h
		g++  -c   process.cc  -O3 -o process.o

parallel.o:	parallel.cc parallel.h scan.h
		g++ -c   parallel.cc  -O3 -o parallel.o

frontend-testscan: scan.o  xxz-state.o  xxz-base.o  xxz-chain.o  iso-state.o iso-base.o  iso-chain.o state.o base.o  chain.o strip.o det.o square.o  bethe.o  exception.o 
		g++   -O3 -g -pg  frontend-testscan.cc exception.o strip.o bethe.o det.o square.o  chain.o base.o state.o   xxz-chain.o  xxz-base.o xxz-state.o   iso-chain.o  iso-base.o   iso-state.o scan.o  -o testscan
		
