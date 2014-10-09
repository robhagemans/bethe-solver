CFLAGS = -D_GLIBCXX_DEBUG -std=c++0x
TARGET = ..


libs:		$(TARGET)/roots.o $(TARGET)/symmetric-string.o $(TARGET)/nonsymmetric-string.o $(TARGET)/central-string.o  $(TARGET)/solver.o $(TARGET)/simplescan.o  $(TARGET)/roots-real.o $(TARGET)/roots-imag.o
#
clean:
		rm $(TARGET)/*.o 
		

$(TARGET)/roots.o:	roots.cc  basics.h roots.h 
		g++  $(CFLAGS) -c   roots.cc  -o $(TARGET)/roots.o   

$(TARGET)/roots-real.o:	roots-real.cc  basics.h roots.h 
		g++  $(CFLAGS) -c   roots-real.cc  -o $(TARGET)/roots-real.o   

$(TARGET)/roots-imag.o:	roots-imag.cc  basics.h roots.h 
		g++  $(CFLAGS) -c   roots-imag.cc  -o $(TARGET)/roots-imag.o   


$(TARGET)/symmetric-string.o: symmetric-string.cc  basics.h roots.h takahashi-string.h symmetric-string.h 
		g++  $(CFLAGS) -c   symmetric-string.cc   -o $(TARGET)/symmetric-string.o  

$(TARGET)/central-string.o: central-string.cc  basics.h roots.h takahashi-string.h symmetric-string.h 
		g++  $(CFLAGS) -c   central-string.cc   -o $(TARGET)/central-string.o  


$(TARGET)/nonsymmetric-string.o: nonsymmetric-string.cc  basics.h roots.h takahashi-string.h nonsymmetric-string.h 
		g++  $(CFLAGS) -c   nonsymmetric-string.cc  -o $(TARGET)/nonsymmetric-string.o  

$(TARGET)/solver.o:	solver.cc  basics.h roots.h takahashi-string.h symmetric-string.h nonsymmetric-string.h solver.h  
		g++  $(CFLAGS) -c   solver.cc   -o $(TARGET)/solver.o  

$(TARGET)/simplescan.o:   simplescan.cc basics.h takahashi-string.h
		g++  $(CFLAGS) -c   simplescan.cc -o $(TARGET)/simplescan.o   		


