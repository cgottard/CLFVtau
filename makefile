all: ev.exe

ev.exe: nominal.o
	g++ -std=c++11 modifiedERC.cpp nominal.o -o ev_mod.exe `root-config --cflags --glibs`
nominal.o: nominal.C
	g++ -std=c++11 -c nominal.C -o nominal.o `root-config --cflags --glibs`

clean:
	rm nominal.o ev_mod.exe
