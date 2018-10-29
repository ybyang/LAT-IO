all:
	g++ -fopenmp -I.. main.C -lz
	./a.out ../tests/test.lat | head
#	diff test.lat*

clean:
	-rm a.out
