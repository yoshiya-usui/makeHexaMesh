CXX           = icpc
CC            = icpc
CXXFLAGS      = -O3 \
                -DNDEBUG
DEST          = ./
LIBS          = 
OBJS          = MeshGeneration.o main.o
PROGRAM       = MeshGeneration

all:            $(PROGRAM)

$(PROGRAM):     $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -o $(PROGRAM)

clean:;		rm -f *.o *~ $(PROGRAM)
