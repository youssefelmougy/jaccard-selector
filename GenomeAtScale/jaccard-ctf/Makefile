include config.mk

FCXX=$(CXX) $(CXXFLAGS)

all: jaccard

file_reader.o: file_reader.h file_reader.cxx
	$(FCXX) -c file_reader.cxx $(INCLUDE_PATH) $(LIB_PATH) $(LIBS)

jaccard: file_reader.o jaccard.cxx
	$(FCXX) -o jaccard jaccard.cxx file_reader.o $(INCLUDE_PATH) $(LIB_PATH) $(LIBS)

clean:
	rm -f file_reader.o jaccard
