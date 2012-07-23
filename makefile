CC=g++
LOP=-O2 -o

MAIN=./Library/lfm
TAG=lfm


$(MAIN).o :
	$(CC) $(LOP) $(TAG) $(MAIN).cpp


