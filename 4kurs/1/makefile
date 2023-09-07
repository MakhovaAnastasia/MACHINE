
all:
# g++ -g - generate debug info 
# g++ -o имя_исполняемого файла - задаем исполняемый файл, который надо указать в настройках 
# запуска отладчика vscode
#	g++ Solve.cpp ReadMatrix.cpp PrintMatrix.cpp main.cpp -g -o Proga.out
#ifndef CallfromVSCode
#	./Proga.out 
#endif
CC=g++
CFLAGS=-c -g -O3 -Wall
LDFLAGS=
SOURCES=main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=1

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@