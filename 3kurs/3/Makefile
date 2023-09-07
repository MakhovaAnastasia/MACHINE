
#all:
# g++ -g - generate debug info 
# g++ -o имя_исполняемого файла - задаем исполняемый файл, который надо указать в настройках 
# запуска отладчика vscode
#	g++ Solve.cpp ReadMatrix.cpp PrintMatrix.cpp main.cpp -g -o Proga.out

CC=g++ -pthread
CFLAGS=-c -g -Wall  -O3
LDFLAGS=
SOURCES=Solve.cpp synchronize.cpp get_time.cpp ReadMatrix.cpp PrintMatrix.cpp main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=hello

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
