# Nazwa pliku wykonywalnego, zależna od systemu operacyjnego
ifeq ($(OS),Windows_NT)
    EXECUTABLE = program.exe
    RM = del /Q
    RUN_CMD = $(EXECUTABLE)
    CLEAR_CMD = cls
else
    EXECUTABLE = program
    RM = rm -f
    RUN_CMD = ./$(EXECUTABLE)
    CLEAR_CMD = clear
endif

# Kompilator
CXX = g++

# Flagi kompilatora
CXXFLAGS = -w -g

# Pliki źródłowe
SRCS = main.cpp \
       matrix.cpp \
       ode_solver.cpp \
       opt_alg.cpp \
       solution.cpp \
       user_funs.cpp

# Pliki nagłówkowe (opcjonalne)
HEADERS = matrix.h \
          ode_solver.h \
          opt_alg.h \
          solution.h \
          user_funs.h

# Pliki obiektowe
OBJS = $(SRCS:.cpp=.o) 

# Zasada domyślna - kompilacja pliku wykonywalnego
all: $(EXECUTABLE)

# Reguła kompilacji pliku wykonywalnego
$(EXECUTABLE): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ && $(CLEAR_CMD)

# Reguła kompilacji plików .o z plików .cpp
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $<

# Dodanie reguły uruchamiania programu
run: $(EXECUTABLE)
	$(RUN_CMD) && $(RM) *.o $(EXECUTABLE)

# Czyszczenie plików .o i pliku wykonywalnego
clean:
	$(RM) *.o *.csv $(EXECUTABLE)

# Phony targets to avoid conflicts with files of the same name
.PHONY: all clean run
