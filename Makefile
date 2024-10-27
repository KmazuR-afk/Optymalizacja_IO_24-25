# Nazwa pliku wykonywalnego
EXECUTABLE = program

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
	$(CXX) $(CXXFLAGS) -o $@ $^

# Reguła kompilacji plików .o z plików .cpp
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $<

# Dodanie reguły uruchamiania programu
run: $(EXECUTABLE)
	./$(EXECUTABLE) && rm -f *.o program
# Czyszczenie plików .o i pliku wykonywalnego
clean:
	rm -f *.o *.csv program

# Phony targets to avoid conflicts with files of the same name
.PHONY: all clean run
