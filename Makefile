# Nazwa pliku wykonywalnego
#w Windowsie
EXECUTABLE = program.exe
#w Linuxie
#EXECUTABLE = program

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
#w Windowsie
	$(EXECUTABLE) && del /Q *.o $(EXECUTABLE)
#w Linuxie
#      ./$(EXECUTABLE) && rm -f *.o $(EXECUTABLE)

# Czyszczenie plików .o i pliku wykonywalnego
clean:
#w windowsie
	del /Q *.o *.csv $(EXECUTABLE)
#w Linuxie
#       rm -f *.o *.csv $(EXECUTABLE)
# Phony targets to avoid conflicts with files of the same name
.PHONY: all clean run
