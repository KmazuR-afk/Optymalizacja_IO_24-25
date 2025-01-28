# Algorytmy optymalizacyjne

Projekt pisany podczas realizowania przedmiotu Optymalziacja w ramach programu studiów Inżynieria Obliczeniowa na uczelni AGH. Do implementacji została wykorzystane dostarczone przez prowadącego klasy matrix(podstawowe operacje na macierzach), solution 
(przechowująca rozwiązanie funkcji celu i umożliwiająca wyznaczenie jej wartości, gradientu i hesjanu) oraz  ode_solver (prosty solver układów równań) - nie zostały udostępnione w repozytorium, ponieważ nie były przez nas edytowane. Podczas zajęć zaimplemenowaliśmy:

## Spis treści
- [Metody bezgradientowe](#metody-bezgradientowe)
  - [Metoda ekspansji](#metoda-ekspansji)
  - [Metoda Fibonacciego](#metoda-fibonacciego)
  - [Metoda interpolacji Lagrange'a](#metoda-interpolacji-lagrangea)
  - [Metoda Hooke'a-Jeevesa](#metoda-hookea-jeevesa)
  - [Metoda Rosenbrocka](#metoda-rosenbrocka)
  - [Metoda sympleksu Neldera-Meada i podstawowe funkcje kary](#metoda-sympleksu-neldera-meada-i-podstawowe-funkcje-kary)
- [Metody gradientowe](#metody-gradientowe)
  - [Metoda najszybszego spadku](#metoda-najszybszego-spadku)
  - [Metoda gradientów sprzężonych](#metoda-gradient%C3%B3w-sprz%C4%99%C5%BConych)
  - [Metoda Newtona](#metoda-newtona)
- [Optymalizacja wielokryterialna](#optymalizacja-wielokryterialna)
- [Algorytm ewolucyjny](#algorytm-ewolucyjny)
- [Uruchamianie](#uruchamianie)

Pominięte zostały opisy wszystkich metod, ale załączone są zdjęcia wykresów i odnośniki do odpowiedniego miejsca w kodzie odnoszącego się do implementacji każdej z metod w pliku [opt_alg.cpp](https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/opt_alg.cpp). Metody używane były zarówno do optymalizacji funkcji matematycznych oraz problemów rzeczywistych - w poniższym README.md zamieszczone są ewentualnie schematy, problemu rzeczywistego, a w celu dogłebniejszej analizy zachęca się do przeczytania zarówno [main.cpp](https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/main.cpp), jak i [user_funs.cpp](https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/user_funs.cpp). Do projektu dołączony jest prosty [Makefile](https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/makefile) umożliwiający kompilacje i uruchomienie na systemach Windows i Linux
- aby przeczytać więcej przejdź do sekcji [Uruchamianie](#uruchamianie)

## Metody bezgradientowe
Najprostsze z metod wyznaczania najmniejszego rozwiązania. Są najbardziej złożone obliczeniowo, ale najmniej skomplikowane - wyznaczają rozwiązania posługując sie wawrtościami funkcji w badanych punktach. W pierwszej kolejności
implementowane były funkcje jednej zmiennej (metody ekspansji, Fibonacciego, interpolacji Lagrange'a), wielu zmiennych (metoda Hooke'a-Jeevesa, Rosenbrocka) i wielu zmiennych z ograniczeniami (metoda sympleksu Neldera-Meada i podstawowe funkcje kary).

### Metoda ekspansji

### Metoda Fibonacciego

### Metoda interpolacji Lagrange'a

### Metoda Hooke'a-Jeevesa

### Metoda Rosenbrocka

### Metoda sympleksu Neldera-Meada i podstawowe funkcje kary

## Metody gradientowe

### Metoda najszybszego spadku

### Metoda gradientów sprzężonych

### Metoda Newtona

## Optymalizacja wielokryterialna

### Metoda Powell'a

## Algorytm ewolucyjny

## Uruchamianie

W celu przeprowadzenia wyłącznie kompilacji:
```bash
make
```
W celu przeprowadzenia kompilacji, uruchomienia programu i usunięcia po zakończeniu plików *.o:
```bash
make run
```
W celu usunięcia plików *.o i *.csv wygenerowanych w wyniku działania programu:
```bash
make clean
```

Aby zmieniać uruchamiane funkcje należy przeanalizować [main.cpp](https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/main.cpp) i edytować go zgodnie ze swoimi potrzebami.
