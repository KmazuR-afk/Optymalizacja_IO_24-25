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

Pominięte zostały opisy wszystkich metod, ale załączone są zdjęcia wykresów i odnośniki do odpowiedniego miejsca w kodzie odnoszącego się do implementacji każdej z metod w pliku [opt_alg.cpp](https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/opt_alg.cpp). Metody używane były zarówno do optymalizacji funkcji matematycznych oraz problemów rzeczywistych - w poniższym README.md zamieszczone są ewentualnie schematy, problemu rzeczywistego, a w celu dogłebniejszej analizy zachęca się do przeczytania zarówno [main.cpp](https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/main.cpp), jak i [user_funs.cpp](https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/user_funs.cpp). Do projektu dołączony jest prosty [Makefile](https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/Makefile) umożliwiający kompilacje i uruchomienie na systemach Windows i Linux - aby przeczytać więcej przejdź do sekcji [Uruchamianie](#uruchamianie)

## Metody bezgradientowe
Najprostsze z metod wyznaczania najmniejszego rozwiązania. Są najbardziej złożone obliczeniowo, ale najmniej skomplikowane - wyznaczają rozwiązania posługując sie wawrtościami funkcji w badanych punktach. W pierwszej kolejności
implementowane były funkcje jednej zmiennej (metody ekspansji, Fibonacciego, interpolacji Lagrange'a), wielu zmiennych (metoda Hooke'a-Jeevesa, Rosenbrocka) i wielu zmiennych z ograniczeniami (metoda sympleksu Neldera-Meada i podstawowe funkcje kary).

### Metoda ekspansji
[Odsyłacz do kodu](https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/opt_alg.cpp#L34-L91)
### Metoda Fibonacciego
[Odsyłacz do kodu](https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/opt_alg.cpp#L93-L140)
### Metoda interpolacji Lagrange'a
[Odsyłacz do kodu](https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/opt_alg.cpp#L142-L216)
### Metoda Hooke'a-Jeevesa
[Odsyłacz do kodu](https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/opt_alg.cpp#L218-L313)
<p align="center">
  <img src="https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/img/Screenshot%20from%202025-02-13%2019-06-03.png" alt="Wykres kroków zawężania przedziału" height="500">
</p>

### Metoda Rosenbrocka
[Odsyłacz do kodu](https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/opt_alg.cpp#L315-L411)
<p align="center">
  <img src="https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/img/Screenshot%20from%202025-02-13%2019-06-21.png" alt="Wykres kroków zawężania przedziału" height="500">
</p>

### Metoda sympleksu Neldera-Meada i podstawowe funkcje kary
[Odsyłacz do kodu](https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/opt_alg.cpp#L413-L561)
<p align="center">
  <img src="https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/img/Screenshot%20from%202025-02-13%2019-28-21.png" alt="Wykres kroków zawężania przedziału" height="500">
</p>

## Metody gradientowe

### Metoda najszybszego spadku
[Odsyłacz do kodu](https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/opt_alg.cpp#L563-L603)
### Metoda gradientów sprzężonych
[Odsyłacz do kodu](https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/opt_alg.cpp#L606-L654)
### Metoda Newtona
[Odsyłacz do kodu](https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/opt_alg.cpp#L656-L697)
<p align="center">
  <img src="https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/img/CDSDN1.png" alt="Wykres kroków zawężania przedziału" height="500">
</p>
<p align="center">
  <img src="https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/img/CDSDN2.png" alt="Wykres kroków zawężania przedziału" height="500">
</p>
<p align="center">
  <img src="https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/img/CDSDN3.png" alt="Wykres kroków zawężania przedziału" height="500">
</p>

## Optymalizacja wielokryterialna

### Metoda Powell'a
[Odsyłacz do kodu](https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/opt_alg.cpp#L741-L791)
<p align="center">
  <img src="https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/img/Screenshot%20from%202025-02-13%2019-42-24.png" alt="Wykres kroków zawężania przedziału" height="500">
</p>

## Algorytm ewolucyjny
[Odsyłacz do kodu](https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/opt_alg.cpp#L794-L887)
<p align="center">
  <img src="https://github.com/KmazuR-afk/Optymalizacja_IO_24-25/blob/main/img/Screenshot%20from%202025-02-13%2019-12-52.png" alt="Wykres kroków zawężania przedziału" height="500">
</p>

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
