#pragma once

#include"ode_solver.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);

matrix ff1T(matrix, matrix = NAN, matrix = NAN);
matrix ff1L(matrix, matrix = NAN, matrix = NAN);
matrix ff1R(matrix, matrix = NAN, matrix = NAN);
matrix df1R(double, matrix, matrix = NAN, matrix = NAN);

matrix ff2T(matrix, matrix = NAN, matrix = NAN);
matrix ff2R(matrix , matrix = NAN,matrix = NAN);
matrix df2R(double, matrix, matrix = NAN, matrix = NAN);

matrix ff3T(matrix);
matrix ff3T1(matrix, matrix = NAN, matrix = NAN);
matrix ff3T2(matrix, matrix = NAN, matrix = NAN);
matrix ff3R(matrix , matrix = NAN,matrix = NAN);
matrix df3R(double, matrix, matrix = NAN, matrix = NAN);

matrix ff4T(matrix, matrix = NAN, matrix = NAN);
matrix gf4T(matrix, matrix = NAN, matrix = NAN);
matrix Hf4T(matrix, matrix = NAN, matrix = NAN);
double h_thet(matrix, matrix);
matrix ff4R(matrix, matrix, matrix);
matrix gf4R(matrix, matrix, matrix);
double Ptheta(matrix, matrix, matrix);