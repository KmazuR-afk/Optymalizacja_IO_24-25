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