/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#define _USE_MATH_DEFINES
#include"opt_alg.h"
#include"math.h"
#include<cstdlib>
#include<cstdio>
#include <cmath>

#define PI 3.14


void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		lab3();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	system("pause");
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;
	int Nmax = 10000;
	matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
	solution opt;
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Wahadlo
	Nmax = 1000;
	epsilon = 1e-2;
	lb = 0;
	ub = 5;
	double teta_opt = 1;
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(opt.x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
	ofstream Sout("symulacja_lab0.csv");
	Sout << hcat(Y[0], Y[1]);
	Sout.close();
	Y[0].~matrix();
	Y[1].~matrix();
}

void lab1()
{

	matrix ud1(1, 1, -100), ud2(1, 1, 100);
	double* eskp = new double[2];
	double x0 = rand() % 201 - 100, d = 1.0, alfa =5.29;
	int Nmax = 10000;
	int l = 0;
	double eps = 1e-2;
	double gamma = 1e-200;
	solution fibo;
	solution lagi;
	ofstream expToFile("./eksp.csv");
	ofstream fibToFile("./fib.csv");
	ofstream lafToFile("./lag.csv");
	ofstream slagToFile("./slag.csv");
	ofstream sfimToFile("./sfib.csv");

	for (int i = 0; i < 100; i++) {
		eskp = expansion(ff1L, x0, d, alfa, Nmax);
		expToFile << x0 << ";" << eskp[0] << ";" << eskp[1] << ";" << solution::f_calls << endl;
		solution::clear_calls();
		x0 = rand() % 201 - 100;

		fibo = fib(ff1L, eskp[0], eskp[1], eps, NULL, NULL);
		fibToFile << m2d(fibo.x) << ";" << m2d(fibo.y) << ";" << solution::f_calls << endl;
		solution::clear_calls();

		lagi = lag(ff1L, eskp[0], eskp[1], eps, gamma, Nmax, NULL, NULL);
		lafToFile << m2d(lagi.x) << ";" << m2d(lagi.y) << ";" << solution::f_calls << endl;
		solution::clear_calls();
	}
	expToFile.close();
	fibToFile.close();
	lafToFile.close();

	//zad1_bezpara
	/*fibo = fib(ff1L, -100, 100, eps);
	lagi = lag(ff1L, -100, 100, eps, gamma, Nmax);
	cout << m2d(fibo.x) << ";" << m2d(fibo.y) << ";" << solution::f_calls << endl;
	cout << m2d(lagi.x) << ";" << m2d(lagi.y) << ";" << solution::f_calls << endl;*/
	//symulacje
	double epss = 0.00000000001;
	double gammas = 0.000001;
	fibo = fib(ff1R, 0.0001, 0.01,epss,ud1,ud2);
	//cout << m2d(fibo.x) << ";" << m2d(fibo.y) << ";" << solution::f_calls << endl;
	lagi= lag(ff1R, 0.0001, 0.01, epss,gammas, Nmax, ud1, ud2);
	//cout << m2d(lagi.x) << ";" << m2d(lagi.y) << ";" << solution::f_calls << endl;
	matrix Y0 = matrix(3, new double[3] {5, 1, 20});
	matrix* Yf = solve_ode(df1R, 0, 1, 2000, Y0, ud1, fibo.x(0));
	sfimToFile << Yf[1] << endl;
	matrix* Yl= solve_ode(df1R, 0, 1, 2000, Y0, ud1, lagi.x(0));
	slagToFile << Yl[1] << endl;
	sfimToFile.close();
	slagToFile.close();
}

void lab2()
{
	srand(time(NULL));

	int Nmax=1e6;
	matrix x0(2,1);
	double s = 0.01 + (rand() % 200) / 100.0;

	//s = 0.59; //do wykresu

	double epsilon = 1e-3, alpha_hj = 0.6, alpha_rs = 1.3, beta = 0.4;
	matrix sol = matrix(2,1,0.0);
	matrix s_wek = matrix(2, 1, s);
	cout << s << endl;

	ofstream HJToFile("./hj.csv");
	ofstream RSToFile("./rs.csv");

	/*for (int i = 0; i < 100; i++) {       //tabele 1-2
		sol = 2 * rand_mat(2, 1) - 1;

		solution solv_hj = HJ(ff2T, sol, s, alpha_hj, epsilon, Nmax);
		int ck = solution::f_calls;
		HJToFile << sol(0) << ";" << sol(1) << ";" << solv_hj.x(0) << ";" << solv_hj.x(1) << ";" << solv_hj.y(0) << ";" << solution::f_calls << ";" << endl;
		solution::clear_calls();

		solution solv_rs = Rosen(ff2T, sol, s_wek, alpha_rs, beta, epsilon, Nmax);
		ck = solution::f_calls;
		RSToFile << sol(0) << ";" << sol(1) << ";" << solv_rs.x(0) << ";" << solv_rs.x(1) << ";" << solv_rs.y(0) << ";" << solution::f_calls << ";" << endl;
		solution::clear_calls();
		
	}*/

	/*sol(0) = -0.973332;  //do wykresu
	sol(1) = -0.0686693;
	solution solv_hj = HJ(ff2T, sol, s, alpha_hj, epsilon, Nmax);
	solution solv_rs = Rosen(ff2T, sol, s_wek, alpha_rs, beta, epsilon, Nmax);*/

	HJToFile.close();
	RSToFile.close();


	///////////////////////////////////////////////////SYMULACJA///////////////////////////////////////////////

	//x0 = 10 * rand_mat(2, 1);
	x0 = matrix(2, 1, 5.0);
	s = 0.1;
	s_wek = matrix(2, 1, s);
	matrix ud1(2, new double[2] {PI, 0});

	//////////////////////////////////SPRAWDZENIE_WG_INSTRUKCJI////////////////////////////////////////////////
	//matrix Y = ff2R(x0, ud1);
	//cout << Y << endl;
	
	//wydruk w kolejsności x1_0, x2_0, x1_k, x2_k, y_k
	//std::cout<< x0(0) << ";" << x0(1) << ";" << hj_result.x(0) << ";" << hj_result.x(1) << ";" << hj_result.y(0) << ";" << ";min;";


	///////////////////////////////////SYMULACJA_WL//////////////////////////////////////////////////////////
	solution::clear_calls();
	x0 = 10 * rand_mat(2, 1);

	cout << "Rozpoczynam zapisywanie do pliku Symulacja.csv..." << endl;

	ofstream Hooke_RToFile("./Hooke_RToFile.csv");
	ofstream Rosen_BToFile("./Rosen_BToFile.csv");

	s = 0.1;
	matrix s0 = matrix(2, 1, s);

	cout << x0 << endl << endl;

	solution Hooke_R = HJ(ff2R, x0, s, alpha_hj, epsilon, Nmax, ud1);

	matrix Y = ff2R(Hooke_R.x, ud1);
	cout << Hooke_R.x << endl;
	cout << Y << endl;
	cout << solution::f_calls << endl;
	solution::clear_calls();

	matrix Y0 = matrix(2, 1);
	matrix* YzHooke = solve_ode(df2R, 0, 0.1, 100, Y0, ud1, Hooke_R.x);
	Hooke_RToFile << YzHooke[1] << endl;

	solution Rosen_B = Rosen(ff2R, x0, s0, alpha_rs, beta, epsilon, Nmax, ud1);

	Y = ff2R(Rosen_B.x, ud1);
	cout << Rosen_B.x << endl;
	cout << Y << endl;
	cout << solution::f_calls << endl;

	matrix* YzRosen = solve_ode(df2R, 0, 0.1, 100, Y0, ud1, Rosen_B.x);
	Rosen_BToFile << YzRosen[1] << endl;

	Hooke_RToFile.close();
	Rosen_BToFile.close();
}

void lab3()
{
	srand(time(NULL));
	//Część tych parametrow trzeba incjować w pen
	double gamma = 2.0;  //z prezentacji
	double beta = 0.5;  //z prezentacji
	double delta = 0.5;  //z prezentacji

	double alfa = 1.0;  //rozszerzenie sympleksu
	double s = 0.4;  //rozmiar początkowy sympleksu - około 0.1 zakresu poszukwiań
	double c = 1.0;  //współczynnik kary początkowy
	double dc1 = 2.0;  //współczynnik kary zewnętrznej 
	double dc2 = 0.5;  //współczynnik kary wewnętrznej
	double epsilon = 1e-3;  //z instrukcji
	int Nmax = 1e6;

	double a[3];
	a[0] = 4.0;
	a[1] = 4.4934;
	a[2] = 5.0;
	matrix sol = matrix(2, 1, 0.0);

	/*s1 = pen(ff3T1, sol, 1, 2, epsilon, Nmax, a[0]);
	int ile_c = solution::f_calls;
	solution::clear_calls();
	s2 = pen(ff3T2, sol, 1, 0.5, epsilon, Nmax, a[0]);
	int ile_c2 = solution::f_calls;
	solution::clear_calls();
	cout << s1.x(0) << "  " << s1.x(1) << "  " << s1.y << "  " << ile_c <<endl;
	cout << s2.x(0) << "  " << s2.x(1) << "  " << s2.y << "  " << ile_c2 << endl;*/

	//Plik CSV- testowanie funkcji celu
	/*ofstream TestToFile("./test_3.csv");
	
	for (int i = 0; i < 100; i++) {
		sol = rand_mat(2, 1); //losowanie z przedziału [0,1] dla x1 i x2
		double ogr1 = sqrt(pow(a[0], 2) - 1);
		sol(0) = sol(0) * (ogr1 - 1); //x1 od 0 do ogr1 -1 
		sol(0) = 1 + sol(0); //x1 od 1 do ogr1
		double ogr2 = sqrt(pow(a[0], 2) - pow(sol(0), 2)); //dopasuj x2 względem wylosowanego x1
		sol(1) = sol(1) * (ogr2 - 1); //x2 od 0 do ogr2 -1 
		sol(1) = 1 + sol(1); //x2 od 1 do ogr2
		//Koniec losowania, dalej wywołania funkcji

		solution s1, s2;
		s1 = pen(ff3T1, sol, c, dc1, epsilon, Nmax, matrix(a[0]));
		int fc_s1 = solution::f_calls;
		solution::clear_calls();

		s2 = pen(ff3T2, sol, c, dc2, epsilon, Nmax, matrix(a[0]));
		int fc_s2 = solution::f_calls;
		solution::clear_calls();

		TestToFile << sol(0) << ";" << sol(1) << ";" << s1.x(0) << ";" << s1.x(1) << ";" << norm(s1.x) << ";" << s1.y << ";" << fc_s1 << ";" << s2.x(0) << ";" << s2.x(1) << ";" << norm(s2.x) << ";" << s2.y << ";" << fc_s2 << endl;
	}

	TestToFile.close();*/
	//////////////////////////////////SPRAWDZENIE_WG_INSTRUKCJI////////////////////////////////////////////////

	matrix x0 = matrix(2, new double[2]{0,0});
	cout<<x0<<endl;
	solution NM = pen(ff3R, x0, 1., 2, epsilon, Nmax, NAN, c);
	cout << NM << endl;

	matrix Y0;
	matrix* Y;

	Y0 = matrix(4, new double[4] {0, NM.x(0), 100, 0});
	Y = solve_ode(df3R, 0, 0.01, 7, Y0, NAN, NM.x(1));

	cout << Y[1] << endl;



}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}

