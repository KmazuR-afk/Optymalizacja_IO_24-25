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
		lab6();
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

	solution s1 = pen(ff3T1, sol, 1, 2, epsilon, Nmax, a[0]);
	int ile_c = solution::f_calls;
	solution::clear_calls();
	solution s2 = pen(ff3T2, sol, 1, 0.5, epsilon, Nmax, a[0]);
	int ile_c2 = solution::f_calls;
	solution::clear_calls();
	cout << s1.x(0) << "  " << s1.x(1) << "  " << s1.y << "  " << ile_c <<endl;
	cout << s2.x(0) << "  " << s2.x(1) << "  " << s2.y << "  " << ile_c2 << endl;

	//Plik CSV- testowanie funkcji celu
	ofstream TestToFile("./test_3.csv");
	
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

	TestToFile.close();
	//////////////////////////////////SPRAWDZENIE_WG_INSTRUKCJI////////////////////////////////////////////////
	/*matrix Y0;
	matrix* Y;
	Y0 = matrix(4, new double[4] {0,5, 100, 0});
	Y = solve_ode(df3R, 0, 0.01, 7, Y0, NAN, 10);
	cout << Y[1] << endl;//po przeskrolowaniu widzimy 21,53; 7,65557; 50,1628; -20,885; oraz 41,4079; 9,43686; -0,00666495; -22,3301;*/
	//////////////////////////////////Symulacja optymalizacja////////////////////////////////////////////////
	ofstream Symulacje("sym_3.csv");
	matrix x0(2,new double[2]{-10,-15});
	solution sv=pen(ff3R,x0,1,2,epsilon,Nmax,NAN,c);
	Symulacje<<sv<<endl;
	Symulacje.close();
	
	Symulacje.open("sym_4.csv");
	matrix Y0;
	matrix* Y;
	Y0 = matrix(4, new double[4] {0,sv.x(0), 100, 0});
	Y = solve_ode(df3R, 0, 0.01, 7, Y0, NAN, sv.x(1));
	Symulacje<<Y[1]<<endl;
	Symulacje.close();

}

void lab4()
{
	double epsilon = 0.001;
	int Nmax = 10000;
	double s[3] = { 0.05, 0.12, -1.0};
	solution sd, cg, newton;
	int sd_f_calls, sd_g_calls, cg_f_calls, cg_g_calls, newton_f_calls, newton_g_calls, newton_H_calls;
	matrix x0;


	ofstream myfile;
	myfile.open("lab_4_t1.csv");
	for(int j=0;j<3;j++){
		for (int i = 0; i < 100; i++) {
			double x_[] = { (double)((rand() % 200) / 10.) - 10., (double)((rand() % 200) / 10.) - 10. };
			x0 = matrix(2, x_);

			sd = SD(ff4T, gf4T, x0, s[j], epsilon, Nmax);
			sd_f_calls = solution::f_calls;
			sd_g_calls = solution::g_calls;

			cg = CG(ff4T, gf4T, x0, s[j], epsilon, Nmax);
			cg_f_calls = solution::f_calls;
			cg_g_calls = solution::g_calls;

			newton = Newton(ff4T, gf4T, Hf4T, x0, s[j], epsilon, Nmax);
			newton_f_calls = solution::f_calls;
			newton_g_calls = solution::g_calls;
			newton_H_calls = solution::H_calls;


			myfile << x0(0) << ";" << x0(1) << ";"
				<< sd.x(0) << ";" << sd.x(1) << ";" << sd.y(0) << ";" << sd_f_calls << ";" << sd_g_calls << ";"
				<< cg.x(0) << ";" << cg.x(1) << ";" << cg.y(0) << ";" << cg_f_calls << ";" << cg_g_calls << ";"
				<< newton.x(0) << ";" << newton.x(1) << ";" << newton.y(0) << ";" << newton_f_calls << ";" << newton_g_calls << ";" << newton_H_calls << endl;
			
		}
	}
	/*double x_[] = { (double)((rand() % 200) / 10.) - 10., (double)((rand() % 200) / 10.) - 10. };
	x0 = matrix(2, x_);
	myfile<<"\t\tsd\n";
	sd = SD(ff4T, gf4T, x0, s[2], epsilon, Nmax,NAN,NAN,&myfile);
	myfile<<"\t\tgradienty\n";
	cg = CG(ff4T, gf4T, x0, s[2], epsilon, Nmax,NAN,NAN,&myfile);
	myfile<<"\t\tNewton\n";
	newton = Newton(ff4T, gf4T, Hf4T, x0, s[2], epsilon, Nmax,NAN,NAN,&myfile);
	myfile.close();*/

	//Część rzeczywista
	//Test porpawności funkcji rzeczywistych
	/*matrix X(3, 100);
	ifstream daneX("XData.txt");
	daneX >> X;

	matrix Y(1, 100);
	ifstream daneY("YData.txt");
	daneY >> Y;

	matrix theta = matrix(3, 1, 0.0);
	matrix J = ff4R(theta, X, Y);
	cout << J << endl;

	matrix dJ = gf4R(theta, X, Y);
	cout << endl << dJ << endl;*/

	//Rozwiązanie problemu rzeczywistego
	/*matrix X(3, 100);
	ifstream daneX("XData.txt");
	daneX >> X;

	matrix Y(1, 100);
	ifstream daneY("YData.txt");
	daneY >> Y;

	matrix theta = matrix(3, 1, 0.0);
	double epsilon = 0.000001;
	int Nmax = 1000000;

	
	double krok[] = {0.01, 0.001, 0.0001};
	
	solution rozw;
	rozw = CG(ff4R, gf4R, theta, krok[2], epsilon, Nmax, X, Y);
	cout << rozw.x(0) << ";" << rozw.x(1) << ";" << rozw.x(2) << ";" << rozw.y(0) << ";" << solution::g_calls;
	//Dolicz prawdopodobieństwo
	double prawd = Ptheta(rozw.x, X, Y);
	cout << ";" << prawd << endl;*/
}

void lab5()
{
	////Parametry dla funkcji
 //   double epsilon = 1e-4;      
 //   int Nmax = 10000;

 //   matrix x0 = matrix(2, 1, 0.0);
	//matrix x1 = matrix(2, 1, 0.0);
	//matrix x2 = matrix(2, 1, 0.0);

 //   //Parametry funkcji celu
	//double wart_a[3] = {1.0, 10.0, 100.0};

 //   matrix ud_a1(2, 1);             
 //   ud_a1(0) = 0; //w
 //   ud_a1(1) = wart_a[0]; //a

	//matrix ud_a2(2, 1);
	//ud_a2(0) = 0; //w
	//ud_a2(1) = wart_a[1]; //a  

	//matrix ud_a3(2, 1);
	//ud_a3(0) = 0; //w
	//ud_a3(1) = wart_a[2]; //a  
 // 
	//ofstream file1("./wyniki_a1.csv");
	//ofstream file2("./wyniki_a2.csv");
	//ofstream file3("./wyniki_a3.csv");

	//for (double w = 0; w < 1.01; w = w + 0.01) {
	//	//Wartości wspólne
	//	ud_a1(0) = w;
	//	ud_a2(0) = w;
	//	ud_a3(0) = w;
	//	cout << w << endl;

	//	x0 = rand_mat(2,1);
	//	x0(0) = 20 * x0(0) - 10;
	//	x0(1) = 20 * x0(1) - 10;
	//	x1 = x0;
	//	x2 = x0;

	//	solution rozw = Powell(ff5T, x0, epsilon, Nmax, ud_a1);
	//	file1 << x0(0) << ";" << x0(1) << ";" << rozw.x(0) << ";" << rozw.x(1) << ";" << rozw.y(0) << ";" << rozw.y(1) << ";" << solution::f_calls << endl;
	//	solution::clear_calls();

	//	solution rozw1 = Powell(ff5T, x1, epsilon, Nmax, ud_a2);
	//	file2 << x1(0) << ";" << x1(1) << ";" << rozw1.x(0) << ";" << rozw1.x(1) << ";" << rozw1.y(0) << ";" << rozw1.y(1) << ";" << solution::f_calls << endl;
	//	solution::clear_calls();

	//	solution rozw2 = Powell(ff5T, x2, epsilon, Nmax, ud_a3);
	//	file3 << x2(0) << ";" << x2(1) << ";" << rozw2.x(0) << ";" << rozw2.x(1) << ";" << rozw2.y(0) << ";" << rozw2.y(1) << ";" << solution::f_calls << endl;
	//	solution::clear_calls();
	//}

	//file1.close();
	//file2.close();
	//file3.close();

	
	//Parametry dla funkcji
	//dobrane przez nas
	double epsilon = 1e-4;
	int Nmax = 10000;
	ofstream file("./wyniki_r.csv");

	matrix ud_a(1, 1);
	//ud_a(0) = 0; //w

	

	for (double w = 0; w < 1.01; w = w + 0.01) {
		//nadawaj kolejne wartości w co 0.01
		ud_a(0) = w;

		matrix x0 = matrix(2, 1, 0.0);
		//losowanie zmiennych l oraz d
		x0 = rand_mat(2, 1);
		//l w przedziale od 0.2 do 1.0 (w metrach)
		x0(0) = 0.8 * x0(0) + 0.2;
		//d w przedziale od 0.01 do 0.05 (w metrach)
		x0(1) = 0.04 * x0(1) + 0.01;

		//wywołanie oraz zapis do pliku
		solution rozw = Powell(ff5R, x0, epsilon, Nmax, ud_a);
		file << x0(0) << ";" << x0(1) << ";" << rozw.x(0) << ";" << rozw.x(1) << ";" << rozw.y(0) << ";" << rozw.y(1) << ";" << solution::f_calls << endl;
		solution::clear_calls();
	}

	file.close();
}

void lab6()
{
	double sigma[5] = { 0.01, 0.1, 1, 10, 100 };
	solution EAf;
	int N = 2; 
	matrix lb(N, 1), ub(N,1);
	lb = matrix(2, 1, -5);
	ub = matrix(2, 1, 5);

	int mi = 20;
	int lambda = 40;
	double epsilon = 1e-5;
	int Nmax = 10000;

	// matrix sigma(2, 1, 1);


	matrix saving(4, 1, 0.);
	ofstream myfile;
	myfile.open("optymalizacja_5_t1.csv");
	for (int a = 0; a < 5; a++) {
		for (int i = 0; i < 100; i++) {
			EAf = EA(ff6T, N, lb, ub, mi, lambda, sigma[a], epsilon, Nmax);
			saving(0) = EAf.x(0);
			saving(1) = EAf.x(1);
			saving(2) = EAf.y(0);
			saving(3) = solution::f_calls;
			string min;
			if (solution::f_calls > 10000)
				min = "NIE";
			else
				min = "TAK";
			myfile << trans(saving) << ";" << min << "\n";
			cout << EAf << endl;
		}
	}
	myfile.close();
}

