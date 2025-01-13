#include"user_funs.h"
#include <cmath>
#define M_PI 3.14159265359

matrix ff0T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2] { m2d(x), 0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
	int n = get_len(Y[0]);
	double teta_max = Y[1](0, 0);
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);
	double m = 1, l = 0.5, b = 0.5, g = 9.81;
	double I = m * pow(l, 2);
	dY(0) = Y(1);
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;
	return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	y = pow(x, 2);
	return y;
}
matrix ff1L(matrix x, matrix ud1, matrix ud2) {
	matrix y = -cos(0.1 * x(0)) * exp(-pow(0.1 * x(0) - 2 * M_PI, 2)) + 0.002 * pow(0.1 * x(0), 2);
	return y;

}
matrix ff1R(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	matrix Y0 = matrix(3, new double[3] { 5, 1, 20 });
	matrix* Y = solve_ode(df1R, 0, 1, 2000, Y0, ud1, x);
	int n = get_len(Y[0]);
	double max = Y[1](0, 2);
	for (int i = 0; i < n; i++) {
		if (max < Y[1](i, 2))
			max = Y[1](i, 2);
	}

	y = abs(max - 50);
	return y;
}

matrix df1R(double t, matrix Y, matrix ud1, matrix ud2) {
	double a = 0.98;
	double b = 0.63;
	double g = 9.81;

	double P_A = 0.5;
	double V_A = 5.0;
	double T_A = 90;

	double P_B = 1;
	double D_B = 0.00365665;

	double F_in = 0.010;
	double T_in = 20;

	double outflow_A, outflow_B;

	if (Y(0) > 0) {
		outflow_A = a * b * m2d(ud2) * sqrt(2.0 * g * Y(1) / P_B);
	}
	else {
		outflow_A = 0.0;
	}
	if (Y(1) > 0) {
		outflow_B = a * b * D_B * sqrt(2.0 * g * Y(1) / P_B);
	}
	else {
		outflow_B = 0.0;
	}

	matrix dY = matrix(3, 1);
	dY(0) = -outflow_A;
	dY(1) = (outflow_A + F_in - outflow_B);
	dY(2) = (F_in / Y(1)) * (T_in - Y(2)) + (outflow_A / Y(1)) * (T_A - Y(2));
	return dY;
}

matrix ff2T(matrix x, matrix ud1, matrix ud2) {
	matrix y = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * M_PI * x(0)) - cos(2.5 * M_PI * x(1)) + 2;
	return y;
}

matrix ff2R(matrix x, matrix ud1, matrix ud2) {
	matrix Y, Msqr;
	matrix Y0(2, 1);
	matrix* dY = solve_ode(df2R, 0, 0.1, 100, Y0, ud1, x);
	int len = get_len(dY[0]);

	for (int i = 0; i < len; i++) {
		Msqr = pow(x(0) * (ud1(0) - dY[1](i, 0)) + x(1) * (ud1(1) - dY[1](i, 1)), 2);
		Y = Y + (10 * pow(ud1(0) - dY[1](i, 0), 2) + pow(ud1(1) - dY[1](i, 1), 2) + Msqr) * 0.1;
	}
	return Y;
}

matrix df2R(double t, matrix y, matrix ud1, matrix ud2) {
	int l = 1;
	int mr = 1;
	int mc = 5;
	double b = 0.5;
	double I = (mr * pow(l, 2) / 3) + mc * pow(l, 2);
	matrix dY(2, 1);
	dY(0) = y(1);
	dY(1) = (ud2(0) * (ud1(0) - y(0)) + ud2(1) * (ud1(1) - y(1)) - b * y(1)) / I;
	return dY;
}

matrix ff3T(matrix x) {
	matrix out = (sin(3.14 * sqrt(pow(x(0) / 3.14, 2) + pow(x(1) / 3.14, 2))) / (3.14 * sqrt(pow(x(0) / 3.14, 2) + pow(x(1) / 3.14, 2))));
	return out;
}

//Kara zewn�trzna
matrix ff3T1(matrix x, matrix ud1, matrix ud2) {
	matrix out = ff3T(x);
	//nie daje max bo warunki zapewniaja mi > 0
	if (-x(0) + 1 > 0) {
		out = out + ud2 * pow(-x(0) + 1, 2);
	}
	if (-x(1) + 1 > 0) {
		out = out + ud2 * pow(-x(1) + 1, 2);
	}
	if (norm(x) - ud1 > 0) {
		out = out + ud2 * pow(norm(x) - ud1, 2);
	}
	return out;
}

//Kara wewn�trzna
matrix ff3T2(matrix x, matrix ud1, matrix ud2) {
	matrix out = ff3T(x);
	if (-x(0) + 1 <= 0) {
		out = out - ud2 / (-x(0) + 1);
	}
	else {
		out = 1000000000;
	}

	if (-x(1) + 1 <= 0) {
		out = out - ud2 / (-x(1) + 1);
	}
	else {
		out = 1000000000;
	}

	if (norm(x) - ud1 <= 0) {
		out = out - ud2 / (norm(x) - ud1); //?
	}
	else {
		out = 1000000000;
	}

	return out;
}

matrix ff3R(matrix x, matrix ud1, matrix ud2) {
	matrix y;
	matrix Y0(4,new double[4] {0,x(0),100,0});
	matrix* Y=solve_ode(df3R,0,0.01,7,Y0,ud1,x(1));
	int i50=0;
	int i0=0;
	for (int i =0;i<get_len(Y[0]);i++){
		if(abs((Y[1](i,2)-50))<(abs(Y[1](i50,2)-50)))	i50=i;
		if(abs(Y[1](i,2))<abs(Y[1](i0,2)))	i0=i;
	}
	y=-Y[1](i0,0);
	if(abs(x(0))-10>0){//kara dla prędkości
		y=y+ud2*pow(abs(x(0))-10,2);
	}
	if(abs(x(1))-15>0){//kara dla omegi
		y=y+ud2*pow(abs(x(1))-15,2);
	}
	if((abs(Y[1](i50,0))-5)-0.5>0){//kara za nie wpadnięcie do kosza
		y=y+ud2*pow((abs(Y[1](i50,0))-5)-0.5,2);
	}

	return y;
}

matrix df3R(double t, matrix y, matrix ud1, matrix ud2) {
	matrix dY(4,1);//dx,d^2x,dy...
	double g=9.81;//m/s^2
	double m=0.6;//kg
	double C=0.47;
	double ro=1.2;//kg/m^3
	double r=0.12;//m
	double S=M_PI*pow(r,2);
	double Dx=0.5*C*ro*S*y(1)*abs(y(1));
	double Dy=0.5*C*ro*S*y(3)*abs(y(3));
	double Fmx=ro*y(3)*ud2(0)*M_PI*pow(r,3);
	double Fmy=ro*y(1)*ud2(0)*M_PI*pow(r,3);
	dY(0)=y(1);
	dY(2)=y(3);
	dY(1)=-(Dx+Fmx)/m;
	dY(3)=-(Dy+Fmy)/m-g;

	return dY;
}

matrix ff4T(matrix x,matrix ud1,matrix ud2) {
	matrix y;
	if(isnan(ud2(0,0))){
		y=pow((x(0)+2*x(1)-7),2)+pow((2*x(0)+x(1)-5),2);
	}
	else{
		y = ff4T(ud2[0] + x * ud2[1], ud1);
	}
	return y;
}

matrix gf4T(matrix x,matrix ud1,matrix ud2) {
	matrix Xprim(2,1);
	Xprim(0)=10.0*x(0)+8.0*x(1)-34.0;
	Xprim(1)=8.0*x(0)+10.0*x(1)-38.0;
	return Xprim;
}

matrix Hf4T(matrix x, matrix ud1, matrix ud2)
{
	matrix y(2, 2);
	y(0, 0) = 10.0;
	y(0, 1) = 8.0;
	y(1, 1) = 10.0;
	y(1.0) = 8.0;
	return y;
}

double h_thet(matrix thet, matrix x) {
	double ht = (1.0 / (1.0 + exp(m2d(-trans(thet) * x))));
	return ht;
}

matrix ff4R(matrix x, matrix ud1, matrix ud2) {
	int* rozmiar = get_size(ud1); //rozmiar[1] = m
	matrix J;
	double ht = 1;
	for (int i = 0; i < rozmiar[1]; i++) {
		ht = h_thet(x, ud1[i]); //wyliczenie h_theta(x)
		J = J + ud2(0, i) * log(ht) + (1 - ud2(0, i)) * log(1 - ht); //operacja sumy do wyliczenia J(theta), w ud2(0,i) znajduje się i-ty y
	}
	J = -J / rozmiar[1];
	return J;
}

matrix gf4R(matrix x, matrix ud1, matrix ud2) {
	int* rozmiar = get_size(ud1); //rozmiar[1] = m
	matrix dJ;
	double ht = 1;
	for (int i = 0; i < rozmiar[1]; i++) {
		ht = h_thet(x, ud1[i]); //wyliczenie h_theta(x)
		dJ = dJ + (ht - ud2(0, i)) * ud1[i];  //operacja sumy do wyliczenia pochodnej cząstkowej J(theta)
	}
	dJ = dJ / rozmiar[1];
	return dJ;
}

double Ptheta(matrix solwx, matrix ud1, matrix ud2) {
	int zakw = 0;
	int* rozmiar = get_size(ud1);
	for (int i = 0; i < rozmiar[1]; i++) {
		double pred = h_thet(solwx, ud1[i]);
		int wynik = -1;
		if (pred >= 0.5) {
			wynik = 1;
		}
		else {
			wynik = 0;
		}

		if (wynik == ud2[i]) {
			zakw++;
		}
	}
	return (double)zakw / rozmiar[1];
}

matrix ff5T(matrix x, matrix ud1, matrix ud2) {
    matrix y;
    if (isnan(ud2(0, 0))) {
        y = matrix(2, 1);
        y(0) = ud1(1) * (pow(x(0) - 2, 2) + pow(x(1) - 2, 2)); // f1
        y(1) = 1.0 / ud1(1) * (pow(x(0) + 2, 2) + pow(x(1) + 2, 2)); // f2
    } else {
        matrix yt;
        yt = ff5T(ud2[0] + ud2[1]*x, ud1, NAN);//agregacja pierwszego wywolania ekspansji w celu umożliwienia poprawnego wykonania funkcji
        y = ud1(0) * yt(0) + (1 - ud1(0)) * yt(1);
    }
    return y;
}

matrix ff5R(matrix x, matrix ud1, matrix ud2) {
	matrix y;

	//Jeżeli ud2 puste
	if (isnan(ud2(0, 0))) {
		y = matrix(3, 1);

		//inicjalizacja danych z instrukcji z zamianą na SI
		double E = 207e9;
		double P = 1e3;
		double ro = 7800;

		//x(0) - l; x(1) - d
		//do y(0) obliczamy mase belki - pierwsze kryterium optymalizacji
		//wzor m = ro * V // V = pole podstawy * wysokość = pi * r^2 * l = pi * (d/2)^2 * l
		y(0) = ro * (3.14 * x(0) * (pow((x(1)), 2)))/4;

		//do y(1) obliczam ugięcie belki - drugie kryterium optymalziacji
		//wzór na ugięcie u z instrukcji
		y(1) = ((64 * P * pow(x(0), 3)) / (3 * E * 3.14 * pow(x(1), 4)));

		//do y(2) naprężenie
		//wzór na naprężenie sigma z instrukcji
		y(2) = ((32 * P * x(0)) / (3.14 * pow(x(1), 3)));
	}
	else {
		//agregacja pierwszego wywołania
		matrix ytmp, xtmp = ud2[0] + x * ud2[1];
		//cout << xtmp << endl;
		//matrix ytmp;

		ytmp = ff5R(xtmp, ud1, NAN);

		//metoda kryterium wazonego bez normalizacji
		//y = ud1 * (ytmp(0)) + (1 - ud1) * (ytmp(1));
		// 
		//metoda kryterium wazonego + normalizacja
		y = ud1 * (ytmp(0)-0.2)/(4.0-0.2) + (1 - ud1) * (ytmp(1)- 0.00005)/(0.005-0.00005);

		//element kary
		double c = 1e9;
		//sprawdz czy l nie jest za małe
		if (xtmp(0) < 0.2) {
			y = y + c * pow(xtmp(0) - 0.2, 2);
		}
		//sprawdz czy l nie jest za duze
		if (xtmp(0) > 1.0) {
			y = y + c * pow(xtmp(0) - 1.0, 2);
		}
		//sprawdz czy d nie jest za male
		if (xtmp(1) < 0.01) {
			y = y + c * pow(xtmp(1) - 0.01, 2);
		}
		//sprawdz czy d nie jest za duze
		if (xtmp(1) > 0.05) {
			y = y + c * pow(xtmp(1) - 0.05, 2);
		}
		//sprawdz czy ugięcie u nie przekracza maksimum
		if (ytmp(1) > 0.005) {
			y = y + c * pow(ytmp(1) - 0.005, 2);
		}
		//sprawdz czy napręzenie sigma nie przekracza maksimum
		if (ytmp(2) > 300e6) {
			y = y + c * pow(ytmp(2) - 300e6, 2);
		}
	}
	return y;
}


