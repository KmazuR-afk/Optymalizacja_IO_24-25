#include"opt_alg.h"
#include <cmath>

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
			Xopt.fit_fun(ff, ud1, ud2);
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2] { 0, 0 };
		//Tu wpisz kod funkcji
		solution X0(x0), X1(x0 + d);
		X0.fit_fun(ff, ud1, ud2);
		X1.fit_fun(ff, ud1, ud2);
		if (X0.y == X1.y) {
			p[0] = m2d(X0.x);
			p[1] = m2d(X1.x);
			return p;
		}

		if (X1.y > X0.y) {
			d = -d;
			X1.x = X0.x + d;
			X1.fit_fun(ff, ud1, ud2);
			if (X1.y >= X0.y) {
				p[0] = m2d(X1.x);
				p[1] = m2d(X0.x - d);
				return p;
			}
		}
		int i = 0;
		solution Xi = X0;
		while (true) {
			if (X0.f_calls > Nmax) {
				return NULL;
			}
			i = i + 1;
			if (i > 1) {
				Xi = X1;
			}
			X1.x = X0.x + pow(alpha, i) * d;
			X1.fit_fun(ff, ud1, ud2);
			Xi.fit_fun(ff, ud1, ud2);
			if (Xi.y <= X1.y) {
				break;
			}
		}
		solution Xim;
		Xim.x = X0.x + pow(alpha, i - 2) * d;
		if (d > 0) {
			p[0] = m2d(Xim.x);
			p[1] = m2d(X1.x);
			return p;
		}
		p[0] = m2d(X1.x);
		p[1] = m2d(Xim.x);
		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

double oblicz_fib(int k) {
	double x = 0;
	double y = 1;

	for (int i = 0; i < k; i++) {
		y = y + x;
		x = y - x;
	}
	return y;
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		solution a0(a), b0(b), c0, d0;
		int k = 0;
		while (oblicz_fib(k) <= (b - a) / epsilon) {
			k++;
		}
		c0.x = b0.x - ((oblicz_fib(k - 1) / (oblicz_fib(k))) * (b0.x - a0.x));
		d0.x = a0.x + b0.x - c0.x;
		for (int i = 0; i < k - 3; i++) {
			c0.fit_fun(ff, ud1, ud2);
			d0.fit_fun(ff, ud1, ud2);
			if (c0.y < d0.y) {
				a0.x = a0.x;
				b0.x = d0.x;
			}
			else {
				b0.x = b0.x;
				a0.x = c0.x;
			}
			c0.x = b0.x - ((oblicz_fib(k - i - 2) / (oblicz_fib(k - i - 1))) * (b0.x - a0.x));
			d0.x = a0.x + b0.x - c0.x;

		}
		Xopt = c0;
		Xopt.fit_fun(ff, ud1, ud2);
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		int i = 0;
		solution Xa(a), Xb(b), Xc((a + b) / 2), Xd, Xd_prev;
		while (true) {
			Xa.fit_fun(ff, ud1, ud2);
			Xb.fit_fun(ff, ud1, ud2);
			Xc.fit_fun(ff, ud1, ud2);
			matrix l = (Xa.y * (pow(Xb.x, 2) - pow(Xc.x, 2)) + Xb.y * (pow(Xc.x, 2) - pow(Xa.x, 2)) + Xc.y * (pow(Xa.x, 2) - pow(Xb.x, 2)));
			matrix m = (Xa.y * (Xb.x - Xc.x) + Xb.y * (Xc.x - Xa.x) + Xc.y * (Xa.x - Xb.x));
			if (m2d(m) <= 0) {
				cout << i;
				return NULL;
			}
			Xd_prev.x = Xd.x;
			Xd.x = 0.5 * l / m;

			Xd_prev.fit_fun(ff, ud1, ud2);
			Xd.fit_fun(ff, ud1, ud2);
			if (Xa.x < Xd.x && Xd.x < Xc.x) {
				if (Xd.y < Xc.y) {
					Xa.x = Xa.x;
					Xb.x = Xc.x;
					Xc.x = Xd.x;
				}
				else {
					Xa.x = Xd.x;
					Xc.x = Xc.x;
					Xb.x = Xb.x;
				}
			}

			else if (Xd.x == Xc.x) {
				if (fabs(m2d(Xd.x - Xc.x)) < gamma) {
					break;
				}
			}
			else {
				if (Xc.x < Xd.x && Xd.x < Xb.x) {
					Xd.fit_fun(ff, ud1, ud2);
					if (Xd.y < Xc.y) {
						Xa.x = Xc.x;
						Xc.x = Xd.x;
						Xb.x = Xb.x;
					}
					else {
						Xa.x = Xa.x;
						Xc.x = Xc.x;
						Xb.x = Xd.x;
					}
				}
				else {
					return NULL;
				}
			}
			i = i + 1;
			if (Xopt.f_calls > Nmax) {
				return NULL;
			}
			if (m2d(Xb.x - Xa.x) < epsilon || fabs(m2d(Xd.x - Xd_prev.x)) < gamma) {
				break;
			}
		}
		Xopt.x = Xd.x;
		Xopt.fit_fun(ff, ud1, ud2);
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		solution XB(x0);
		Xopt=HJ_trial(ff,x0,s);
		Xopt.fit_fun(ff);
		XB.fit_fun(ff);
		while(true){
			if(Xopt.y<XB.y){
				solution X_B;
				while(true){
					X_B=XB;
					XB=Xopt;
					Xopt.x=2*XB.x-X_B.x;
					Xopt=HJ_trial(ff,Xopt,s);
					XB.fit_fun(ff);
					Xopt.fit_fun(ff);
					if(solution::f_calls>Nmax){
						cout<<"przekroczono maksymalną ilość wywołania funkcji\n";
						Xopt.flag=0;
						return NULL;
					}
					else if (Xopt.y>=XB.y)
					{
						break;
					}	
				}
				Xopt=XB;
			}
			else{
				s*=alpha;
			}
			if(solution::f_calls>Nmax){
				cout<<"przekroczono maksymalną ilość wywołania funkcji\n";
				Xopt.flag=0;
				return NULL;
			}
			if(s<epsilon){
				Xopt.flag=1;
				break;
			}
		}
		Xopt.fit_fun(ff);
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		int dim=get_dim(XB);
		solution X;
		matrix e=ident_mat(dim);
		XB.fit_fun(ff);
		for(int j=0;j<dim-1;j++){
			X.x=XB.x+s*e[j];
			X.fit_fun(ff);
			if(X.y<XB.y){
				XB=X;
			}
			else{
				X.x=XB.x-2*s*e[j];
				X.fit_fun(ff);
				if(X.y<XB.y){
					XB=X;
				}
			}
		}


		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
