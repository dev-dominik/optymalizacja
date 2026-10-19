#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	// Zmienne wejœciowe:
	// ff - wskaŸnik do funkcji celu
	// N - liczba zmiennych funkcji celu
	// lb, ub - dolne i górne ograniczenie
	// epslion - zak³¹dana dok³adnoœæ rozwi¹zania
	// Nmax - maksymalna liczba wywo³añ funkcji celu
	// ud1, ud2 - user data
	try
	{
		solution Xopt;
		Xopt = rand_mat(N);									// losujemy macierz Nx1 stosuj¹c rozk³ad jednostajny na przedziale [0,1]
		for (int i = 0; i < N; ++i)							// przeskalowywujemy rozwi¹zanie do przedzia³u [lb, ub]
			Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
		Xopt.fit_fun(ff, ud1, ud2);							// obliczamy wartoœæ funkcji celu
		Xopt.ud = (Xopt.y);									// wartoœæ funckji celu zapamiêtujemy w macierzy ud
		if (Xopt.y < epsilon)								// sprawdzamy 1. kryterium stopu
		{
			Xopt.flag = 1;
			return Xopt;
		}
		solution X;
		while (true)
		{
			X = rand_mat(N);									// losujemy macierz Nx1 stosuj¹c rozk³ad jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)							// przeskalowywujemy rozwi¹zanie do przedzia³u [lb, ub]
				X.x(i) = (ub(i) - lb(i)) * X.x(i) + lb(i);
			X.fit_fun(ff, ud1, ud2);							// obliczmy wartoœæ funkcji celu
			if (X.y < Xopt.y)
			{
				Xopt = X;
				if (Xopt.y < epsilon)							// sprawdzmy 1. kryterium stopu
				{
					Xopt.flag = 1;								// flaga = 1 ozancza znalezienie rozwi¹zanie z zadan¹ dok³adnoœci¹
					break;
				}
			}
			if (solution::f_calls > Nmax)						// sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0;									// flaga = 0 ozancza przekroczenie maksymalnej liczby wywo³añ funkcji celu
				break;
			}
			Xopt.ud.add_row(Xopt.y);							// wartoœæ funckji dodajemy jako kolejny wiersz macierzy ud
		}
		Xopt.ud.add_row(Xopt.y);
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
		//kod
		if (d <= 0 || alpha <= 1 || Nmax <= 0)
			throw string("Niepoprawne parametry ekspansji");
		//koniec
		double* p = new double[2] { 0, 0 };
		//kod
		
		auto f = [&](double x) { solution X(x); X.fit_fun(ff, ud1, ud2); return m2d(X.y); };
		double x1 = x0 + d;
		double y0 = f(x0), y1 = f(x1);
		if (y1 == y0) { p[0] = x0; p[1] = x1; return p; }
		if (y1 > y0)
		{
			d = -d; x1 = x0 + d; y1 = f(x1);
			if (y1 >= y0)
			{
				p[0] = min(x1, x0 - d); p[1] = max(x1, x0 - d);
				return p;
			}
		}
		double x_prev = x0, x_curr = x1, y_curr = y1;
		for (int i = 1; solution::f_calls < Nmax; ++i)
		{
			double x_next = x0 + pow(alpha, i) * d;
			double y_next = f(x_next);
			if (y_next >= y_curr)
			{
				p[0] = min(x_prev, x_next);
				p[1] = max(x_prev, x_next);
				return p;
			}
			x_prev = x_curr; x_curr = x_next; y_curr = y_next;
		}
		delete[] p;
		throw string("Przekroczono Nmax w metodzie ekspansji");
		
		//return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//kod
		if (epsilon <= 0 || a >= b)
			throw string("Niepoprawny przedzial lub epsilon");

		vector<double> F = { 0, 1 };

		while (F.back() <= (b - a) / epsilon)
		{
			int n = (int)F.size();
			F.push_back(F[n - 1] + F[n - 2]);
		}

		int k = (int)F.size() - 1;

		double ai = a, bi = b;

		double c = bi - F[k - 1] / F[k] * (bi - ai);
		double d = ai + bi - c;

		for (int i = 0; i <= k - 3; i++)
		{
			solution C(c), D(d);

			C.fit_fun(ff, ud1, ud2);
			D.fit_fun(ff, ud1, ud2);

			if (C.y < D.y)
			{
				bi = d;
			}
			else
			{
				ai = c;
			}

			c = bi - F[k - i - 2] / F[k - i - 1] * (bi - ai);
			d = ai + bi - c;
		}

		Xopt.x = c;
		Xopt.fit_fun(ff, ud1, ud2);
		Xopt.flag = 1;
		//koniec 
		
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
		//Tu wpisz kod funkcji
		if (a >= b || epsilon <= 0 || gamma <= 0 || Nmax <= 0)
			throw string("Niepoprawne parametry metody Lagrange'a");
		double c = (a + b) / 2, d = c, d_prev = c;
		matrix historia; bool pierwsza = true;
		for (int i = 0; ; ++i)
		{
			if (pierwsza) { historia = b - a; pierwsza = false; }
			else historia.add_row(b - a);
			if (solution::f_calls + 4 > Nmax) { Xopt.flag = 0; break; }
			solution A(a), B(b), C(c);
			A.fit_fun(ff, ud1, ud2); B.fit_fun(ff, ud1, ud2); C.fit_fun(ff, ud1, ud2);
			double fa = m2d(A.y), fb = m2d(B.y), fc = m2d(C.y);
			double l = fa * (b * b - c * c) + fb * (c * c - a * a) + fc * (a * a - b * b);
			double m = fa * (b - c) + fb * (c - a) + fc * (a - b);
			if (m <= 0 || !isfinite(m)) { Xopt.flag = -1; break; }
			d = 0.5 * l / m;
			if (!isfinite(d) || d <= a || d >= b || d == c) { Xopt.flag = -1; break; }
			solution D(d); D.fit_fun(ff, ud1, ud2);
			if (d < c)
			{
				if (D.y < C.y) { b = c; c = d; }
				else a = d;
			}
			else
			{
				if (D.y < C.y) { a = c; c = d; }
				else b = d;
			}
			if (b - a < epsilon || (i > 0 && abs(d - d_prev) < gamma))
			{
				Xopt.flag = 1; break;
			}
			d_prev = d;
		}
		historia.add_row(b - a);
		Xopt.x = (Xopt.flag == 1) ? d : c;
		Xopt.fit_fun(ff, ud1, ud2);
		Xopt.ud = historia;
		//koniec

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
		//Tu wpisz kod funkcji

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
		//Tu wpisz kod funkcji

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
