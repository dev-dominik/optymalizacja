/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 15.09.2026
*********************************************/

#include"opt_alg.h"

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
		lab1();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	return 0;
}

void lab0()
{
	double epsilon;												// dok³adnoœæ
	int Nmax;													// maksymalna liczba wywo³añ funkcji celu
	matrix lb, ub;												// dolne oraz górne ograniczenie
	solution opt;												// rozwi¹zanie optymalne znalezione przez algorytm
	if (0) // Funkcja testowa
	{
		epsilon = 1e-2;
		Nmax = 10000;
		lb = matrix(2, 1, -5);
		ub = matrix(2, 1, 5);
		matrix a(2, 1);											// dok³adne rozwi¹zanie optymalne
		a(0) = -1;
		a(1) = 2;
		opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);			// wywo³anie procedury optymalizacji
		cout << opt << endl << endl;							// wypisanie wyniku
		solution::clear_calls();								// wyzerowanie liczników
		int* n = get_size(opt.ud);
		cout << n[0] << '\t' << n[1] << endl << endl;
	}
	if (1) // Wahad³o
	{
		Nmax = 1000;
		epsilon = 1e-2;
		lb = 0, ub = 5;
		double teta_opt = 1;									// maksymalne wychylenie wahad³a
		opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);		// wywo³anie procedury optymalizacji
		cout << opt << endl << endl;							// wypisanie wyniku
		solution::clear_calls();								// wyzerowanie liczników

		// Zapis zmian wartoœci funckji celu do pliku
		ofstream Sout("optymalicja_lab0.csv");					// definiujemy strumieñ do pliku .csv
		Sout << opt.ud;											// zapisujemy wyniki w pliku
		Sout.close();											// zamykamy strumieñ

		// Zapis symulacji do pliku csv
		matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz¹tkowe
			MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment si³y dzia³aj¹cy na wahad³o oraz czas dzia³ania
		matrix* Y = solve_ode(df0, 0, 0.1, 20, Y0, NAN, MT);	// rozwi¹zujemy równanie ró¿niczkowe
		Sout.open("symulacja_lab0.csv");						// otwieramy plik .csv
		Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
		Sout.close();											// zamykamy strumieñ
		Y[0].~matrix();											// usuwamy z pamiêci rozwi¹zanie RR
		Y[1].~matrix();
	}
}

void lab1()
{
    double epsilon = 1e-2, gamma = 1e-4;
    int Nmax = 10000;
    double alpha[3] = { 1.5, 2.0, 3.0 };
    mt19937 gen(random_device{}());
    uniform_real_distribution<double> losuj(-100, 100);
    ofstream Sout("wyniki_lab1.csv");
    Sout << "proba;alpha;x0;metoda;a;b;x;y;f_calls;flag\n";
    for (int j = 0; j < 3; ++j)
    {
        for (int i = 0; i < 100; ++i)
        {
            double x0 = losuj(gen);
            solution::clear_calls();
            double* przedzial = expansion(ff1T, x0, 1, alpha[j], Nmax);
            double a = przedzial[0], b = przedzial[1];
            delete[] przedzial;
            solution::clear_calls();
            solution opt = fib(ff1T, a, b, epsilon);
            Sout << i + 1 << ';' << alpha[j] << ';' << x0 << ";fib;" << a << ';' << b << ';'
                << m2d(opt.x) << ';' << m2d(opt.y) << ';' << solution::f_calls << ';' << opt.flag << '\n';
            solution::clear_calls();
            opt = lag(ff1T, a, b, epsilon, gamma, Nmax);
            Sout << i + 1 << ';' << alpha[j] << ';' << x0 << ";lag;" << a << ';' << b << ';'
                << m2d(opt.x) << ';' << m2d(opt.y) << ';' << solution::f_calls << ';' << opt.flag << '\n';
        }
    }
    solution::clear_calls();
    solution opt = fib(ff1T, -100, 100, epsilon);
    Sout << "0;0;0;fib_bez_ekspansji;-100;100;" << m2d(opt.x) << ';' << m2d(opt.y)
        << ';' << solution::f_calls << ';' << opt.flag << '\n';
    ofstream Historia("przedzial_fib_lab1.csv");
    for (int i = 0; i < get_len(opt.ud); ++i) Historia << i << ';' << opt.ud(i) << '\n';
    Historia.close();
    solution::clear_calls();
    opt = lag(ff1T, -100, 100, epsilon, gamma, Nmax);
    Sout << "0;0;0;lag_bez_ekspansji;-100;100;" << m2d(opt.x) << ';' << m2d(opt.y)
        << ';' << solution::f_calls << ';' << opt.flag << '\n';
    Historia.open("przedzial_lag_lab1.csv");
    for (int i = 0; i < get_len(opt.ud); ++i) Historia << i << ';' << opt.ud(i) << '\n';
    Historia.close();
    Sout.close();

    ofstream Real("wyniki_zbiorniki_lab1.csv");
    Real << "metoda;DA_cm2;f_celu;f_calls;flag\n";
    solution::clear_calls();
    solution opt_fib = fib(ff1R, 1, 100, epsilon);
    Real << "fib;" << m2d(opt_fib.x) << ';' << m2d(opt_fib.y) << ';' << solution::f_calls << ';' << opt_fib.flag << '\n';
    solution::clear_calls();
    solution opt_lag = lag(ff1R, 1, 100, epsilon, gamma, Nmax);
    Real << "lag;" << m2d(opt_lag.x) << ';' << m2d(opt_lag.y) << ';' << solution::f_calls << ';' << opt_lag.flag << '\n';
    Real.close();
    matrix Y0(3, new double[3] { 5, 1, 20 });
    matrix* Y = solve_ode(df1, 0, 1, 2000, Y0, NAN, opt_fib.x);
    ofstream Sym("symulacja_lab1.csv");
    Sym << "czas_s;VA_m3;VB_m3;TB_C\n";
    for (int i = 0; i < get_len(Y[0]); ++i)
        Sym << Y[0](i) << ';' << Y[1](i, 0) << ';' << Y[1](i, 1) << ';' << Y[1](i, 2) << '\n';
    Sym.close();
    delete[] Y;
    cout << "Lab1 zakonczone. Wyniki zapisano do plikow CSV." << endl;
    cout << "Zbiorniki fib: " << opt_fib << endl;
    cout << "Zbiorniki lag: " << opt_lag << endl;



}

void lab2()
{

}

void lab3()
{

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
