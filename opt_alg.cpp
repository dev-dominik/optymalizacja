#include"opt_alg.h"

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	// Zmienne wej?ciowe:
	// ff - wska?nik do funkcji celu
	// N - liczba zmiennych funkcji celu
	// lb, ub - dolne i g?rne ograniczenie
	// epslion - zak??dana dok?adno?? rozwi?zania
	// Nmax - maksymalna liczba wywo?a? funkcji celu
	// ud1, ud2 - user data
	try
	{
		solution Xopt;
		Xopt = rand_mat(N);									// losujemy macierz Nx1 stosuj?c rozk?ad jednostajny na przedziale [0,1]
		for (int i = 0; i < N; ++i)							// przeskalowywujemy rozwi?zanie do przedzia?u [lb, ub]
			Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
		Xopt.fit_fun(ff, ud1, ud2);							// obliczamy warto?? funkcji celu
		Xopt.ud = (Xopt.y);									// warto?? funckji celu zapami?tujemy w macierzy ud
		if (Xopt.y < epsilon)								// sprawdzamy 1. kryterium stopu
		{
			Xopt.flag = 1;
			return Xopt;
		}
		solution X;
		while (true)
		{
			X = rand_mat(N);									// losujemy macierz Nx1 stosuj?c rozk?ad jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)							// przeskalowywujemy rozwi?zanie do przedzia?u [lb, ub]
				X.x(i) = (ub(i) - lb(i)) * X.x(i) + lb(i);
			X.fit_fun(ff, ud1, ud2);							// obliczmy warto?? funkcji celu
			if (X.y < Xopt.y)
			{
				Xopt = X;
				if (Xopt.y < epsilon)							// sprawdzmy 1. kryterium stopu
				{
					Xopt.flag = 1;								// flaga = 1 ozancza znalezienie rozwi?zanie z zadan? dok?adno?ci?
					break;
				}
			}
			if (solution::f_calls > Nmax)						// sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0;									// flaga = 0 ozancza przekroczenie maksymalnej liczby wywo?a? funkcji celu
				break;
			}
			Xopt.ud.add_row(Xopt.y);							// warto?? funckji dodajemy jako kolejny wiersz macierzy ud
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

		matrix historia;
		bool pierwsza = true;

		for (int i = 0; i <= k - 3; i++)
		{
			solution C(c), D(d);

			if (pierwsza)
			{
				historia = bi - ai;
				pierwsza = false;
			}
			else
				historia.add_row(bi - ai);

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

		historia.add_row(bi - ai);
		Xopt.ud = historia;
		
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
        solution Xopt;             // Obiekt przechowuj?cy wynik optymalizacji

        Xopt.ud = trans(x0);       // Zapisujemy punkt pocz?tkowy do historii punktów

        solution XB(x0), XB_old, X;
        // XB     - aktualny punkt bazowy (najlepszy znaleziony punkt)
        // XB_old - poprzedni punkt bazowy
        // X      - punkt próbny, który b?dziemy sprawdza?

        XB.fit_fun(ff, ud1, ud2);  // Obliczamy warto?? funkcji celu w punkcie pocz?tkowym

        while (true)              // G?ówna p?tla algorytmu
        {
            // ETAP 1: POSZUKIWANIE PRÓBNE
            // Sprawdzamy, czy mo?na poprawi? wynik, przesuwaj?c si?
            // o krok s wzd?u? poszczególnych osi uk?adu wspó?rz?dnych.
            // Robi to funkcja HJ_trial().

            X = HJ_trial(ff, XB, s, ud1, ud2);

            if (X.y < XB.y)       // Je?eli znaleziono punkt lepszy od bazowego
            {
                // ETAP 2: RUCH ROBOCZY
                // Skoro znale?li?my popraw?, próbujemy kontynuowa?
                // poszukiwanie w kierunku, w którym funkcja mala?a.

                while (true)
                {
                    XB_old = XB;  // Zapami?tujemy poprzedni punkt bazowy
                    XB = X;      // Nowy, lepszy punkt staje si? punktem bazowym

                    // Zapisujemy wspó?rz?dne nowego punktu bazowego do historii
                    // trans() zamienia wektor kolumnowy na wiersz,
                    // dzi?ki czemu ka?dy punkt zajmuje osobny wiersz macierzy.
                    Xopt.ud.add_row(trans(XB.x));

                    // Wyznaczamy nowy punkt w kierunku poprzedniej poprawy:
                    //
                    // X = XB + (XB - XB_old)
                    // X = 2*XB - XB_old
                    //
                    // Przyk?ad:
                    // XB_old = (1, 2)
                    // XB     = (2, 3)
                    // X      = (3, 4)
                    //
                    // Przesuwamy si? dalej w tym samym kierunku.

                    X.x = 2.0 * XB.x - XB_old.x;

                    X.fit_fun(ff, ud1, ud2);  // Obliczamy funkcj? celu w nowym punkcie

                    // Wykonujemy kolejne poszukiwanie próbne wokó? nowego punktu.
                    // Mo?e ono jeszcze bardziej poprawi? znalezione rozwi?zanie.
                    X = HJ_trial(ff, X, s, ud1, ud2);

                    // Je?eli nowy punkt nie jest lepszy od aktualnego bazowego,
                    // ko?czymy ruch roboczy i wracamy do g?ównej p?tli.
                    if (X.y >= XB.y)
                        break;

                    // Sprawdzamy, czy przekroczono dopuszczaln?
                    // liczb? wywo?a? funkcji celu.
                    if (solution::f_calls > Nmax)
                    {
                        Xopt = XB;       // Zwracamy najlepszy znaleziony punkt
                        Xopt.flag = 0;  // Flaga 0 oznacza przekroczenie limitu
                        return Xopt;    // Ko?czymy optymalizacj?
                    }
                }
            }
            else
            {
                // Je?eli poszukiwanie próbne nie przynios?o poprawy,
                // zmniejszamy d?ugo?? kroku.
                //
                // Przyk?ad:
                // s = 0.5, alpha = 0.5
                // nowy krok: s = 0.25
                //
                // Dzi?ki temu mo?emy dok?adniej przeszukiwa?
                // otoczenie aktualnego punktu bazowego.
                s *= alpha;
            }

            // PIERWSZE KRYTERIUM STOPU
            // Je?eli krok jest mniejszy od zadanej dok?adno?ci,
            // uznajemy, ?e znale?li?my rozwi?zanie.
            if (s < epsilon)
            {
                Xopt = XB;       // Zapisujemy najlepszy znaleziony punkt
                Xopt.flag = 1;   // Flaga 1 oznacza prawid?owe zako?czenie
                break;
            }

            // DRUGIE KRYTERIUM STOPU
            // Ko?czymy, je?eli wykorzystali?my limit wywo?a? funkcji celu.
            if (solution::f_calls > Nmax)
            {
                Xopt = XB;
                Xopt.flag = 0;  // Nie osi?gni?to kryterium dok?adno?ci
                break;
            }
        }

        return Xopt;             // Zwracamy wynik optymalizacji
    }
    catch (string ex_info)
    {
        // Przechwytujemy b??d i dodajemy informacj?,
        // ?e wyst?pi? on podczas dzia?ania funkcji HJ().
        throw ("solution HJ(...):\n" + ex_info);
    }
}

// Metoda Hooke'a-Jeevesa (HJ)
// ff      - wska?nik do funkcji celu
// x0      - punkt pocz?tkowy optymalizacji
// s       - pocz?tkowa d?ugo?? kroku poszukiwania
// alpha   - wspó?czynnik zmniejszania kroku (np. 0.5)
// epsilon - dok?adno??, przy której ko?czymy optymalizacj?
// Nmax    - maksymalna liczba wywo?a? funkcji celu
// ud1, ud2 - dodatkowe dane przekazywane do funkcji celu
solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
    try
    {
        int n = get_dim(XB);      // Liczba zmiennych funkcji celu

        matrix D = ident_mat(n);  // Macierz jednostkowa n x n
        // Jej kolumny okre?laj? podstawowe kierunki poszukiwania.
        //
        // Dla dwóch zmiennych:
        //
        //     [1  0]
        // D = [0  1]
        //
        // D[0] = (1, 0) - kierunek pierwszej wspó?rz?dnej
        // D[1] = (0, 1) - kierunek drugiej wspó?rz?dnej

        solution X;               // Punkt próbny

        // Sprawdzamy kolejno ka?dy kierunek poszukiwania.
        for (int i = 0; i < n; ++i)
        {
            // KROK W KIERUNKU DODATNIM
            //
            // Dla i = 0:
            // X = XB + s*(1,0)
            //
            // Dla i = 1:
            // X = XB + s*(0,1)

            X.x = XB.x + s * D[i];

            X.fit_fun(ff, ud1, ud2);  // Obliczamy warto?? funkcji w punkcie próbnym

            // Je?eli warto?? funkcji jest mniejsza,
            // przyjmujemy punkt próbny jako nowy punkt bazowy.
            if (X.y < XB.y)
                XB = X;
            else
            {
                // KROK W KIERUNKU UJEMNYM
                //
                // Je?eli ruch dodatni nie poprawi? wyniku,
                // sprawdzamy przeciwny kierunek.

                X.x = XB.x - s * D[i];

                X.fit_fun(ff, ud1, ud2);

                // Je?eli ruch ujemny poprawi? wynik,
                // aktualizujemy punkt bazowy.
                if (X.y < XB.y)
                    XB = X;
            }

            // Po sprawdzeniu jednego kierunku przechodzimy
            // do kolejnej zmiennej.
            //
            // Wa?ne: je?eli znale?li?my popraw?,
            // kolejne kroki wykonujemy ju? od nowego XB.
        }

        // Zwracamy najlepszy punkt znaleziony
        // podczas poszukiwania próbnego.
        return XB;
    }
    catch (string ex_info)
    {
        throw ("solution HJ_trial(...):\n" + ex_info);
    }
}

// Funkcja pomocnicza metody Hooke'a-Jeevesa
// HJ_trial() realizuje POSZUKIWANIE PRÓBNE.
//
// XB - punkt bazowy, od którego rozpoczynamy poszukiwanie
// s  - d?ugo?? kroku
//
// Dla ka?dej zmiennej sprawdzamy najpierw kierunek dodatni,
// a je?eli nie ma poprawy, kierunek ujemny.
//
// Dla funkcji dwóch zmiennych sprawdzamy kolejno:
// (x1+s, x2), (x1-s, x2), (x1, x2+s), (x1, x2-s)
//
// Nie zawsze sprawdzane s? wszystkie cztery punkty,
// poniewa? po znalezieniu poprawy nie sprawdzamy
// przeciwnego kierunku dla tej samej wspó?rz?dnej.
// Metoda Rosenbrocka (Rosen)
// Jest to metoda bezgradientowa, która oprócz zmiany
// d?ugo?ci kroków potrafi równie? zmienia? kierunki poszukiwania.
//
// W metodzie Hooke'a-Jeevesa poszukujemy wzd?u? sta?ych osi.
// W metodzie Rosenbrocka kierunki mog? zosta? obrócone,
// aby lepiej odpowiada?y kszta?towi optymalizowanej funkcji.
//
// ff      - funkcja celu
// x0      - punkt pocz?tkowy
// s0      - wektor pocz?tkowych d?ugo?ci kroków
// alpha   - wspó?czynnik zwi?kszania kroku po poprawie (np. 2)
// beta    - wspó?czynnik zmniejszania kroku po niepowodzeniu (np. 0.5)
// epsilon - zadana dok?adno??
// Nmax    - maksymalna liczba wywo?a? funkcji celu
// ud1, ud2 - dodatkowe dane funkcji celu

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
    try
    {
        solution Xopt;           // Obiekt przechowuj?cy wynik

        Xopt.ud = trans(x0);     // Zapisujemy punkt pocz?tkowy do historii

        solution XB(x0), X;
        // XB - aktualny najlepszy punkt
        // X  - punkt próbny

        int n = get_dim(XB);     // Liczba zmiennych funkcji celu

        matrix l(n, 1), p(n, 1), s(s0), D = ident_mat(n);
        // l - suma udanych przesuni?? w poszczególnych kierunkach
        // p - licznik nieudanych prób w poszczególnych kierunkach
        // s - aktualne d?ugo?ci kroków (osobny krok dla ka?dej zmiennej)
        // D - macierz kierunków poszukiwania
        //
        // Pocz?tkowo D jest macierz? jednostkow?,
        // czyli poszukujemy wzd?u? osi uk?adu wspó?rz?dnych.
        //
        // W odró?nieniu od HJ, gdzie mamy jeden krok s,
        // tutaj ka?dy kierunek ma w?asny krok s(i).

        XB.fit_fun(ff, ud1, ud2); // Obliczamy funkcj? celu w punkcie pocz?tkowym

        while (true)             // G?ówna p?tla metody Rosenbrocka
        {
            // ETAP 1: POSZUKIWANIE W AKTUALNYCH KIERUNKACH
            //
            // Sprawdzamy kolejno ka?dy z n kierunków.

            for (int i = 0; i < n; ++i)
            {
                // Wyznaczamy nowy punkt:
                //
                // X = XB + s(i)*D[i]
                //
                // D[i] okre?la kierunek,
                // a s(i) okre?la d?ugo?? i zwrot kroku.
                //
                // Ujemne s(i) oznacza ruch w przeciwn? stron?.

                X.x = XB.x + s(i) * D[i];

                X.fit_fun(ff, ud1, ud2); // Obliczamy funkcj? celu

                if (X.y < XB.y)
                {
                    // SUKCES - znaleziono lepszy punkt

                    XB = X;          // Aktualizujemy najlepszy punkt

                    l(i) += s(i);    // Sumujemy udane przesuni?cia
                    // w kierunku i

                    s(i) *= alpha;   // Zwi?kszamy d?ugo?? kroku
                    // w tym kierunku

                    // Przyk?ad:
                    // s(i) = 0.5, alpha = 2
                    // nowy krok: s(i) = 1.0
                    //
                    // Skoro ruch poprawi? wynik,
                    // próbujemy wykonywa? wi?ksze kroki.
                }
                else
                {
                    // NIEPOWODZENIE - nie znaleziono lepszego punktu

                    ++p(i);         // Zwi?kszamy licznik niepowodze?

                    s(i) *= -beta;   // Zmniejszamy d?ugo?? kroku
                    // i zmieniamy jego znak

                    // Przyk?ad:
                    // s(i) = 1.0, beta = 0.5
                    // nowy krok: s(i) = -0.5
                    //
                    // W kolejnej próbie sprawdzimy przeciwny kierunek,
                    // wykonuj?c jednocze?nie krótszy krok.
                }
            }

            // Zapisujemy aktualny najlepszy punkt do historii.
            // Dla funkcji dwóch zmiennych ka?dy wiersz
            // zawiera wspó?rz?dne (x1, x2).
            Xopt.ud.add_row(trans(XB.x));

            // ETAP 2: SPRAWDZENIE MO?LIWO?CI ZMIANY KIERUNKÓW
            //
            // Nowe kierunki wyznaczamy dopiero wtedy,
            // gdy w KA?DYM kierunku wyst?pi?:
            //
            // - udane przesuni?cia (l(i) != 0),
            // - nieudane próby (p(i) != 0).

            bool change = true;

            for (int i = 0; i < n; ++i)
                if (l(i) == 0 || p(i) == 0)
                {
                    // Je?eli w którym? kierunku nie spe?niono
                    // obu warunków, nie zmieniamy jeszcze kierunków.
                    change = false;
                    break;
                }

            if (change)
            {
                // ETAP 3: WYZNACZANIE NOWYCH KIERUNKÓW
                //
                // Budujemy nowy uk?ad kierunków poszukiwania.
                // Wykorzystujemy informacje o udanych przesuni?ciach,
                // zapisane w wektorze l.

                matrix Q(n, n), v(n, 1);

                // Tworzymy pomocnicz? macierz Q.
                // Wype?niamy jej doln? cz??? warto?ciami l(i).
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j <= i; ++j)
                        Q(i, j) = l(i);

                // Przekszta?camy kierunki za pomoc? aktualnej macierzy D.
                // Dzi?ki temu nowe kierunki uwzgl?dniaj? dotychczasowy
                // uk?ad wspó?rz?dnych poszukiwania.
                Q = D * Q;

                // Pierwszy nowy kierunek otrzymujemy przez
                // normalizacj? pierwszego wektora Q.
                //
                // norm() oblicza d?ugo?? wektora.
                // Dzielenie przez d?ugo?? daje wektor jednostkowy.

                v = Q[0] / norm(Q[0]);

                // Zapisujemy pierwszy nowy kierunek w macierzy D.
                D.set_col(v, 0);

                // Kolejne kierunki wyznaczamy metod?
                // ortogonalizacji Grama-Schmidta.
                //
                // Jej zadaniem jest utworzenie wektorów:
                // - prostopad?ych do siebie,
                // - o d?ugo?ci równej 1.
                //
                // Dzi?ki temu otrzymujemy nowy uk?ad kierunków
                // dostosowany do przebiegu dotychczasowych poszukiwa?.

                for (int i = 1; i < n; ++i)
                {
                    matrix temp(n, 1);

                    // Obliczamy sum? rzutów aktualnego wektora Q[i]
                    // na wcze?niej wyznaczone kierunki D[j].
                    //
                    // Odejmuj?c t? sum?, usuwamy sk?adowe
                    // równoleg?e do poprzednich kierunków.

                    for (int j = 0; j < i; ++j)
                        temp = temp + trans(Q[i]) * D[j] * D[j];

                    // Wyznaczamy sk?adow? prostopad?? do poprzednich kierunków
                    // i normalizujemy j? do d?ugo?ci 1.

                    v = (Q[i] - temp) / norm(Q[i] - temp);

                    // Zapisujemy nowy kierunek w odpowiedniej kolumnie D.
                    D.set_col(v, i);
                }

                // Po zmianie kierunków rozpoczynamy kolejny etap poszukiwa?.
                //
                // Przywracamy pocz?tkowe d?ugo?ci kroków.
                s = s0;

                // Zerujemy informacje o udanych przesuni?ciach.
                l = matrix(n, 1);

                // Zerujemy liczniki niepowodze?.
                p = matrix(n, 1);
            }

            // ETAP 4: SPRAWDZENIE KRYTERIÓW STOPU

            // Szukamy najwi?kszej bezwzgl?dnej d?ugo?ci kroku.
            //
            // U?ywamy abs(), poniewa? krok mo?e by? ujemny.
            // Interesuje nas jego d?ugo??, a nie zwrot.

            double max_s = abs(s(0));

            for (int i = 1; i < n; ++i)
                if (max_s < abs(s(i)))
                    max_s = abs(s(i));

            // PIERWSZE KRYTERIUM STOPU
            // Je?eli wszystkie kroki s? mniejsze od epsilon,
            // ko?czymy optymalizacj?.
            //
            // Wystarczy sprawdzi? najwi?kszy krok.
            // Je?eli najwi?kszy jest mniejszy od epsilon,
            // to pozosta?e równie?.

            if (max_s < epsilon)
            {
                Xopt = XB;       // Zwracamy najlepszy znaleziony punkt
                Xopt.flag = 1;   // Flaga 1 - prawid?owe zako?czenie
                break;
            }

            // DRUGIE KRYTERIUM STOPU
            // Ko?czymy, gdy przekroczymy limit
            // wywo?a? funkcji celu.

            if (solution::f_calls > Nmax)
            {
                Xopt = XB;
                Xopt.flag = 0;   // Flaga 0 - przekroczony limit
                break;
            }
        }

        return Xopt;             // Zwracamy wynik optymalizacji
    }
    catch (string ex_info)
    {
        // Obs?uga b??dów wyst?puj?cych podczas dzia?ania Rosen().
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
