#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <iomanip> 

using namespace std;

int liczba = 3; 

double min(double tablica[], int nN)
{
	double tmp = tablica[0];
	for (int i = 1; i < nN; i++)
	{
		if (tmp > tablica[i])
		{
			tmp = tablica[i];
		}
	}
	return tmp;
}

double max(double tablica[], int nN)
{
	double tmp = tablica[0];
	for (int i = 1; i < nN; i++)
	{
		if (tmp < tablica[i])
		{
			tmp = tablica[i];
		}
	}
	return tmp;
}

struct fun_ksztaltu
{
	double** tabela_dnPOksi = new double* [liczba * liczba];
	double** tabela_dnPOeta = new double* [liczba * liczba];
	double** fun_ksztaltu_dlaC = new double* [liczba * liczba];

	double** punkty_calkowania = new double* [4 * liczba];
	double** tablica_funkcja_ksztaltu = new double* [liczba * 4];
};

struct Elem
{
	double** tabela_dnPOdx = new double* [liczba * liczba];
	double** tabela_dnPOdy = new double* [liczba * liczba];

	double*** h = new double** [liczba * liczba];

	double h_koncowe[4][4];

	double hbc[4][4];

	double** jakobian = new double* [liczba * liczba];

	double* wyznacznik = new double[liczba * liczba];

	double l[4];

	double P[4] = { NULL };

	double c[4][4];
};

struct node {
	double x;
	double y;
	int BC;
};

struct element {
	int ID[4];
};

struct grid {
	int nN; //liczba wezlow
	int nE; //liczba elementow
	node* ND; //tablica z nodami
	element* EL; //tablica z elementami
};

struct globaldata { //reszta danych
	int SimulationTime;
	int SimulationStepTime;
	int Conductivity; //wspolczynnik przewodnosci cieplnej
	int Alfa;
	int Tot;
	int InitialTemp;
	int Density;
	int SpecificHeat;
};

struct dane {
	double* P_agr;
	double* PplusCdT;
	double* t0;
};

struct uklad_rownan {
	double** hg;
	double** h;
	double** C_agregowane;
	double** HplusCdT;
};

double n1odKsi(double eta)
{
	return -0.25 * (1 - eta);
}

double n2odKsi(double eta)
{
	return 0.25 * (1 - eta);
}

double n3odKsi(double eta)
{
	return 0.25 * (1 + eta);
}

double n4odKsi(double eta)
{
	return -0.25 * (1 + eta);
}




double n1odEta(double ksi)
{
	return -0.25 * (1 - ksi);
}

double n2odEta(double ksi)
{
	return -0.25 * (1 + ksi);
}

double n3odEta(double ksi)
{
	return 0.25 * (1 + ksi);
}

double n4odEta(double ksi)
{
	return 0.25 * (1 - ksi);
}

double n1(double ksi, double eta)
{
	return 0.25 * (1 - ksi) * (1 - eta);
}

double n2(double ksi, double eta)
{
	return 0.25 * (1 + ksi) * (1 - eta);
}

double n3(double ksi, double eta)
{
	return 0.25 * (1 + ksi) * (1 + eta);
}

double n4(double ksi, double eta)
{
	return 0.25 * (1 - ksi) * (1 + eta);
}


void wczytajDoPliku(grid* siatka, globaldata* data)
{
	string smietnik, x, y, ele1, ele2, ele3, ele4;

	ifstream plik;
	//plik.open("Test1_4_4.txt");
	//plik.open("Test2_4_4_MixGrid.txt");
	plik.open("Test3_31_31_kwadrat.txt");

	plik >> smietnik >> data->SimulationTime;
	plik >> smietnik >> data->SimulationStepTime;
	plik >> smietnik >> data->Conductivity; //wspolczynnik przewodnosci cieplnej
	plik >> smietnik >> data->Alfa;
	plik >> smietnik >> data->Tot;
	plik >> smietnik >> data->InitialTemp;
	plik >> smietnik >> data->Density;
	plik >> smietnik >> data->SpecificHeat;
	plik >> smietnik >> smietnik >> siatka->nN;
	plik >> smietnik >> smietnik >> siatka->nE; //wyczytanie liczby elementow

	//utworzenie tablicy dynamicznej a potem doczepienie jej do obiektu struktury
	siatka->ND = new node[siatka->nN];

	//utworzenie tab dynamicznej a potem doczepienie jej do EL w strukturze grid
	siatka->EL = new element[siatka->nE];

	plik >> smietnik;

	//cout << "Nody:" << endl;
	for (int i = 0; i < siatka->nN; i++) //wprowadzenie nodow z pliku do struktur
	{
		plik >> smietnik >> x >> y; //wyczytuje jedna linijke
		smietnik.pop_back(); //usuwa przecinki
		x.pop_back(); //usuwa przecinki
		siatka->ND[i].x = stod(x); //wprowadza wartosci do struktur
		siatka->ND[i].y = stod(y);
		//cout << siatka.ND[i].x << " " << siatka.ND[i].y << endl; //wypisanie wsystkich nodow
	}

	plik >> smietnik >> smietnik;

	//cout << "Elementy:" << endl;
	for (int i = 0; i < siatka->nE; i++)
	{
		plik >> smietnik >> ele1 >> ele2 >> ele3 >> ele4;
		ele1.pop_back(); //usuwanie przecinkow
		ele2.pop_back();
		ele3.pop_back();

		siatka->EL[i].ID[0] = stoi(ele1); //wprawadzanie wartosci do 4 elementowej tablicy
		siatka->EL[i].ID[1] = stoi(ele2);
		siatka->EL[i].ID[2] = stoi(ele3);
		siatka->EL[i].ID[3] = stoi(ele4);

		//cout << ele1 << " " << ele2 << " " << ele3 << " " << ele4 << endl; //wypisanie elementow
	}

	plik >> smietnik; // czyta "*BC"

	while (1)
	{
		plik >> smietnik; //wczytuje po jednej wartosci plus potencjalny przecinek
		if (smietnik.back() == ',')
		{
			smietnik.pop_back();
			siatka->ND[stoi(smietnik) - 1].BC = 1;
		}
		else //ostatnia wartosc wykrywa sie po tym ze nie ma przecinka wiec wtedy petla sie zatrzymuje
		{
			siatka->ND[stoi(smietnik) - 1].BC = 1;
			break;
		}
	}

	//cout << "BC:" << endl;
	for (int i = 0; i < siatka->nN; i++)
	{
		if (siatka->ND[i].BC != 1) //wypełnianie pozostalosci zerami
		{							//kiedy 2 wejdzie do gry to doda sie kolejny warunek ||
			siatka->ND[i].BC = 0;
		}
		//cout << siatka->ND[i].BC << endl; //wypisanie BC
	}
}

void przydzielWagiiPkt(double* punkty, double* wagi)
{
	if (liczba == 2)
	{
		punkty[0] = -1. / sqrt(3.);
		punkty[1] = 1. / sqrt(3.);

		wagi[0] = 1;
		wagi[1] = 1;
	}
	else if (liczba == 3)
	{
		punkty[0] = -sqrt(3. / 5.);
		punkty[1] = 0.;
		punkty[2] = sqrt(3. / 5.);

		wagi[0] = 5. / 9.;
		wagi[1] = 8. / 9.;
		wagi[2] = 5. / 9.;
	}
	else if (liczba == 4)
	{
		punkty[0] = -0.861136;
		punkty[1] = -0.339981;
		punkty[2] = 0.339981;
		punkty[3] = 0.861136;

		wagi[0] = (18. - sqrt(30.)) / 36.;
		wagi[1] = (18. + sqrt(30.)) / 36.;
		wagi[2] = (18. + sqrt(30.)) / 36.;
		wagi[3] = (18. - sqrt(30.)) / 36.;
	}
}

void zerowanie(uklad_rownan* uklad, dane* d,fun_ksztaltu* funkcje ,grid* siatka, globaldata* data)
{
	//zainicjowanie i wyzerowanie tablicy HG oraz C
	uklad->hg = new double* [siatka->nN];
	uklad->h = new double* [siatka->nN];
	uklad->C_agregowane = new double* [siatka->nN];
	uklad->HplusCdT = new double* [siatka->nN];

	for (int i = 0; i < siatka->nN; i++)
	{
		uklad->hg[i] = new double[siatka->nN + 1];
		uklad->h[i] = new double[siatka->nN + 1];
		uklad->C_agregowane[i] = new double[siatka->nN + 1];
		uklad->HplusCdT[i] = new double[siatka->nN + 1];

		for (int j = 0; j < siatka->nN; j++)
		{
			uklad->hg[i][j] = 0.;
			uklad->h[i][j] = 0.;
			uklad->C_agregowane[i][j] = 0.;
			uklad->HplusCdT[i][j] = 0.;
		}
	}


	d->P_agr = new double[siatka->nN];
	d->PplusCdT = new double[siatka->nN];
	d->t0 = new double[siatka->nN];
	for (int i = 0; i < siatka->nN; i++)
	{
		d->P_agr[i] = 0;
		d->PplusCdT[i] = 0;
		d->t0[i] = data->InitialTemp;
	}

	for (int i = 0; i < liczba * liczba; i++)
	{
		funkcje->tabela_dnPOksi[i] = new double[5];
		funkcje->tabela_dnPOeta[i] = new double[5];
		funkcje->fun_ksztaltu_dlaC[i] = new double[4];
	}
	for (int i = 0; i < 4 * liczba; i++) {
		funkcje->punkty_calkowania[i] = new double[2];
		funkcje->tablica_funkcja_ksztaltu[i] = new double[4];
	}
}

void liczenie_fun_ksztaltu(fun_ksztaltu* funkcje, double* punkty)
{
	for (int j = 0; j < liczba; j++)
	{
		for (int i = 0; i < liczba; i++)
		{
			funkcje->tabela_dnPOksi[j * liczba + i][0] = punkty[(j * liczba + i) / liczba];
			funkcje->tabela_dnPOksi[j * liczba + i][1] = n1odKsi(funkcje->tabela_dnPOksi[j * liczba + i][0]);
			funkcje->tabela_dnPOksi[j * liczba + i][2] = n2odKsi(funkcje->tabela_dnPOksi[j * liczba + i][0]);
			funkcje->tabela_dnPOksi[j * liczba + i][3] = n3odKsi(funkcje->tabela_dnPOksi[j * liczba + i][0]);
			funkcje->tabela_dnPOksi[j * liczba + i][4] = n4odKsi(funkcje->tabela_dnPOksi[j * liczba + i][0]);

			funkcje->tabela_dnPOeta[i * liczba + j][0] = punkty[j];
			funkcje->tabela_dnPOeta[i * liczba + j][1] = n1odEta(funkcje->tabela_dnPOeta[i * liczba + j][0]);
			funkcje->tabela_dnPOeta[i * liczba + j][2] = n2odEta(funkcje->tabela_dnPOeta[i * liczba + j][0]);
			funkcje->tabela_dnPOeta[i * liczba + j][3] = n3odEta(funkcje->tabela_dnPOeta[i * liczba + j][0]);
			funkcje->tabela_dnPOeta[i * liczba + j][4] = n4odEta(funkcje->tabela_dnPOeta[i * liczba + j][0]);


			//DLA MACIERZY C
			funkcje->fun_ksztaltu_dlaC[j * liczba + i][0] = n1(punkty[i], punkty[j]);
			funkcje->fun_ksztaltu_dlaC[j * liczba + i][1] = n2(punkty[i], punkty[j]);
			funkcje->fun_ksztaltu_dlaC[j * liczba + i][2] = n3(punkty[i], punkty[j]);
			funkcje->fun_ksztaltu_dlaC[j * liczba + i][3] = n4(punkty[i], punkty[j]);
		}
	}

	//DO HBC ORAZ P
	for (int i = 0; i < liczba * 4; i++)
	{
		funkcje->tablica_funkcja_ksztaltu[i][0] = n1(funkcje->punkty_calkowania[i][0], funkcje->punkty_calkowania[i][1]);
		funkcje->tablica_funkcja_ksztaltu[i][1] = n2(funkcje->punkty_calkowania[i][0], funkcje->punkty_calkowania[i][1]);
		funkcje->tablica_funkcja_ksztaltu[i][2] = n3(funkcje->punkty_calkowania[i][0], funkcje->punkty_calkowania[i][1]);
		funkcje->tablica_funkcja_ksztaltu[i][3] = n4(funkcje->punkty_calkowania[i][0], funkcje->punkty_calkowania[i][1]);
	}
}

void liczenie_warunkow_brzeg(fun_ksztaltu* funkcje, double* punkty)
{
	//LICZENIE WARUNKOW BRZEGOWYCH
	for (int i = 0; i < liczba; i++)
	{
		funkcje->punkty_calkowania[i][0] = punkty[i];
		funkcje->punkty_calkowania[i][1] = -1.;
	}

	for (int i = liczba; i < liczba * 2; i++)
	{
		funkcje->punkty_calkowania[i][0] = 1.;
		funkcje->punkty_calkowania[i][1] = punkty[i - liczba];
	}

	for (int i = liczba * 3 - 1, j = liczba * 2; i >= liczba * 2; i--, j++)
	{
		funkcje->punkty_calkowania[j][0] = punkty[i - liczba * 2];
		funkcje->punkty_calkowania[j][1] = 1.;
	}

	for (int i = liczba * 4 - 1, j = liczba * 3; i >= liczba * 3; i--, j++)
	{
		funkcje->punkty_calkowania[j][0] = -1.;
		funkcje->punkty_calkowania[j][1] = punkty[i - liczba * 3];
	}
}

void zerowanie_wewnatrz_petli(int ele, Elem* macierz)
{
	for (int i = 0; i < liczba * liczba; i++)
	{
		macierz[ele].tabela_dnPOdx[i] = new double[4];
		macierz[ele].tabela_dnPOdy[i] = new double[4];


		macierz[ele].h[i] = new double* [4];

		macierz[ele].jakobian[i] = new double[4];

		for (int j = 0; j < 4; j++)
		{

			macierz[ele].h[i][j] = new double[4];
		}
	}

	for (int k = 0; k < liczba * liczba; k++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int i = 0; i < 4; i++)
			{

				macierz[ele].h_koncowe[j][i] = 0;
				macierz[ele].hbc[j][i] = 0;
				macierz[ele].c[j][i] = 0;
			}
			macierz[ele].jakobian[k][j] = 0;
		}
	}
}

void liczenie_jakobianu(int ele, Elem* macierz, fun_ksztaltu* funkcje, grid* siatka)
{
	//LICZENIE JAKOBIANU
	for (int j = 0; j < liczba * liczba; j++)
	{
		for (int i = 0; i < 4; i++)
		{
			macierz[ele].jakobian[j][0] += funkcje->tabela_dnPOeta[j][i + 1] * siatka->ND[siatka->EL[ele].ID[i] - 1].y;
			macierz[ele].jakobian[j][1] += funkcje->tabela_dnPOksi[j][i + 1] * siatka->ND[siatka->EL[ele].ID[i] - 1].y;
			macierz[ele].jakobian[j][2] += funkcje->tabela_dnPOeta[j][i + 1] * siatka->ND[siatka->EL[ele].ID[i] - 1].x;
			macierz[ele].jakobian[j][3] += funkcje->tabela_dnPOksi[j][i + 1] * siatka->ND[siatka->EL[ele].ID[i] - 1].x;
		}
		macierz[ele].wyznacznik[j] = macierz[ele].jakobian[j][0] * macierz[ele].jakobian[j][3] - macierz[ele].jakobian[j][1] * macierz[ele].jakobian[j][2];
		//cout << macierz[ele].wyznacznik[j] << " ";
	}

	for (int j = 0; j < liczba * liczba; j++) //PRZEMNAZANIE JAKOBIANU PRZEZ 1/WYZNACZNIK
	{
		macierz[ele].jakobian[j][0] = macierz[ele].jakobian[j][0] * (1. / macierz[ele].wyznacznik[j]);
		macierz[ele].jakobian[j][1] = -macierz[ele].jakobian[j][1] * (1. / macierz[ele].wyznacznik[j]);
		macierz[ele].jakobian[j][2] = -macierz[ele].jakobian[j][2] * (1. / macierz[ele].wyznacznik[j]);
		macierz[ele].jakobian[j][3] = macierz[ele].jakobian[j][3] * (1. / macierz[ele].wyznacznik[j]);
	}
}

void get_macierz_h(int ele, Elem* macierz, fun_ksztaltu* funkcje, globaldata data, double* wagi)
{
	//MNOŻENIE MACIERZY

	for (int j = 0; j < liczba * liczba; j++)
	{
		for (int i = 0; i < 4; i++)
		{
			macierz[ele].tabela_dnPOdx[j][i] = macierz[ele].jakobian[j][0] * funkcje->tabela_dnPOksi[j][i + 1] + macierz[ele].jakobian[j][1] * funkcje->tabela_dnPOeta[j][i + 1];
			macierz[ele].tabela_dnPOdy[j][i] = macierz[ele].jakobian[j][2] * funkcje->tabela_dnPOksi[j][i + 1] + macierz[ele].jakobian[j][3] * funkcje->tabela_dnPOeta[j][i + 1];
		}
	}

	//UZYSKIWANIE MACIERZY H

	for (int k = 0; k < liczba * liczba; k++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int i = 0; i < 4; i++)
			{
				macierz[ele].h[k][j][i] = data.Conductivity * macierz[ele].wyznacznik[k] * (macierz[ele].tabela_dnPOdx[k][j] * macierz[ele].tabela_dnPOdx[k][i] + macierz[ele].tabela_dnPOdy[k][j] * macierz[ele].tabela_dnPOdy[k][i]);
			}
		}
	}

	//LICZENIE RAZEM Z WAGAMI
	for (int j = 0; j < 4; j++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int n = 0; n < liczba; n++)
			{
				for (int m = 0; m < liczba; m++)
				{
					macierz[ele].h_koncowe[j][i] += macierz[ele].h[(n * liczba) + m][j][i] * wagi[n] * wagi[m];
				}
			}
		}
	}
}

void oblicz_l(int ele, Elem* macierz,grid siatka)
{
	//LICZENIE L - DLUGOSCI BOKU
	macierz[ele].l[0] = sqrt(pow(siatka.ND[siatka.EL[ele].ID[1] - 1].x - siatka.ND[siatka.EL[ele].ID[0] - 1].x, 2) + pow(siatka.ND[siatka.EL[ele].ID[1] - 1].y - siatka.ND[siatka.EL[ele].ID[0] - 1].y, 2));
	macierz[ele].l[1] = sqrt(pow(siatka.ND[siatka.EL[ele].ID[2] - 1].x - siatka.ND[siatka.EL[ele].ID[1] - 1].x, 2) + pow(siatka.ND[siatka.EL[ele].ID[2] - 1].y - siatka.ND[siatka.EL[ele].ID[1] - 1].y, 2));
	macierz[ele].l[2] = sqrt(pow(siatka.ND[siatka.EL[ele].ID[3] - 1].x - siatka.ND[siatka.EL[ele].ID[2] - 1].x, 2) + pow(siatka.ND[siatka.EL[ele].ID[3] - 1].y - siatka.ND[siatka.EL[ele].ID[2] - 1].y, 2));
	macierz[ele].l[3] = sqrt(pow(siatka.ND[siatka.EL[ele].ID[0] - 1].x - siatka.ND[siatka.EL[ele].ID[3] - 1].x, 2) + pow(siatka.ND[siatka.EL[ele].ID[0] - 1].y - siatka.ND[siatka.EL[ele].ID[3] - 1].y, 2));
}

void oblicz_hbc_i_p(int ele, Elem* macierz, grid siatka, double* wagi, globaldata data, fun_ksztaltu funkcje)
{
	//TWORZENIE MACIERZY H DLA ELEMENTOW PRZED AGREGACJA
	for (int n = 0; n < 4; n++) //odpowiada za liczenie kazdej sciany z osobna
	{
		if (n != 3)
		{
			if (siatka.ND[siatka.EL[ele].ID[n] - 1].BC == 1 && siatka.ND[siatka.EL[ele].ID[n + 1] - 1].BC == 1)
			{
				for (int k = 0; k < liczba; k++) //odpowiada za liczenie kazdego wezla osobno
				{
					for (int j = 0; j < 4; j++) //j i i odpowiadaja za liczenie pojedynczej macierzy
					{
						for (int i = 0; i < 4; i++)
						{
							macierz[ele].hbc[j][i] += (macierz[ele].l[n] / 2) * wagi[k] * data.Alfa * funkcje.tablica_funkcja_ksztaltu[n * liczba + k][j] * funkcje.tablica_funkcja_ksztaltu[n * liczba + k][i];
						}
						macierz[ele].P[j] += data.Alfa * wagi[k] * funkcje.tablica_funkcja_ksztaltu[(n * liczba) + k][j] * data.Tot * macierz[ele].l[n] / 2;
					}
				}
			}
		}
		else
		{
			if (siatka.ND[siatka.EL[ele].ID[3] - 1].BC == 1 && siatka.ND[siatka.EL[ele].ID[0] - 1].BC == 1)
			{
				for (int k = 0; k < liczba; k++) //odpowiada za liczenie kazdego wezla osobno
				{
					for (int j = 0; j < 4; j++) //j i i odpowiadaja za liczenie pojedynczej macierzy
					{
						for (int i = 0; i < 4; i++)
						{
							macierz[ele].hbc[j][i] += (macierz[ele].l[n] / 2) * wagi[k] * data.Alfa * funkcje.tablica_funkcja_ksztaltu[n * liczba + k][j] * funkcje.tablica_funkcja_ksztaltu[n * liczba + k][i];
						}
						macierz[ele].P[j] += data.Alfa * wagi[k] * funkcje.tablica_funkcja_ksztaltu[(n * liczba) + k][j] * data.Tot * macierz[ele].l[n] / 2;
					}
				}
			}

		}
	}
}

void oblicz_c(int ele,Elem* macierz, globaldata data, fun_ksztaltu funkcje, double* wagi)
{
	//LICZENIE C
	for (int n = 0; n < liczba; n++)
	{
		for (int k = 0; k < liczba; k++)
		{
			for (int j = 0; j < 4; j++)
			{
				for (int i = 0; i < 4; i++)
				{
					macierz[ele].c[j][i] += data.Density * data.SpecificHeat * funkcje.fun_ksztaltu_dlaC[n * liczba + k][j] * funkcje.fun_ksztaltu_dlaC[n * liczba + k][i] * wagi[n] * wagi[k] * macierz[ele].wyznacznik[n * liczba + k];
				}
			}
		}
	}
}

void agregacja(int ele, Elem* macierz, uklad_rownan* uklad, grid siatka,dane* d)
{
	//DODAWANIE WARTOSCI DO DUZEJ MACIERZY HG - AGREGACJA
	for (int j = 0; j < 4; j++)
	{
		for (int i = 0; i < 4; i++)
		{
			uklad->hg[siatka.EL[ele].ID[j] - 1][siatka.EL[ele].ID[i] - 1] += macierz[ele].h_koncowe[j][i] + macierz[ele].hbc[j][i];
			uklad->C_agregowane[siatka.EL[ele].ID[j] - 1][siatka.EL[ele].ID[i] - 1] += macierz[ele].c[j][i];
		}
		d->P_agr[siatka.EL[ele].ID[j] - 1] += macierz[ele].P[j]; //AGREGACJA WEKTORA P
	}
}

void oblicz_zmod_h(grid siatka,uklad_rownan* uklad, globaldata data)
{
	//LICZENIE ZMODYFIKOWANEGO H
	for (int j = 0; j < siatka.nN; j++)
	{
		for (int i = 0; i < siatka.nN; i++)
		{
			uklad->HplusCdT[j][i] = uklad->hg[j][i] + uklad->C_agregowane[j][i] / data.SimulationStepTime;
		}
	}
}

void oblicz_zmod_p(grid siatka, uklad_rownan* uklad, globaldata data, dane* d)
{
	//LICZENIE ZMODYFIKOWANEGO P
	for (int j = 0; j < siatka.nN; j++)
	{
		for (int i = 0; i < siatka.nN; i++)
		{
			d->PplusCdT[j] += (uklad->C_agregowane[j][i] / data.SimulationStepTime) * d->t0[i];
		}
		d->PplusCdT[j] += d->P_agr[j];
	}
}

void oblicz_uklad_rownan(int liczba_rownan, uklad_rownan* uklad, dane* d, grid siatka, int step)
{
	//LICZENIE UKLADU ROWNAN
	double mnoznik, licznik;
	for (int i = 0; i < liczba_rownan - 1; i++)
	{
		for (int j = i + 1; j < liczba_rownan; j++)
		{
			mnoznik = uklad->HplusCdT[j][i] / uklad->HplusCdT[i][i];
			for (int k = i; k < liczba_rownan + 1; k++)
			{
				uklad->HplusCdT[j][k] = uklad->HplusCdT[j][k] - (uklad->HplusCdT[i][k] * mnoznik);
			}
		}
	}

	for (int i = liczba_rownan - 1; i >= 0; i--)
	{
		licznik = uklad->HplusCdT[i][liczba_rownan];
		for (int j = i + 1; j < liczba_rownan; j++)
		{
			licznik = licznik - (uklad->HplusCdT[i][j] * d->t0[j]);
		}
		d->t0[i] = licznik / uklad->HplusCdT[i][i];
	}

	cout << endl << "Temperatury dla t = " << step << "\t";
	cout << "Min: " << min(d->t0, siatka.nN) << "\t";
	cout << "Max: " << max(d->t0, siatka.nN) << endl;

	for (int i = 0; i < siatka.nN; i++)
	{
		d->PplusCdT[i] = 0;
	}
}

int main()
{
	string nazwa, x, y;
	string ele1, ele2, ele3, ele4;
	grid siatka;
	globaldata data;
	uklad_rownan uklad;
	dane d;
	fun_ksztaltu funkcje;

	cout.precision(6);

	wczytajDoPliku(&siatka, &data);

	Elem* macierz = new Elem[siatka.nE];

	double* punkty = new double[liczba];
	double* wagi = new double[liczba];

	przydzielWagiiPkt(punkty, wagi);
	
	zerowanie(&uklad, &d, &funkcje, &siatka, &data);

	liczenie_warunkow_brzeg(&funkcje, punkty);

	liczenie_fun_ksztaltu(&funkcje, punkty);

	for (int ele = 0; ele < siatka.nE; ele++)
	{
		zerowanie_wewnatrz_petli(ele, macierz);

		liczenie_jakobianu(ele, macierz, &funkcje, &siatka);

		get_macierz_h(ele, macierz, &funkcje, data, wagi);

		oblicz_l(ele, macierz, siatka);

		oblicz_hbc_i_p(ele, macierz, siatka, wagi, data, funkcje);

		oblicz_c(ele, macierz, data, funkcje, wagi);

		agregacja(ele, macierz, &uklad, siatka, &d);
	}

	for (int step = data.SimulationStepTime; step <= data.SimulationTime; step = step + data.SimulationStepTime)
	{
		oblicz_zmod_h(siatka, &uklad, data);

		oblicz_zmod_p(siatka, &uklad, data, &d);

		for (int i = 0; i < siatka.nN; i++)
		{
			uklad.HplusCdT[i][siatka.nN] = d.PplusCdT[i];
		}

		int liczba_rownan = siatka.nN;

		oblicz_uklad_rownan(liczba_rownan, &uklad, &d, siatka, step);
	}
	
}
