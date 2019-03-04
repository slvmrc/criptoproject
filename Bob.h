// CriptTools.h
#ifndef CRIPT_TOOLS_H
#define CRIPT_TOOLS_H
#pragma warning(disable: 4996)
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <iomanip>
#include <bitset>
#include <thread>
using namespace std;

#define MyInt int 

struct  Bloc {
	char c;
	int nr;
};
/*variabile globale*/
char caracter[100] = { 0 };
int N = 0;
//====================
int cmmdc(MyInt x, MyInt y) {
	MyInt rest;
	do {
		rest = x % y;
		x = y;
		y = rest;
	} while (rest != 0);
	return x;
}

/*extindem operatorul modulo (%) si pentru numere negative*/
int modulo(MyInt k, MyInt n) {
	if (k<0)k = n - (-k) % n;
	if (k >= n) return k % n;
	return k;
}
/*va returna -1 daca numarul nu este inversabil*/
int invers(MyInt a, MyInt n) {
	MyInt q, r, x0 = 1, x1 = 0, copy_n = n;
	a = modulo(a, n);
	while (n != 0)
	{
		r = n;
		q = a / n;
		n = a % n;
		a = r;

		r = x1;
		x1 = x0 - q * x1;
		x0 = r;
	}
	if (a == 1)
		return modulo(x0, copy_n);
	return -1;
}

/*aloca memorie matrice*/
void creare(MyInt**& matrice, int lin, int col) {
	matrice = new MyInt*[lin];
	for (int i = 0; i < lin; i++) {
		matrice[i] = new MyInt[col];
	}
}
/*aloca memorie vector*/
void creare(MyInt*& vector, int dim) {
	vector = new MyInt[dim];
}
/*afiseaza o matrice*/
void scrie(MyInt**& matrice, int lin, int col) {
	for (int i = 0; i < lin; i++) {
		for (int j = 0; j < col; j++)
			cout << matrice[i][j] << " ";
		cout << endl;
	}
}
/*afiseaza un vector*/
void scrie(MyInt*& vector, int dim) {
	for (int i = 0; i < dim; i++) {
		cout << vector[i] << " ";
	}
}
/*citeste o matrice*/
void citeste(MyInt**& matrice, int lin, int col) {
	cout << "Dati elementele matricei (" << lin << "x" << col << ")" << endl;
	for (int i = 0; i < lin; i++) {
		for (int j = 0; j < col; j++)
			cin >> matrice[i][j];
	}
}
/*citeste un vector*/
void citeste(MyInt*& vector, int dim) {
	cout << "Dati elementele vectorului (" << dim << "elemente)" << endl;
	for (int i = 0; i < dim; i++) {
		cin >> vector[i];
	}
}
/*eliberare memorie matrice*/
void stergere(MyInt**& matrice, int lin, int col) {
	for (int i = 0; i < lin; i++) {
		delete[] matrice[i];
	}
	delete[] matrice;
}
/*eliberare memorie vector*/
void stergere(MyInt*& vector, int dim) {
	delete[] vector;
}
/*modulo pentru matrice*/
void modulo(int**& matrice, int lin, int col, int n) {
	for (int i = 0; i<lin; i++)
		for (int j = 0; j < col; j++)
			matrice[i][j] = modulo(matrice[i][j], n);
}

/*calculeaza suma a doi vectori*/
void suma(int *vector1, int *vector2, int dim, int *rezultat) {
	for (int i = 0; i<dim; ++i)
		rezultat[i] = vector1[i] + vector2[i];
}

/*calculeaza suma a doua matrice*/
void suma(int **matrice1, int **matrice2, int lin, int col, int **rezultat) {
	for (int i = 0; i<lin; ++i)
		for (int j = 0; j<col; ++j)
			rezultat[i][j] = matrice1[i][j] + matrice2[i][j];
}

/*calculeaza produsul a doua matrice dim1xdim2  dim2xdim3*/
void produs(int **matrice1, int **matrice2, int dim1, int dim2, int dim3, int **rezultat) {
	int** temp;
	creare(temp, dim1, dim3);
	for (int i = 0; i < dim1; i++)
		for (int j = 0; j < dim3; j++) {
			temp[i][j] = 0;
			for (int k = 0; k < dim2; k++)
				temp[i][j] += matrice1[i][k] * matrice2[k][j];
		}
	for (int i = 0; i < dim1; i++)
		for (int j = 0; j < dim3; j++)
			rezultat[i][j] = temp[i][j];
	stergere(temp, dim1, dim3);
}

/*calculeaza produsul unei matrice linxcol cu un vector*/
void produs(int **matrice, int *vector, int lin, int col, int *rezultat) {
	int* temp;
	creare(temp, col);
	for (int i = 0; i < lin; i++) {
		temp[i] = 0;
		for (int j = 0; j < col; j++)
			temp[i] += matrice[i][j] * vector[j];
	}
	for (int i = 0; i < col; i++)
		rezultat[i] = temp[i];
	stergere(temp, col);
}

/*calculeaza minorul corespunzator pentru lin, col din matricea matrice*/
void calcul_minor(int **matrice, int lin, int col, int n, int **rezultat) {
	for (int i = 0; i < n - 1; i++) {
		for (int j = 0; j < n - 1; j++) {
			if (i < lin) {
				if (j < col)
					rezultat[i][j] = matrice[i][j];
				else
					rezultat[i][j] = matrice[i][j + 1];
			}
			else {
				if (j < col)
					rezultat[i][j] = matrice[i + 1][j];
				else
					rezultat[i][j] = matrice[i + 1][j + 1];
			}
		}
	}
}

/*calculeaza valoarea determinantului matricei matrice*/
int calcul_det(int** matrice, int dim) {
	if (dim <= 1) return matrice[0][0];
	int **a = new int *[dim - 1];
	for (int i = 0; i< dim - 1; i++) {
		a[i] = new int[dim - 1];
	}
	int S = 0;
	for (int i = 0; i < dim; i++) {
		calcul_minor(matrice, 0, i, dim, a);
		S += matrice[0][i] * (i % 2 ? -1 : 1)*calcul_det(a, dim - 1);
	}
	for (int i = 0; i < dim - 1; i++)
		delete[] a[i];
	return S;
}

/*calculeaza inversa unei matrice in Zn*/
void invers(int **matrice, int dim, int n, int **rezultat) {
	int d = invers(calcul_det(matrice, dim), n);
	if (d < 0 || dim == 1) {
		rezultat[0][0] = d;
		return;
	}
	int **a, **temp;
	creare(a, dim - 1, dim - 1);
	creare(temp, dim, dim);
	for (int i = 0; i<dim; i++)
		for (int j = 0; j < dim; j++) {
			calcul_minor(matrice, j, i, dim, a);
			temp[i][j] = modulo(d*((i + j) % 2 ? -1 : 1)*calcul_det(a, dim - 1), n);
		}
	for (int i = 0; i<dim; i++)
		for (int j = 0; j < dim; j++) {
			rezultat[i][j] = temp[i][j];
		}
	stergere(a, dim - 1, dim - 1);
	stergere(temp, dim, dim);

}

int a_la_b_mod_c(int a, int b, int c) {
	int p = 1;
	a %= c;
	while (b>0)
	{
		if (b % 2)
			p = (p*a) % c;
		a = (a*a) % c;
		b /= 2;
	}
	return p;
}

int prim(int x) {
	int i;
	if ((x == 1) || (x == 2)) return 1;
	if (x % 2 == 0) return 0;
	for (i = 3; i <= sqrt((double)x); i += 2)
		if (x%i == 0) return 0;
	return 1;
}

//da un numar prim diferit de nr. situat intre min si max.
//Generam numere astfel incat produsul a doua asemenea numere sa nu depaseasca domeniul pentru variabile te tip int
int da_prim(int min, int max, int nr) {
	int k, nrIncercari = 10;
	if (min<0 || max<0 || min>max)return-1;
	if (min == max)
		if (prim(min))
			return min;
		else
			return -1;
	for (int i = 0; i<nrIncercari; i++) {
		k = (int)sqrt((double)rand());
		k = min + k % (max - min);
		while (!prim(k) || k == nr)k++;
		if (k <= max)return k;
	}
	k = min;
	while (!prim(k) || k == nr)k++;
	if (k <= max)return k;
	return -1;
}

int val_pol(int *coef, int grad, int x, int p) {
	int S = 0;
	for (int i = grad; i >= 0; i--)
		S = modulo((coef[i] + S * x), p);
	return S;
}

long este_patrat_perfect(long x) {
	long temp = (long)sqrt((double)x);
	if (temp*temp == x)return temp;
	return 0;
}

void factorizare(int n, int &p, int &q) {
	long s_patrat, t;
	t = (long)(sqrt((double)n) + 1);
	s_patrat = t * t - n;
	while (!este_patrat_perfect(s_patrat) && (t <= n)) {
		t++;
		s_patrat = t * t - n;
	}
	if (este_patrat_perfect(s_patrat) >= 0) {
		s_patrat = (long)sqrt((double)s_patrat);
		p = t + s_patrat;
		q = t - s_patrat;
	}
	else {
		p = -1; q = -1;
	}
}

//pentru logaritmul discret
int log_d(int g, int y, int p) {
	int i, j, *v, a, b;
	a = (int)sqrt((double)p);
	v = new int[a + 1];
	v[0] = 1;
	v[1] = a_la_b_mod_c(g, a, p);

	for (i = 2; i <= a; i++)
		v[i] = (v[i - 1] * v[1]) % p;
	b = modulo(y, p);
	for (i = 0; i<a; i++) {
		for (j = 0; j <= a; j++)
			if (b == v[j])return modulo(j*a - i, p - 1);
		b = (b*g) % p;
	}
	return 0;
}
/*citeste caracterele alfabetului precum si numarul de caractere din alfabet*/
void citeste_alfabet() {
	N = 0;
	ifstream in("alfabet.txt");
	if (!in.good())
		perror("Fisier inexistent");
	char c;
	while (in >> noskipws >> c) {
		caracter[N] = c;
		N++;
	}
	if (N == 0)
		cout << "Alfabetul dat are 0 caractere" << endl;
}

/*da codul(indexul) caracterului c in tabloul caracter */
int da_cod(char c) {
	for (int i = 0; i<N; i++)
		if (caracter[i] == c)return i;
	return -1;
}

/*da caracterul corespunzator codului dat*/
char da_caracter(int cod) {
	return caracter[modulo(cod, N)];
}

/*transforma in baza 10 un text scris cu alfabetul dat*/
int in_baza_10(char* text) {
	int i = 0, rez = 0;
	while (text[i] != '\0')
		rez = rez * N + da_cod(text[i++]);
	return rez;
}

/*transforma un numar scris in baza 10 intr-un text scris cu alfabetul dat*/
void din_baza_10(int nr, char text[]) {
	int i = 0;
	while (nr>0) {
		text[i++] = da_caracter(nr%N);
		nr /= N;
	}
	text[i] = '\0';
	_strrev(text);
}


/*analiza frecventelor pentru blocuri de lungime 1*/
void frecvente(ifstream &in, Bloc bloc[]) {
	for (int i = 0; i<N; i++) {
		bloc[i].c = da_caracter(i);
		bloc[i].nr = 0;
	}
	char c;
	while (in >> noskipws >> c) {
		int cod = da_cod(c);
		if (cod >= 0)
			bloc[cod].nr++;
	}
	int ordonat = 0;
	while (!ordonat) {
		ordonat = 1;
		for (int i = 0; i<N - 1; i++) {
			if (bloc[i].nr<bloc[i + 1].nr) {
				ordonat = 0;
				Bloc temp = bloc[i];
				bloc[i] = bloc[i + 1];
				bloc[i + 1] = temp;
			}
		}
	}
}

/*pentru lucrul cu numere mari*/
template<unsigned int N>
class BigInt {
private:
	bitset<N> bits;

public:
	BigInt() : bits() { }
	BigInt(const BigInt& x) : bits(x.bits) { }
	BigInt(unsigned long x) {
		int n = 0;
		while (x) {
			bits[n++] = x & (unsigned long)1;
			x >>= 1;
		}
	}
	BigInt(int x) :BigInt((unsigned long)x) {}
	explicit BigInt(const bool b[], unsigned int dim) {
		bits.reset();
		if (dim > N)dim = N;
		for (int i = 0; i < N && i < dim; i++)bits[dim - 1 - i] = b[i];
	}
	explicit BigInt(const bitset<N>& x) : bits(x) { }
	BigInt(const char* x)
	{
		int indexBit = 0;
		char p[1 + N / 3];
		strncpy(p, x, N / 3);//luam sa incapa
		p[N / 3] = '\0';
		while (strlen(p)) {
			bits[indexBit++] = p[strlen(p) - 1] % 2;
			int temp = 0, aux = 0;
			for (int i = 0; i < strlen(p); i++) {
				aux = p[i] - '0' + temp * 10;
				p[i] = '0' + aux / 2;
				temp = aux % 2;
			}
			if (p[0] == '0')
				strcpy(p, p + 1);
		}
	}

	static void divide(BigInt x, BigInt y, BigInt& q, BigInt& r) {
		if (y == 0) throw std::domain_error("division by zero undefined");
		q = r = 0;
		if (x == 0) return;
		if (x == y) {
			q.bits[0] = 1;
			return;
		}
		r = x;
		if (x < y) return;

		// determinam bitul cel mai semnificativ in x si y
		unsigned int sig_x;
		for (int i = N - 1; i >= 0; i--) {
			sig_x = i;
			if (x[i]) break;
		}
		unsigned int sig_y;
		for (int i = N - 1; i >= 0; i--) {
			sig_y = i;
			if (y[i]) break;
		}
		// aliniem x si y
		unsigned int n = (sig_x - sig_y);
		y <<= n;
		n += 1;
		// algoritmul: deplasare si scadere 
		while (n--) {
			if (y <= r) {
				q.bits[n] = true;
				r -= y;
			}
			y >>= 1;
		}
	}
	static BigInt sqrt(BigInt a) {
		BigInt x0 = a, x1 = (a + 1) / 2;
		while (x1 < x0) {
			x0 = x1;
			x1 = (x1 + a / x1) / 2;
		}
		return x0;
	}
	static bool isSqrt(BigInt a) {
		BigInt x0 = a, x1 = (a + 1) / 2;
		while (x1 < x0) {
			x0 = x1;
			x1 = (x1 + a / x1) / 2;
		}
		return x0 * x0 == a;
	}
	static BigInt pow(BigInt a, BigInt b) {
		if (b == 0) return 1;
		BigInt temp = pow(a, b >> 1);
		if (b[0] == 0) return temp * temp;
		return temp * temp * a;
	}
	static unsigned int log(BigInt a) { // log_2(a)
		for (unsigned int i = N - 1; i >= 0; i--)
			if (a.bits[i] != 0)return i;
		throw std::domain_error("log invalid argument");
	}
	void factorizare(BigInt &p, BigInt &q) {
		BigInt n = *this;
		BigInt s_patrat, t;
		t = sqrt(n) + 1;
		s_patrat = t * t - n;
		while (!isSqrt(s_patrat) && (t <= n)) {
			t++;
			s_patrat = t * t - n;
		}
		if (isSqrt(s_patrat) >= 0) {
			s_patrat = sqrt(s_patrat);
			p = t + s_patrat;
			q = t - s_patrat;
		}
		else {
			p = 0; q = 0;
		}
	}
	void scrie(ostream& out) {
		unsigned char c = 0;
		unsigned long i;
		for (i = 0; i < N; i++) {
			c <<= 1;
			if (bits[N - i - 1])c++;
			if (i % 8 == 7) {
				out << c;
				c = 0;
			}
		}
		if (i % 8)out << c;
	}
	bool citeste(istream& in) {
		unsigned char c = 0;
		unsigned long i = 0;
		bits.reset();
		while (i<N && in >> noskipws >> c) {
			for (int j = 7; i<N && j >= 0; j--) {
				bits[N - 1 - i++] = c & (1 << j);
			}
		}
		return !in.eof();
	}
	static bool este_supercrescator(BigInt v[], const int dim) {
		if (dim <= 0)return true;
		BigInt S = v[0];
		for (int i = 1; i < dim; i++) {
			if (S >= v[i])return false;
			S += v[i];
		}
		return true;
	}
	static bool rezolva_rucsac(BigInt V, BigInt v[], bool eps[], int dim) {
		//punem rezultatul in eps
		for (int i = 0; i < dim; i++)
			eps[i] = 0;
		for (int i = dim - 1; i >= 0; i--) {
			if (V >= v[i]) {
				eps[i] = 1;
				V -= v[i];
				if (V == 0)return true;
			}
		}
		return V == 0;
	}
	BigInt invers(BigInt n) {
		BigInt q, r, x0 = 1, x1 = 0, copy_n = n, a = *this % n;
		while (n != 0)
		{
			r = n;
			q = a / n;
			n = a % n;
			a = r;

			r = x1;
			x1 = q * x1;
			while (x0 < x1)x0 += copy_n;

			x1 = x0 - x1;
			x0 = r;
		}
		if (a == BigInt(1))
			return x0 % copy_n;
		return 0;
	}
	bitset<N> Bits() {
		return bits;
	}
	bool operator[](int n) const {
		return bits[n];
	}
	operator string() const {
		int nrBlocuri = 1 + N / 3;
		bitset<4> bloc[1 + N / 3];

		for (int i = N - 1; i >= 0; i--) {
			for (int j = 0; j < nrBlocuri; j++) {
				if (bloc[j].to_ulong() >= 5)bloc[j] = bloc[j].to_ulong() + 3;
			}
			for (int j = 0; j < nrBlocuri - 1; j++) {
				bloc[j] <<= 1;
				bloc[j][0] = bloc[j + 1][3];
			}
			bloc[nrBlocuri - 1] <<= 1;
			bloc[nrBlocuri - 1][0] = bits[i];
		}
		string rez = "";
		bool inceput = false;
		for (int i = 0; i < nrBlocuri; i++) {
			if (bloc[i].to_ulong() > 0)inceput = true;
			if (inceput)rez.append(to_string(bloc[i].to_ulong()));
		}
		return inceput ? rez : "0";
	}
	BigInt& operator<<=(unsigned int n) {
		bits <<= n;
		return *this;
	}
	BigInt& operator>>=(unsigned int n) {
		bits >>= n;
		return *this;
	}
	BigInt operator++(int) {
		BigInt temp = *this;
		operator++();
		return temp;
	}
	BigInt operator--(int) {
		BigInt temp = *this;
		operator--();
		return temp;
	}
	BigInt& operator++() {
		bool temp = false;
		bool suma = ~bits[0] ^ temp;
		temp = bits[0] || temp || (bits[0] && temp);
		bits[0] = suma;
		for (int i = 1; i < N; i++) {
			suma = bits[i] ^ temp;
			temp = bits[i] && temp;
			bits[i] = suma;
		}
		return *this;
	}
	BigInt& operator--() {
		bool temp = !bits[0];
		bits[0] = ~bits[0];
		for (int i = 1; i < N; i++) {
			if (temp) {
				temp = !bits[i];
				bits[i] = !bits[i];
			}
			else {
				temp = 0;
			}
		}
		return *this;
	}
	BigInt& operator+=(const BigInt& x) {
		bool temp = false;
		for (int i = 0; i < N; i++) {
			bool suma = (bits[i] ^ x.bits[i]) ^ temp;
			temp = (bits[i] && x.bits[i]) || (bits[i] && temp) || (x.bits[i] && temp);
			bits[i] = suma;
		}
		return *this;
	}
	BigInt& operator-=(const BigInt& x) {
		bool borrow = false;
		for (int i = 0; i < N; i++) {
			if (borrow) {
				if (bits[i]) {
					bits[i] = x.bits[i];
					borrow = x.bits[i];
				}
				else {
					bits[i] = !x.bits[i];
					borrow = true;
				}
			}
			else {
				if (bits[i]) {
					bits[i] = !x.bits[i];
					borrow = false;
				}
				else {
					bits[i] = x.bits[i];
					borrow = x.bits[i];
				}
			}
		}
		return *this;
	}
	BigInt& operator*=(const BigInt& x) {
		BigInt temp = *this;
		*this = 0;
		if (temp.bits.count() < x.bits.count()) {
			for (unsigned int i = 0; i < N; i++)
				if (temp[i]) *this += x << i;
		}
		else {
			for (unsigned int i = 0; i < N; i++)
				if (x[i]) *this += (temp << i);
		}
		return *this;
	}
	BigInt& operator/=(const BigInt& x) {
		BigInt temp;
		divide(*this, x, *this, temp);
		return *this;
	}
	BigInt& operator%=(const BigInt& x) {
		BigInt temp;
		divide(*this, x, temp, *this);
		return *this;
	}
	BigInt operator~() const {
		return ~bits;
	}
	BigInt& operator&=(BigInt x) {
		bits &= x.bits;
		return *this;
	}
	BigInt& operator|=(BigInt x) {
		bits |= x.bits;
		return *this;
	}
	BigInt& operator^=(BigInt x) {
		bits ^= x.bits;
		return *this;
	}
	friend BigInt operator<<(BigInt x, unsigned int n) {
		return x <<= n;
	}
	friend BigInt operator >> (BigInt x, unsigned int n) {
		return x >>= n;
	}
	friend BigInt operator+(BigInt x, const BigInt& y) {
		return x += y;
	}
	friend BigInt operator-(BigInt x, const BigInt& y) {
		return x -= y;
	}
	friend BigInt operator*(BigInt x, const BigInt& y) {
		return x *= y;
	}
	friend BigInt operator/(BigInt x, const BigInt& y) {
		return x /= y;
	}
	friend BigInt operator%(BigInt x, const BigInt& y) {
		return x %= y;
	}
	friend BigInt operator^(BigInt x, const BigInt& y) {
		return x ^= y;
	}
	friend BigInt operator&(BigInt x, const BigInt& y) {
		return x &= y;
	}
	friend BigInt operator|(BigInt x, const BigInt& y) {
		return x |= y;
	}
	friend bool operator==(const BigInt& x, const BigInt& y) {
		return x.bits == y.bits;
	}
	friend bool operator!=(const BigInt& x, const BigInt& y) {
		return x.bits != y.bits;
	}
	friend bool operator>(const BigInt& x, const BigInt& y) {
		for (int i = N - 1; i >= 0; i--) {
			if (x[i] && !y[i]) return true;
			if (!x[i] && y[i]) return false;
		}
		return false;
	}
	friend bool operator<(const BigInt& x, const BigInt& y) {
		for (int i = N - 1; i >= 0; i--) {
			if (x[i] && !y[i]) return false;
			if (!x[i] && y[i]) return true;
		}
		return false;
	}
	friend bool operator>=(const BigInt& x, const BigInt& y) {
		for (int i = N - 1; i >= 0; i--) {
			if (x[i] && !y[i]) return true;
			if (!x[i] && y[i]) return false;
		}
		return true;
	}
	friend bool operator<=(const BigInt& x, const BigInt& y) {
		for (int i = N - 1; i >= 0; i--) {
			if (x[i] && !y[i]) return false;
			if (!x[i] && y[i]) return true;
		}
		return true;
	}
	friend istream& operator >> (istream& in, BigInt& n) {
		in >> n.bits;
		return in;
	}
	friend ostream& operator << (ostream& out, const BigInt& n) {
		out << (string)n;
		return out;
	}
	friend class CitesteBinar;
};

class CitesteBinar {
	ifstream in;
	long dimMax, dimCitit;
	unsigned char c;
	int index;
public:
	CitesteBinar(char* sursa) {
		in = ifstream(sursa, std::ios::binary);
		c = 0;
		dimMax = dimCitit = 0;
		//citim lungimea fisierului
		if (in) {
			in.seekg(0, in.end);
			dimMax = in.tellg();
			in.seekg(0, in.beg);
		}
		index = -1;
	}
	int citeste(bool bts[], int max) {
		int i = 0;
		while (i<max && index >= 0) {
			bts[i++] = c & (1 << index--);
		}
		while (i<max && in >> noskipws >> c) {
			dimCitit++;
			for (index = 7; index >= 0 && i<max; index--) {
				bts[i++] = c & (1 << index);
			}
		}
		return i;
	}
	template<unsigned int N>
	int citeste(BigInt<N> &V) {
		V.bits.reset();
		int i = 0;
		while (i<N && index >= 0) {
			V.bits[N - 1 - i++] = c & (1 << index--);
		}
		while (i<N && in >> noskipws >> c) {
			dimCitit++;
			for (index = 7; index >= 0 && i<N; index--) {
				V.bits[N - 1 - i++] = c & (1 << index);
			}
		}
		return i;
	}
	void scrieProcent() {
		int procent = (dimMax == 0) ? 100 : 100 * dimCitit / dimMax;
		printf("\r%3d%% [%.*s%*s]", procent, procent / 2, "||||||||||||||||||||||||||||||||||||||||||||||||||", 50 - procent / 2, "");
		fflush(stdout);
	}
	void close() {
		if (in)in.close();
	}
};
//criptare Merkle-Hellman
template<int N>
void criptareMH(char* sursa, char* destinatie, BigInt<N> w[], int nr_elem) {
	CitesteBinar cb(sursa);
	ofstream out(destinatie, std::ios::binary);
	int k;
	bool* bts = new bool[nr_elem];
	while ((k = cb.citeste(bts, nr_elem)) > 0) {
		BigInt<N> S = 0;
		for (int i = 0; i < k; i++)
			if (bts[i])S += w[i];
		S.scrie(out);
		cb.scrieProcent();
	}
	delete[] bts;
	cb.close();
	out.close();
}


class SHA1 {
	unsigned long l1, l2;//pentru lungimea textului text (numarul de biti). Nr. total de biti va fi l1*2^32+l2<2^64
	unsigned long *h;
public:
	SHA1() {
		l1 = l2 = 0;
		h = new unsigned long[5];
		h[0] = 0x67452301;//valorile initiale
		h[1] = 0xEFCDAB89;
		h[2] = 0x98BADCFE;
		h[3] = 0x10325476;
		h[4] = 0xC3D2E1F0;
	}
	void aduna_lungime(unsigned long valoare) {//se aduna valoare la valoarea pastrata in l1 si l2
		unsigned long x = 0xFFFFFFFF - l2;
		if (x >= valoare) {//daca valoare incape in l2
			l2 += valoare;
		}
		else {
			l1++;
			l2 = valoare - x - 1;
		}
	}
	unsigned long leftrotate(unsigned long x, int n) {
		return (x << n) | (x >> (32 - n));
	}
	/*	unsigned long xor(unsigned long x, unsigned long y) {//se poate folosi in locul operatorului ^
	return (x | y)&(~x | ~y);
	}*/
	void transforma_h(unsigned long* w) {//w are 80 elemente
		unsigned long a, b, c, d, e, f, k, temp;
		a = h[0];
		b = h[1];
		c = h[2];
		d = h[3];
		e = h[4];
		for (int i = 16; i<80; i++)
			w[i] = leftrotate(w[i - 3] ^ w[i - 8] ^ w[i - 14] ^ w[i - 16], 1);
		for (int i = 0; i<80; i++) {
			if (0 <= i && i<20) {
				f = (b&c) | ((~b)&d);
				k = 0x5A827999;
			}
			else if (20 <= i && i<40) {
				f = b ^ c^d;
				k = 0x6ED9EBA1;
			}
			else if (40 <= i && i<60) {
				f = (b&c) | (b&d) | (c&d);
				k = 0x8F1BBCDC;
			}
			else if (60 <= i && i<80) {
				f = b ^ c^d;
				k = 0xCA62C1D6;
			}
			temp = leftrotate(a, 5) + f + e + k + w[i];
			e = d;
			d = c;
			c = leftrotate(b, 30);
			b = a;
			a = temp;
		}
		h[0] = h[0] + a;
		h[1] = h[1] + b;
		h[2] = h[2] + c;
		h[3] = h[3] + d;
		h[4] = h[4] + e;
	}

	void unsigned_long_to_bin(unsigned long l, char* rez) {
		for (int i = 7; i >= 0; i--) { //doar un for de la 31 la 0 //cont sf de sir, da lg diferita
			for (int j = 3; j >= 0; j--)
				rez[3 - j + i * 4] = ((l % 16) & (1 << j)) ? '1' : '0';
			l /= 16;
		}
	}
	unsigned long* Valoare(const char* sursa) { //apelata in while
		ifstream in;
		in.open(sursa, std::ios_base::binary);
		in.seekg(0, ios::end);
		int n = in.tellg();
		in.seekg(0, ios::beg);
		char* res = new char[n];
		in.read(res, n);
		in.close();
		return Valoare(res, n);
		//delete res;
	}
	char* ValoareBin(char* text, int len) {
		char* rez = new char[161]; //160 caract {0,1}; //atentie, pb if sirul comparat n-are sf de sir
		unsigned long* hash = Valoare(text, len);
		for (int i = 0; i < 5; i++)
			unsigned_long_to_bin(hash[i], rez + i * 32); //scrie din 32 in 32 de biti
		rez[160] = '\0';
		return rez;
	}
	unsigned long* Valoare(char* text, int len) { //cripteaza exact textul-sa testam validitatea cheilor
		unsigned long* w = new unsigned long[80]; //o apelam in for
		unsigned char* c = new unsigned char[4];//din 4 char facem un double
		int i = 0, index = 0;
		while (index<len) {
			c[i % 4] = text[index++];
			i++;
			aduna_lungime(8);//am citit 8 biti
			if (i % 4 == 0) {
				w[i / 4 - 1] = (c[0] << 24) | (c[1] << 16) | (c[2] << 8) | c[3];
				if (i % 64 == 0) {
					transforma_h(w);
					i = 0;
				}
			}
		}
		//completam corespunzator
		c[i % 4] = 1 << 7;
		i++;
		while ((i % 4)>0) {
			c[i % 4] = 0;
			i++;
		}
		w[i / 4 - 1] = (c[0] << 24) | (c[1] << 16) | (c[2] << 8) | c[3];
		if (i <= 56) {//incap in acest bloc cei 64 biti in care trecem lungimea blocului
			while (i<56) {
				i += 4;
				w[i / 4 - 1] = 0;
			}
			w[14] = l1;
			w[15] = l2;
			transforma_h(w);
		}
		else {//nu incap cei 64 biti, vor fi trecuti in blocul urmator
			while (i<64) {
				i += 4;
				w[i / 4 - 1] = 0;
			}
			transforma_h(w);
			for (i = 4; i<56; i += 4)
				w[i / 4 - 1] = 0;
			w[14] = l1;
			w[15] = l2;
			transforma_h(w);
		}
		return h;
	}
};
#endif
