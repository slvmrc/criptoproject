#include"Bob2.h"

//bulinele 3, 6, 7, 8

typedef struct {
	int n, de;
} RSA_KEY;

void RSA_key_gen(int MIN, int MAX, RSA_KEY &p_key, RSA_KEY &s_key) {
	if (MAX > 40000)MAX = 40000;
	srand((int)time(NULL));
	long p, q, phi;
	p = da_prim(sqrt(MIN), sqrt(MAX), 0);
	q = da_prim(sqrt(MIN), sqrt(MAX), p);
	p_key.n = s_key.n = p * q;
	phi = (p - 1)*(q - 1);
	p_key.de = 3 + rand() % (phi - 3);
	s_key.de = invers(p_key.de, phi);
	while (s_key.de<0) {//pentru ca e sa fie inversabil
		p_key.de++;
		s_key.de = invers(p_key.de, phi);
	}
}

void RSA_cript_decript(const char* sursa, const char* destinatie, int lc, int ls, RSA_KEY &key) {
	int i, m;
	char *c = new char[lc >= ls ? lc : ls];
	ifstream in(sursa);
	ofstream out(destinatie);
	i = m = 0;
	while (in >> noskipws >> c[i]) {
		m = m * N + da_cod(c[i]);
		if (i == lc - 1) {
			m = a_la_b_mod_c(m, key.de, key.n);
			i = ls - 1;
			while (m>0) {
				c[i] = da_caracter(m%N);
				m = m / N;
				i--;
			}
			while (i >= 0)c[i--] = da_caracter(0);
			for (i = 0; i<ls; i++)out << c[i];
			i = 0;
		}
		else i++;
	}
	out.close();
	in.close();
	delete[] c;
}

void RSA_cript_cu_semnatura(const char* semnatura, const char* sursa, const char* destinatie, int lc, int ls, RSA_KEY
	p_key, RSA_KEY& s_p_key) {//semnatura de max 8 caractere
	const int DIM = 48 + 1;//dimensiunea blocului din final (cu semnatura)
	int i, m;
	char c[DIM];
	RSA_KEY s_s_key;
	RSA_key_gen(pow(N, lc) - 1, pow(N, ls), s_p_key, s_s_key);//generam cheile pentru semnatura
	ifstream in(sursa);
	ofstream temp("temp.txt");
	while (in.get(c, DIM)) {
		temp << c;
	}
	in.close();
	SHA1 sha1;
	unsigned long* rezultat = sha1.Valoare(sursa);
	sprintf(c, "%8s%08lx%08lx%08lx%08lx%08lx", semnatura, rezultat[0], rezultat[1], rezultat[2],
		rezultat[3], rezultat[4]);
	for (int i = 0; i + lc < DIM; i += lc) {
		int m = 0, j;
		char bloc[10];
		for (j = 0; j<lc; j++)m = m * N + da_cod(c[i + j]);
		m = a_la_b_mod_c(m, s_s_key.de, s_s_key.n);
		j = ls - 1;
		while (m>0) {
			bloc[j] = da_caracter(m%N);
			m = m / N;
			j--;
		}
		while (j >= 0)bloc[j--] = da_caracter(0);
		for (j = 0; j<ls; j++)temp << bloc[j];
		j = 0;
	}
	in.close();
	temp.close();
	RSA_cript_decript("temp.txt", "destinatie.txt", lc, ls, p_key);
	remove("temp.txt");//stergem fisierul temporar
}

void RSA_verifica_semnatura(const char* sursa, int DIM_SEMNATURA, int lc, int ls, RSA_KEY s_p_key) {
	char* c = new char[DIM_SEMNATURA + 1], *semnatura = new char[DIM_SEMNATURA + 1];
	int index = 0;
	FILE *in = fopen(sursa, "r");
	fseek(in, 0, SEEK_END);
	long lungime = ftell(in);
	fseek(in, (lungime - DIM_SEMNATURA), SEEK_SET);
	for (int i = 0; i<DIM_SEMNATURA; i++)c[i] = fgetc(in);
	fclose(in);
	for (int i = 0; i + ls < DIM_SEMNATURA; i += ls) {
		int m = 0, j;
		char bloc[10];
		for (j = 0; j<ls; j++)m = m * N + da_cod(c[i + j]);
		m = a_la_b_mod_c(m, s_p_key.de, s_p_key.n);
		j = lc - 1;
		while (m>0) {
			bloc[j] = da_caracter(m%N);
			m = m / N;
			j--;
		}
		while (j >= 0)bloc[j--] = da_caracter(0);
		for (j = 0; j<lc; j++) semnatura[index++] = bloc[j];
		j = 0;
	}
	semnatura[index] = '\0';
	cout << semnatura << endl;
	delete[] c;
	delete[] semnatura;
}

int main() {

		int n = 6731, e = 4393;
		int p, q;
		factorizare(n, p, q);
		cout << p << " " << q << endl;
		int phi_n = (p - 1)*(q - 1);
		int d = invers(e, phi_n); //cheia secreta pt RSA
		cout << d;
		
		citeste_alfabet();
		RSA_KEY p_key, s_key;
		p_key.n = s_key.n = 31313; //obt in urma apelului de fct rsa_key_gen
		p_key.de = 19777;
		s_key.de = 10033;
		int lc = 2, ls = 3;

		//Am pus urmatoarele linii de cod in comentariu pentru a nu mai genera noi chei, dar prima data am rulat programul cu ele ! 
		//RSA_key_gen(pow(N, lc) - 1, pow(N, ls), p_key, s_key); 
		//out << "\nCheia publica--> Ke=(" << p_key.n << "," << p_key.de << ")";
		//out << "\nCheia secreta--> Kd=(" << s_key.n << "," << s_key.de << ")" << endl;
		//ofstream out("chei.txt");

		RSA_cript_decript("sursa.txt", "cript.txt", lc , ls, p_key);
		RSA_cript_cu_semnatura("sevilla", "sursa.txt", "destinatie.txt", lc, ls, p_key, s_key);
}
