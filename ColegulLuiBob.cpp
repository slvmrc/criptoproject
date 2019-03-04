#include"ColegulLuiBob.h"

//bulinele 4 si 5

char string_din_nr(int nr) {//10->A,11->B
	if (nr <= 9 && nr >= 0) return (char)(nr + '0');//a=10,B=11...
	if (nr <= 35 && nr >= 10) return (char)(nr - 10 + 'A');
		return '0';
}
char* in_baza_doi(int nr)
{
	int dim = log2(nr) + 1;
	char* temp = new char[dim];
	char i = 0;
	while (nr>0) {
		temp[i++] = string_din_nr(nr%2);
		nr /= 2;
	}
	char* rezultat = new char[dim];
	rezultat[i] = '\0';
	for (int j = 0; j<i; j++)
		rezultat[j] = temp[i - j - 1];
	return rezultat;
}

bool versha1(const char* sursa) { //pt fisiserul xorat, verifica pt fiecare linie daca exista cheia de decript (sha1 etc)
	
	string linie;
	char* bits_128 = new char[128];
	char* bits_160 = new char[160];
	char* ki = new char[128];
	char* sha1keyi = new char[160];
	char* rez = new char[160];
	ifstream s(sursa);
	s.seekg(0, ios::beg);
	char* b1c = new char[129];
	while (!s.eof()) {
		getline(s, linie);
		if (linie == "") {
			s.close(); break;
		}
		ki = strcpy(bits_128, linie.substr(0, 128).c_str());
		ki[129] = '\0';
		sha1keyi = strcpy(bits_160, linie.substr(128, 160).c_str());
		SHA1 sha1;
		rez = sha1.ValoareBin(ki, 128);
		if (strncmp(sha1keyi, rez, 160) == 0) {
			ofstream out("B1coleg.txt");
			out << ki;
			out.close();
			s.close();
			//cout << ki << endl;
			return 1;
		}
	}
	s.close();
	return 0;
}

void decriptFile(int i) //xoreaza msg2 cu potentiala cheie a colegului lui Bob
{
	string linie;
	ifstream s("msg2.txt");
	ofstream d("msg2xored.txt");
	int dim =log2(i) + 1;
	char* coleg = new char[dim];
	coleg = in_baza_doi(i);
	//cout << coleg;

	while (!s.eof())
	{	
		getline(s, linie);
		int contor = 0;
		for (int i = 0; i < 288; i++)
		{
			d << ((int)(linie[i] - '0') ^ (int)(coleg[contor]-'0'));
			contor++;
			if (contor == dim) contor = 0;
			
		}
		d << endl;
	}
	s.close();
	d.close();
}

bitset<128> cb() { //pune cheia de decriptare a fisierului intr-un bitset
	bitset<128> col;
	string linie;
	ifstream s("B1coleg.txt");
	getline(s, linie);
	int j = 0;
	for (int i = 127; i >= 0; i--) {
		col[i] = (int)(linie[j] - '0'); //trans char in nr
		j++;
	}
	return col;
}

void decriptPic(char const* sursa, char const* destinatie, bitset<128> ch)
{
	ifstream s(sursa, ios::binary | ios::ate); //ate=at end
	ofstream d(destinatie, ios::binary);
	ifstream::pos_type pst = s.tellg(); //aflam dim
	s.seekg(0, ios::beg); //readuce pointerul la inceput

	char * reznou = new char[pst]; //vect cu fis initial
	s.read(reznou, pst);

	char * bn = new char[pst];
	memset(bn, 0, pst); //ne asig ca toti bitii sunt 0
	int i = 127;
	for(int p=0; p<pst; p++)
	{
		for(int k=7; k>=0; k--)
		{
			bn[p] |= ((((reznou[p] >> k) & 1) ^ ch[i--]) << k); //xor intre cheie si fis binar, pus in bn in seturi de cate 8 biti
			if (i == -1) i = 127;
		}
	}
	d.write(bn, pst);
	s.close();
	s.close();
}

int main()
{
	for (int i=1; i<2903; i++)
	{
		decriptFile(i);
		versha1("msg2xored.txt"); //=> cheia de decriptare a pozei
		if (versha1("msg2xored.txt")) { //=> cheia de comunicare cu serverul
			//cout << in_baza_doi(i);
			ofstream out("cheieptserver.txt");
			out << in_baza_doi(i);
		}
	}
	bitset<128> col;
	col = cb();
	//cout << col;
	decriptPic("cript_F2.jpg", "F2.jpg", col);

	system("pause");
}
