#include"Bob.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstddef>

//bulinele 1 si 2

bitset<64> KB() {
	bitset<64> bob;
	string linie;
	ifstream b("KB.txt");
	getline(b, linie);
	int j = 0;
	for (int i = 63; i >= 0; i--) {
		bob[i] = (int)(linie[j] - '0'); //trans char in nr
		j++;
	}
	return bob;
}

bitset<128> KF() {
	bitset<128> b1;
	string linie;
	ifstream b("B1.txt");
	getline(b, linie);
	int j = 0;
	for (int i = 127; i >= 0; i--) {
		b1[i] = (int)(linie[j] - '0'); //trans char in nr
		j++;
	}
	return b1;
}

void decriptFile(bitset<64> bob)
{
	string linie;
	ifstream s("msg1.txt");
	ofstream d("msg1xored.txt");
	while (!s.eof())
	{	
		getline(s, linie);
		int contor = 63;
		for (int i = 0; i < 288; i++)
		{
			d<< ((int)(linie[i]-'0') ^ bob[contor]);
			contor--;
			if (contor == -1) contor = 63;
		}
		d << endl;
	}
	s.close();
	d.close();
}

void versha1(){
	string linie;
	char* bits_128 = new char[128];
	char* bits_160 = new char[160];
	char* ki = new char[128];
	char* sha1keyi = new char[160];
	char* rez = new char[160];
	ifstream s("msg1xored.txt");
	ofstream d("B1.txt");
	while (!s.eof()) {	
		getline(s, linie);
		ki=strcpy(bits_128,linie.substr(0, 128).c_str());
		ki[129] = '\0';
		sha1keyi=strcpy(bits_160,linie.substr(128, 160).c_str());
		SHA1 sha1;
		rez = sha1.ValoareBin(ki, 128);
		if (strncmp(sha1keyi, rez, 15) == 0)
			d << ki;
	}
	s.close();
	d.close();
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

void main()
{
	bitset<64> bob;
	bob = KB();
	decriptFile(bob); //Xoram msg1
	versha1(); //Aflam cheia de decriptare a pozei din msg1xored, deducem ca primul fisier ii era adresat lui Bob !

	bitset<128> b1;
	b1 = KF();

	decriptPic("cript_F1.png", "F1.png", b1);

	system("pause");
}