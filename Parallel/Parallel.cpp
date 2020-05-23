#include "pch.h"
#include <iostream>
#include <string>
#include <fstream>
#include <chrono>
#include <omp.h>
#include <mpi.h>
#define PrimeBase 31

using namespace std;

long long* powers;

void simpleRealization(string text) {
	long long hash = 0;
	long long len = text.length();
	auto start = chrono::system_clock::now();
	for (int i = 0; i < len; ++i)
	{
		hash += (text[i] - 'a' + 1) * powers[i];
	}
	auto end = chrono::system_clock::now();
	chrono::duration<double> seconds = end - start;
	double secondsCount = seconds.count();
	cout << "SIMPLE: hash is " << hash << ", seconds is " << secondsCount << ", text length is " << len << "\n";
}

void openMPRealization(string text, int threadNum) {
	long long hash = 0;
	long long len = text.length();
	auto start = omp_get_wtime();
	#pragma omp parallel for num_threads(threadNum) reduction(+: hash) 
	for (int i = 0; i < len; ++i)
	{
		int code = (text[i] - 'a' + 1);
		hash += code * powers[i];
	}
	auto end = omp_get_wtime();
	double seconds = end - start;
	cout << "OpenMP: hash is " << hash << ", seconds is " << seconds << ", text length is " << len << "\n";
}

void mpiRealiation(int* argc, char** argv[], string text) {
	int size, rank;
	long long hash = 0, tot;
	long long len = text.length();
	MPI_Init(argc, argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
	auto start = chrono::system_clock::now();
	for (int i = 0; i < len; ++i)
	{
		int code = (text[i] - 'a' + 1);
		hash += code * powers[i];
	}
	auto end = chrono::system_clock::now();
	MPI_Reduce(&hash, &tot, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		chrono::duration<double> seconds = end - start;
		double secondsCount = seconds.count();
		cout << "MPI: hash is " << hash << ", seconds is " << secondsCount << ", text length is " << len << "\n";
	}
	MPI_Finalize();
}


string initText() {
	long long powerPrime = 1;
	ifstream inFile;
	inFile.open("C:\\Users\\Lenovo\\source\\repos\\Parallel\\Parallel\\TextSmall.txt");
	if (!inFile) {
		cerr << "Unable to open file Text.txt";
		exit(1);
	}
	string text = "", t;
	while (inFile >> t)
		text += t;
	inFile.close();
	powers = new long long[text.length()];
	for (int i = 0; i < text.length(); i++) {
		powers[i] = powerPrime;
		powerPrime *= PrimeBase;
	}
	return text;
}

int main(int argc, char* argv[])
{
	string text = initText();
	// simple realization 
	simpleRealization(text);
	// OpenMP realization
	openMPRealization(text, 1);
	// MPI realization
	mpiRealiation(&argc, &argv, text);
}

