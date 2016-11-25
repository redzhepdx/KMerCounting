//#include "stdafx.h"
#include "BloomFilter.h"
#include "merprocess.h"
#include "kseq.h"
#include <zlib.h>
#include <iostream>
#include <cmath>
#include <ctime>

KSEQ_INIT(gzFile, gzread);

using namespace std;




int main(int argc,char** argv)
{
	char *filename = new char[1024];
	size_t MerSize;
	size_t MostFreqMerCount;
	strcpy(filename, "ERR266411_1.first1000.fastq");
	MerSize = 30;
	MostFreqMerCount = 25;
	/*
	if (argc != 0) {
		strcpy(filename, argv[1]);
		MerSize = std::atoi(argv[2]);
		MostFreqMerCount = std::atoi(argv[2]);
		std::cout << filename << endl;
		std::cout << MerSize << endl;
		std::cout << MostFreqMerCount << endl;
	}
	else {
		strcpy(filename,"ERR266411_1.first1000.fastq");
		MerSize = 30;
		MostFreqMerCount = 25;
	}
	*/
	std::clock_t startTime = std::clock();
	std::cout << "!!!!!!!!!!!!!!!Start!!!!!!!!!!!!!!!!!!!!!" << endl;

	// your fastq/a file has to be read by zlib
	gzFile fp;

	vector <string> seqs;
	
	unordered_map <string, size_t> T_hashMap;
	
	BloomFilter bloom_filter(NUMBEROFBITS, NUMBEROFHASHES);

	const char* big_input_fastq = "ERR055763.filt.fastq";

	const char* small_input_fastq = "ERR266411_1.first1000.fastq";
	
	fp = gzopen(filename, "r");

	//  initializing fastq/a stream
	kseq_t* seq = kseq_init(fp);

	unsigned counter = 0;

	// read fastq/a
	while (kseq_read(seq) >= 0) {

		seqs.push_back(seq->seq.s);
		
		counter++;
		if (counter == BUFFER) {

			for (size_t i = 0; i < seqs.size(); i++) {
				ProcessInfo(seqs.at(i), T_hashMap, bloom_filter, MerSize);
			}

			for (size_t i = 0; i < seqs.size(); i++) {
				CountReadKmers(seqs.at(i), T_hashMap, bloom_filter , MerSize);
			}
			
			seqs.clear();
			counter = 0;
		}

	}

	vector<string> mostFrequentMers;
	vector<size_t> frequencies;
	
	FindMaxFrequencyKMersInHash(mostFrequentMers, frequencies, T_hashMap, MostFreqMerCount);

	double duration = (std::clock() - startTime) / (double)CLOCKS_PER_SEC;

	std::cout << "duration : " << duration << endl;
	std::cout << "Size of Hash Table : " << T_hashMap.size() << endl;
	std::cout << "Press Enter To See 25 Most Frequent K-MERS" << endl;
	getchar();

	for (size_t i = 0; i < mostFrequentMers.size(); i++) {
		std::cout << mostFrequentMers.at(i) << "-> Frequency : " << frequencies.at(i) << endl;
	}


	std::cout<<"!!!!!!DONE!!!!!!!!";
	getchar();


	
	// clean up
	kseq_destroy(seq);
	gzclose(fp);
	
}
