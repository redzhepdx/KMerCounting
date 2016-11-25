#pragma once
#ifndef _MERPROCESS_
#define _MERPROCESS_

#include "BloomFilter.h"
#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>

#define BUFFER 1000
#define KMER_SIZE 20
#define NUMBEROFBITS 718879379
#define NUMBEROFHASHES 17
#define MAX_FREQ_KMER_COUNT 25


using namespace std;

static string ReverseComplement(string str) {

	string newStr = str;

	for (size_t i = 0; i < str.length(); i++) {

		if (str[i] == 'A') {
			newStr.at(str.length() - i - 1) = 'T';
		}
		else if (str[i] == 'T') {
			newStr.at(str.length() - i - 1) = 'A';
		}
		else if (str[i] == 'G') {
			newStr.at(str.length() - i - 1) = 'C';
		}
		else if (str[i] == 'C') {
			newStr.at(str.length() - i - 1) = 'G';
		}
	}

	return newStr;
}

static string LexicographicalComprison(string str_1, string str_2) {

	return str_1.compare(str_2) >= 0 ? str_2 : str_1;
}

static void ProcessInfo(string seq,
	unordered_map<string,
	size_t>  &hashMap,
	BloomFilter &bloom_filter,
	size_t MerSize) {

	for (size_t i = 0; i < seq.length() - MerSize + 1; i++) {

		string substr = seq.substr(i, MerSize);

		substr = LexicographicalComprison(substr, ReverseComplement(substr));

		if (bloom_filter.possiblyContains(substr.c_str(), substr.length())) {
			if (!hashMap.count(substr)) {
				hashMap[substr] = 0;
			}
		}
		else {
			bloom_filter.add(substr.c_str(), substr.length());
		}
	}

}

static void CountReadKmers(string seq,
	unordered_map<string,
	size_t> &hashMap,
	BloomFilter &bloom_filter,
	size_t MerSize) {
	for (size_t i = 0; i < seq.length() - MerSize + 1; i++) {

		string substr = seq.substr(i, MerSize);

		char *substrChar = new char[substr.size() + 1];

		std::strcpy(substrChar, LexicographicalComprison(substr, ReverseComplement(substr)).c_str());

		if (hashMap.count(substrChar)) {
			hashMap[substrChar]++;
		}
	}
}

static void FindMaxFrequencyKMersInHash(vector<string> &mostFrequentMers,
	vector<size_t> &frequencies,
	unordered_map<string, size_t> &hashMap,
	size_t MostFreqMerCount) {

	for (size_t i = 0; i < MostFreqMerCount; i++) {

		size_t max = 0;
		string maxFreqMer;
		for (auto keyVal : hashMap) {
			if (keyVal.second > max && keyVal.second > 1) {
				max = keyVal.second;
				maxFreqMer = keyVal.first;
			}
		}
		hashMap[maxFreqMer] = 0;

		mostFrequentMers.push_back(maxFreqMer);

		frequencies.push_back(max);
	}


}
	

#endif