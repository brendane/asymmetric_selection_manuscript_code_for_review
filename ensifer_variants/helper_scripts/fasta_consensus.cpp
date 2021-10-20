/*
 * Produce a consensus sequence for DNA alignments. Breaks ties at random
 * and does not include gaps and Ns in allele counts.
 *
 * g++ -Wall -O3 -o fasta_consensus -lz fasta_consensus.cpp 
 *
 * Requires kseq.h.
 *
 * fasta_consensus <input.fasta> [<record name>] [<random seed>]
 */

#include<cstdlib>
#include <ctype.h>
#include<iostream>
#include<unordered_map>
#include<string>
#include<stdexcept>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include<vector>

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

using std::cerr;
using std::cout;
using std::endl;
using std::rand;
using std::runtime_error;
using std::string;
using std::unordered_map;
using std::vector;

char getConsensus(unordered_map<char, int>& alleles);

int main(int argc, char ** argv) {
    string seqName = "consensus";
    if(argc > 2) {
        seqName = string(argv[2]);
    }
    srand(time(0));
    if(argc > 3) {
        srand(atoi(argv[3]));
    }

	gzFile fp;
    if(strcmp(argv[1], "-") == 0) {
        fp = gzdopen(fileno(stdin), "r");
    } else {
        fp = gzopen(argv[1], "r");
    }

    vector<string> seqs;
    kseq_t *seq = kseq_init(fp);
    unsigned n = 0;
    unsigned length = 0;
	while (true) {
        if(kseq_read(seq) < 0) break;
        if(n == 0) {
            length = seq->seq.l;
        } else {
            if(seq->seq.l != length) {
                throw runtime_error("sequences are not aligned");
            }
        }
        n++;
        seqs.push_back(string(seq->seq.s));
    }
    kseq_destroy(seq);
	gzclose(fp);

    cout << ">" << seqName << endl;
    unordered_map<char, int> alleles;
    char a;
    for(unsigned j=0; j<length; j++) {
        alleles.clear();
        for(unsigned i=0; i<seqs.size(); i++) {
            a = toupper(seqs[i][j]);
            switch(a) {
                case 'G':
                    alleles['G']++;
                    break;
                case 'C':
                    alleles['C']++;
                    break;
                case 'A':
                    alleles['A']++;
                    break;
                case 'T':
                    alleles['T']++;
                    break;
                case 'U':
                    alleles['U']++;
                    break;
                case 'M':
                    alleles['A']++;
                    alleles['C']++;
                    break;
                case 'R':
                    alleles['A']++;
                    alleles['G']++;
                    break;
                case 'W':
                    alleles['A']++;
                    alleles['T']++;
                    break;
                case 'S':
                    alleles['C']++;
                    alleles['G']++;
                    break;
                case 'Y':
                    alleles['C']++;
                    alleles['T']++;
                    break;
                case 'K':
                    alleles['T']++;
                    alleles['G']++;
                    break;
                case 'V':
                    alleles['A']++;
                    alleles['C']++;
                    alleles['G']++;
                    break;
                case 'H':
                    alleles['A']++;
                    alleles['C']++;
                    alleles['T']++;
                    break;
                case 'D':
                    alleles['A']++;
                    alleles['T']++;
                    alleles['G']++;
                    break;
                case 'B':
                    alleles['C']++;
                    alleles['T']++;
                    alleles['G']++;
                    break;
                default:
                    break;
            }
        }
        cout << getConsensus(alleles);
    }
    cout << endl;

	return 0;
}


char getConsensus(unordered_map<char, int>& alleles) {
    int c = 0;
    char a = 'N';
    vector<char> best;

    for(unordered_map<char, int>::iterator it = alleles.begin(); it != alleles.end(); ++it) {
        if(it->second > c) {
            c = it->second;
        }
    }
    
    if(c > 0) {
        for(unordered_map<char, int>::iterator it = alleles.begin(); it != alleles.end(); ++it) {
            if(it->second == c) {
                best.push_back(it->first);
            }
        }
        if(best.size() == 1) {
            a = best[0];
        } else {
            a = best[rand()/((RAND_MAX + 1u)/ best.size())];
        }
    }

    return a;
}
