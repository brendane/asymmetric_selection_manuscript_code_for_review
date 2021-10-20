#include<iostream>
#include<set>
#include<string>
#include<vector>

#include <ctype.h>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

using std::cout;
using std::endl;
using std::set;
using std::string;
using std::vector;

const char* chr = "chrUn";

int main(int argc, char ** argv) {
	gzFile fp;
    if(string(argv[1]) == "-") {
        fp = gzdopen(fileno(stdin), "r");
    } else {
        fp = gzopen(argv[1], "r");
    }

    vector<string> seqs;
    vector<string> names;
    kseq_t *seq = kseq_init(fp);
	while (true) {
        if(kseq_read(seq) < 0) break;
        seqs.push_back(string(seq->seq.s));
        names.push_back(string(seq->name.s));
    }

    int length = seqs[0].size();
    cout << "##fileformat=VCFv4.2" << endl;
    cout << "##contig=<ID=" << chr << ",length=" << length << ">" << endl;
    cout << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
    cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for(size_t i=0; i<seqs.size(); i++) {
        cout << "\t" << names[i];
    }
    cout << endl;
    vector<char> alleles;
    vector<char> alleleCodes;
    for(int j=0; j<length; j++) {
        alleles.clear();
        alleleCodes.clear();
        for(size_t i=0; i<seqs.size(); i++) {
            char a = toupper(seqs[i][j]);
            alleles.push_back(a);
            if(a == 'A' || a == 'G' || a == 'C' || a == 'T') {
                bool found = false;
                for(size_t k=0; k < alleleCodes.size(); k++) {
                    if(alleleCodes[k] == a) {
                        found = true;
                        break;
                    }
                }
                if(!found) {
                    alleleCodes.push_back(a);
                }
            }
        }

        if(alleleCodes.size() > 0) {
            cout << chr << "\t" << j+1 << "\t.\t" << alleleCodes[0] << "\t";
            if(alleleCodes.size() > 1) {
                for(size_t k=1; k<alleleCodes.size(); k++) {
                    if(k > 1) cout << ",";
                    cout << alleleCodes[k];
                }
            } else {
                cout << "<*>";
            }
            cout << "\t.\t.\t.\tGT";

            for(size_t i=0; i<seqs.size(); i++) {
                char a = alleles[i];
                bool found = false;
                int kk = 0;
                for(size_t k=0; k<alleleCodes.size(); k++) {
                    if(alleleCodes[k] == a) {
                        found = true;
                        kk = k;
                        break;
                    }
                }
                if(!found) {
                    cout << "\t.";
                } else {
                    cout << "\t" << kk;
                }
            }
            cout << endl;
        }
    }


    kseq_destroy(seq);
	gzclose(fp);
	return 0;
}
