// compile with std=c++11 -O3 -Wall
// Updated 2020-02-14 to stop looking for clusters if the minor allele
// frequency difference is too big; results in substantial speedup.
// Don't forget -O3, it makes the program run much faster.
#include<algorithm>
#include<cmath>
#include<fstream>
#include<iostream>
#include<sstream>
#include<string>
#include<utility>
#include<vector>

#include<math.h>
#include<stdlib.h>

using std::cerr;
using std::cout;
using std::endl;
using std::getline;
using std::ifstream;
using std::pair;
using std::sort;
using std::string;
using std::stringstream;
using std::vector;

struct Record {
    float maf;
    unsigned pos;
    string chr;
    vector<float> data;
    string rs;
    unsigned nChrs;
    long assignment; 
    bool seed;

    Record():assignment(-1){}
};

// To sort by MAF in descending order
bool recordRevMafSort(const Record &a, const Record &b) {
    return (a.maf > b.maf);
}

// Heterozygotes are coded as 0.5 in the input for this script -
// may not produce the ideal result (I haven't thought much
// about it) - this is really written for haploids.
void getFreqs(const Record &x, const Record &y,
        float &pa, float &pb, float &pab) {
    unsigned ab = 0;
    unsigned a = 0;
    unsigned b = 0;
    unsigned c = 0;
    for(unsigned i=0; i < x.data.size(); i++) {
        if(std::isnan(x.data[i]) || std::isnan(y.data[i])) {
            continue;
        }
        c++;
        if(x.data[i] > 0) a++;
        if(y.data[i] > 0) b++;
        if(x.data[i] > 0 && y.data[i] > 0) ab++;
    }
    pa  = ((float)a)/c;
    pb  = ((float)b)/c;
    pab = ((float)ab)/c;
}

// Returns R2 or -1 if it is not possible for R2 to be >= thresh.
// Formula for determining possible R2 comes from Wray 2005.
float calcR2(const Record &x, const Record &y, float thresh=0,
        bool check_maf=true) {
    float pa, pb, pab;
    if(check_maf) {
        if(y.maf < ( (thresh * x.maf) / (1 - x.maf * (1 - thresh)) ) ||
                y.maf >  ( x.maf / (x.maf * (1 - thresh) + thresh) )) {
            return -1.0;
        }
    }
    getFreqs(x, y, pa, pb, pab);
    return pow((pab - (pa * pb)), 2) / (pa * (1-pa) * pb * (1-pb));
}


inline void assignGroups(vector<Record> &data, unsigned max_genotyped, float threshold) {
    // First group using the more highly genotyped markers as seeds.
    // Calculate min and max possible minor allele frequency that would
    // allow a match for each seed, and use this to quickly skip variants
    // that would not match.

    float r2;
    float min_maf;
    float max_maf;
    long assign = 0;

    // loop over number genotyped (for use as seeds)
    for(unsigned g = max_genotyped; g > 0; g--) {

        // loop over variants; reject those that are not genotyped at the
        // current desired level; when one is found that passes and is
        // not assigned to a group, choose it as a seed.
        for(unsigned j = 0; j < data.size(); j++) {
            if(data[j].assignment != -1) continue;
            if(data[j].nChrs < g) continue;
            assign++;
            data[j].assignment = assign;
            data[j].seed = true;
            min_maf = (threshold * data[j].maf) / (1 - data[j].maf * (1 - threshold));
            max_maf = data[j].maf / (data[j].maf * (1 - threshold) + threshold);

            // loop over all variants to potentially group with current seed
            // seed = data[j], test = data[k]
            for(unsigned k = 0; k < data.size(); k++) {
                if(k == j) continue;
                if(data[k].assignment != -1) continue;
                // If the data are sorted by MAF in descending order,
                // then either the top two checks together or the third
                // check by itself seems to give the same answer. Of course,
                // if the data are not sorted, only the third check is
                // appropriate. Doesn't seem to make much difference
                // in run time.
                if(data[k].maf < min_maf) break;    // if sorted descending by MAF
                if(data[k].maf > max_maf) continue; // if sorted by MAF
                //if(data[k].maf < min_maf || data[k].maf > max_maf) continue; // if not sorted by maf
                r2 = calcR2(data[j], data[k], threshold, false);
                if(r2 >= threshold) {
                    data[k].assignment = assign;
                    data[k].seed = false;
                }
            }
        }
    }
    if(data[data.size()-1].assignment == -1) {
        assign++;
        data[data.size()-1].assignment = assign;
        data[data.size()-1].seed = true;
    }
}



int main(int argc, char* argv[]) {
    if(argc < 3) {
        cerr << "Error: wrong number of arguments" << endl;
        cerr << "r2_groups <r2 threshold> [use rs] <input file>" << endl;
        return 1;
    }

    bool use_rs = false;
    if(argc == 4) {
        use_rs = (bool) atoi(argv[2]);
    }

    float threshold = atof(argv[1]);

    vector<Record> data;
    ifstream instream(argv[argc-1]);
    string line;
    float aa;
    bool multiAllele;
    unsigned nChrs;

    unsigned i = 0;
    unsigned max_genotyped = 0;
    while(std::getline(instream, line)) {
        multiAllele = false;
        nChrs = 0;
        if(i != 0) {
            Record rec;
            rec.rs = "";
            unsigned c = 0;
            float refCount = 0.;
            string column;
            stringstream lineStream(line);
            bool na = false;
            while(getline(lineStream, column, '\t')) {
                na = false;
                if(c == 0) {
                    rec.chr = column;
                } else if(c == 1) {
                    rec.pos = atoi(column.c_str());
                } else if(c == 2 && use_rs) {
                    rec.rs = column;
                } else {
                    if(column == "NA" || column == "NaN" || column == "nan") {
                        na = true;
                        aa = std::nanf("");
                    } else {
                        aa = floor(atof(column.c_str()));
                    }
                    if(!na && aa > 1) {
                        multiAllele = true;
                        break;
                    }
                    if(!na) {
                        nChrs++;
                        refCount += aa;
                    }
                    rec.data.push_back(aa);
                }
                c++;
            }
            if(multiAllele) {
                continue;
            } else {
                // No variation at this site
                if(refCount == 0.0 || refCount == nChrs)
                    continue;
                rec.maf = refCount / nChrs;
                if(rec.maf > 0.5)
                    rec.maf = 1.0 - rec.maf;
                rec.nChrs = nChrs;
                data.push_back(rec);
                if(nChrs > max_genotyped)
                    max_genotyped = nChrs;
            }
        }
        i++;
    }

    sort(data.begin(), data.end(), recordRevMafSort);

    // Calculate R2 and group the variants
    assignGroups(data, max_genotyped, threshold);

    // Write the output
    if(use_rs) {
        cout << "chrom\tpos\trs\tmaf\tn\tgroup\tseed" << endl;
        for(unsigned j = 0; j < data.size(); j++) {
            cout << data[j].chr << "\t" << data[j].pos << "\t" <<
                data[j].rs << "\t" <<
                data[j].maf << "\t" << data[j].nChrs << "\t" <<
                data[j].assignment << "\t" << 
                (int)(data[j].seed) << endl;
        }
    } else {
        cout << "chrom\tpos\tmaf\tn\tgroup\tseed" << endl;
        for(unsigned j = 0; j < data.size(); j++) {
            cout << data[j].chr << "\t" << data[j].pos << "\t" <<
                data[j].maf << "\t" << data[j].nChrs << "\t" <<
                data[j].assignment << "\t" <<
                (int)(data[j].seed) << endl;
        }
    }
}
