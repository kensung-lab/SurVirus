#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <cassert>

using namespace std;

int main(int argc, char* argv[]) {

	ifstream reads_in(argv[1]), kmers_in(argv[2]);

	vector<string> kmers;
	string s;
	while (kmers_in >> s) {
		kmers.push_back(s);
	}	

	while (reads_in >> s) {
		bool found = false;
		for (string& kmer : kmers) {
//			for (int i = 0; i < s.length(); i++) if (s[i] == 'N') s[i] = 'A';
			if (strstr(s.data(), kmer.data()) != NULL) {
				found = true;
				break;
			}
		}
		if (!found) cout << s << endl;
		assert(found);
	}

}
