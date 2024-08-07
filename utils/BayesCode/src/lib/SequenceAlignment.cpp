#include "SequenceAlignment.hpp"
#include <fstream>
#include "BiologicalSequences.hpp"
#include "Random.hpp"
#include "StringStreamUtils.hpp"
using namespace std;

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//     SequenceAlignment
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

int Int(string s) { return atoi(s.c_str()); }

double Double(string s) { return atof(s.c_str()); }

vector<double> SequenceAlignment::GetEmpiricalFreq() const {
    vector<double> in(GetNstate(), 0);
    int n = 0;
    for (int i = 0; i < GetNtaxa(); i++) {
        for (int j = 0; j < GetNsite(); j++) {
            if (GetState(i, j) != unknown) {
                in[GetState(i, j)]++;
                n++;
            }
        }
    }
    for (int i = 0; i < GetNstate(); i++) { in[i] /= n; }
    return in;
}

void SequenceAlignment::ToStream(ostream &os) const {
    os << Ntaxa << '\t' << GetPrintNsite() << '\n';
    for (int i = 0; i < Ntaxa; i++) {
        os << taxset->GetTaxon(i) << '\t';
        for (int j = 0; j < Nsite; j++) { os << statespace->GetState(GetState(i, j)); }
        os << '\n';
    }
}

void SequenceAlignment::ToFasta(ostream &os) const {
    for (int i = 0; i < Ntaxa; i++) {
        os << '>' << taxset->GetTaxon(i) << '\n';
        for (int j = 0; j < GetPrintNsite(); j++) { os << statespace->GetState(GetState(i, j)); }
        os << '\n';
    }
}

int SequenceAlignment::GetNonMissingTriplet() const {
    int nsite = 0;
    for (int j = 2; j < Nsite - 3; j += 3) {
        if (NoMissingColumn(j - 1) && NoMissingColumn(j) && NoMissingColumn(j + 1)) { nsite++; }
    }
    return nsite;
}

void SequenceAlignment::ToStreamTriplet(ostream &os) const {
    int nsite = 0;
    for (int j = 2; j < Nsite - 3; j += 3) {
        if (NoMissingColumn(j - 1) && NoMissingColumn(j) && NoMissingColumn(j + 1)) { nsite++; }
    }

    os << Ntaxa << '\t' << 3 * nsite << '\n';
    int max = 0;
    for (int i = 0; i < Ntaxa; i++) {
        int l = taxset->GetTaxon(i).length();
        if (max < l) { max = l; }
    }

    for (int i = 0; i < Ntaxa; i++) {
        os << taxset->GetTaxon(i);
        for (unsigned int j = 0; j < 5 + max - taxset->GetTaxon(i).length(); j++) { os << ' '; }
        for (int j = 2; j < Nsite - 3; j += 3) {
            if (NoMissingColumn(j - 1) && NoMissingColumn(j) && NoMissingColumn(j + 1)) {
                os << statespace->GetState(GetState(i, j - 1));
                os << statespace->GetState(GetState(i, j));
                os << statespace->GetState(GetState(i, j + 1));
            }
        }
        os << '\n';
    }
    os << '\n';
}

FileSequenceAlignment::FileSequenceAlignment(string filename) { ReadDataFromFile(filename, 0); }

FileSequenceAlignment::FileSequenceAlignment(std::istream &is) {}

int FileSequenceAlignment::ReadDataFromFile(string filespec, int forceinterleaved) {
    string tmp;
    ifstream is(filespec.c_str());
    if (!is) {
        cerr << "error : cannot find data file " << filespec << '\n';
        cerr << "\n";
        exit(1);
    }
    is >> tmp;
    try {
        if (tmp == "#NEXUS") {
            ReadNexus(filespec);
            return 1;
        }
        if (tmp == "#SPECIALALPHABET") {
            ReadSpecial(filespec);
            return 1;
        } else {
            // cerr << "-- [SequenceAlignment] Alignment file uses Phylip format" <<
            // endl;
            if (forceinterleaved == 0) {
                int returnvalue = TestPhylipSequential(filespec);
                if (returnvalue != 0) {
                    // cerr << "-- [SequenceAlignment] Alignment file is sequential" <<
                    // endl;
                    ReadPhylipSequential(filespec);
                    return 1;
                }
            }
            // cerr << "-- [SequenceAlignment] Alignment file is interleaved";
            int returnvalue = TestPhylip(filespec, 1);
            if (returnvalue != 0) {
                // cerr << ", taxon names repeated" << endl;
                ReadPhylip(filespec, 1);
                return 1;
            }
            // cerr << ", taxon names not repeated" << endl;
            TestPhylip(filespec, 0);
            ReadPhylip(filespec, 0);
            return 1;
        }
    } catch (...) {
        exit(1);
        return 0;
    }
    return 1;
}

int FileSequenceAlignment::ReadNexus(string filespec) {
    ifstream theStream(filespec.c_str());
    try {
        GoPastNextWord(theStream, "dimensions");
        GoPastNext(theStream, '=');
        theStream >> Ntaxa;
        GoPastNext(theStream, '=');
        theStream >> Nsite;
        GoPastNextWord(theStream, "format");
        GoPastNextWord(theStream, "datatype");
        GoPastNext(theStream, '=');
        string type;
        theStream >> type;

        if (EquivalentStrings(type, "protein") != 0) {
            statespace = new ProteinStateSpace();
        } else if (EquivalentStrings(type, "dna") != 0) {
            statespace = new DNAStateSpace();
        } else if (EquivalentStrings(type, "rna") != 0) {
            statespace = new RNAStateSpace();
        } else {
            cerr << "error cannot recognise data type\n";
            cerr << type << "\n";
            exit(1);
        }

        Data.assign(Ntaxa, std::vector<int>(Nsite, 0));
        std::vector<std::string> SpeciesNames(Ntaxa, "");

        GoPastNextWord(theStream, "Matrix");

        int l = 0;
        while (l < Nsite) {
            int m = 0;
            for (int i = 0; i < Ntaxa; i++) {
                string temp;
                theStream >> temp;
                while (temp == "[") {
                    unsigned char c;
                    c = 'i';
                    while (c != ']') { c = theStream.get(); }
                    theStream >> temp;
                }

                if (l == 0) {
                    SpeciesNames[i] = temp;
                } else {
                    if (temp != SpeciesNames[i]) {
                        cerr << "error when reading tree base: " << temp << '\t' << SpeciesNames[i]
                             << '\n';
                        exit(1);
                    }
                }

                unsigned char c;
                int k = l;
                do {
                    c = theStream.get();
                    if (c == '[') {
                        while (c != ']') { c = theStream.get(); }
                        c = theStream.get();
                    }
                    if ((c != ' ') && (c != '\t') && (c != '\n') && (c != 13)) {
                        if (c == '(') {
                            Data[i][k] = unknown;
                            while (c != ')') { theStream >> c; }
                        } else if (c == '{') {
                            Data[i][k] = unknown;
                            while (c != '}') { theStream >> c; }
                        } else {
                            ostringstream s;
                            s << c;
                            Data[i][k] = statespace->GetState(s.str());
                        }
                        k++;
                    }
                } while ((!theStream.eof()) && (c != '\n') && (c != 13));
                if (theStream.eof()) {
                    if (i < Ntaxa - 1) {
                        cerr << "error : found " << i << " taxa instead of " << Ntaxa
                             << " in datafile\n";
                        exit(1);
                    }
                }
                if (m == 0) {
                    m = k;
                } else {
                    if (m != k) {
                        cerr << "error when reading nexus : " << m << '\t' << k << '\n';
                        cerr << "taxa : " << i << '\t' << SpeciesNames[i] << '\n';
                        if (m > k) {
                            while (k != m) {
                                Data[i][k] = unknown;
                                k++;
                            }
                        }
                    }
                }
            }
            l = m;
        }
        taxset = new TaxonSet(SpeciesNames);
    } catch (...) {
        cerr << "error while reading data file\n";
        return 0;
    }
    return 1;
}

// ---------------------------------------------------------------------------
//     ReadSpecial()
// ---------------------------------------------------------------------------

int FileSequenceAlignment::ReadSpecial(string filename) {
    int returnvalue = 0;
    try {
        ifstream theStream(filename.c_str());

        string tmp;
        theStream >> tmp;
        theStream >> Ntaxa;
        theStream >> Nsite;
        theStream >> tmp;
        // cerr << tmp << '\n';
        int Nstate = tmp.length();
        auto Alphabet = new char[Nstate];
        int NAlphabetSet = Nstate + 5;
        auto AlphabetSet = new char[NAlphabetSet];
        // cerr << "alphabet size : " << Nstate << '\n';
        // cerr << "alphabet : ";
        for (int i = 0; i < Nstate; i++) {
            Alphabet[i] = tmp[i];
            AlphabetSet[i] = tmp[i];
            // cerr << Alphabet[i] << ' ';
        }
        // cerr << '\n';

        AlphabetSet[Nstate] = '?';
        AlphabetSet[Nstate + 1] = '-';
        AlphabetSet[Nstate + 2] = '*';
        AlphabetSet[Nstate + 3] = 'X';
        AlphabetSet[Nstate + 4] = 'x';

        statespace = new GenericStateSpace(Nstate, Alphabet, NAlphabetSet, AlphabetSet);

        delete[] Alphabet;
        delete[] AlphabetSet;

        Data.assign(Ntaxa, std::vector<int>(Nsite, 0));
        std::vector<std::string> SpeciesNames(Ntaxa, "");

        int ntaxa = 0;
        string temp;
        while ((!theStream.eof()) && (ntaxa < Ntaxa)) {
            theStream >> temp;
            SpeciesNames[ntaxa] = temp;
            int nsite = 0;

            char c;
            do {
                c = theStream.get();
                if ((!theStream.eof()) && (c != ' ') && (c != '\n') && (c != '\t') && (c != 13)) {
                    if (c == '(') {
                        Data[ntaxa][nsite] = unknown;
                        while (c != ')') { theStream >> c; }
                    } else if (c == '{') {
                        Data[ntaxa][nsite] = unknown;
                        while (c != '}') { theStream >> c; }
                    } else {
                        int p = 0;
                        while ((p < NAlphabetSet) && (c != AlphabetSet[p])) { p++; }
                        if (p == NAlphabetSet) {
                            cout << "error: does not recognise character. taxon " << ntaxa << '\t'
                                 << SpeciesNames[ntaxa] << "  site  " << nsite << '\t' << c << '\n';
                            exit(1);
                        }
                        if (p >= Nstate) {
                            Data[ntaxa][nsite] = unknown;
                        } else {
                            for (int l = 0; l < Nstate; l++) {
                                if (c == Alphabet[l]) { Data[ntaxa][nsite] = l; }
                            }
                        }
                    }
                    nsite++;
                }
            } while ((!theStream.eof()) && (nsite < Nsite));
            ntaxa++;
        }
        taxset = new TaxonSet(SpeciesNames);
    } catch (...) {
        cerr << "error while reading data file\n";
        return 0;
    }
    return returnvalue;
}

// ---------------------------------------------------------------------------
//     ReadPhylip()
// ---------------------------------------------------------------------------

int FileSequenceAlignment::TestPhylipSequential(string filespec) {
    ifstream theStream(filespec.c_str());
    try {
        string temp;
        theStream >> temp;
        if (IsInt(temp) == 0) {
            cerr << "error when reading data\n";
            cerr << "data should be formatted as follows:\n";
            cerr << "#taxa #sites\n";
            cerr << "name1 seq1.....\n";
            cerr << "name2 seq2.....\n";
            cerr << "...\n";
            cerr << '\n';
            exit(1);
        }
        Ntaxa = atoi(temp.c_str());
        theStream >> temp;
        if (IsInt(temp) == 0) {
            cerr << "error when reading data\n";
            cerr << "data should be formatted as follows:\n";
            cerr << "#taxa #sites\n";
            cerr << "name1 seq1.....\n";
            cerr << "name2 seq2.....\n";
            cerr << "...\n";
            cerr << '\n';
            exit(1);
        }
        Nsite = atoi(temp.c_str());

        std::vector<std::string> SpeciesNames(Ntaxa, "");

        int AAcomp = 1;
        int DNAcomp = 1;
        int RNAcomp = 1;

        int ntaxa = 0;
        while ((!theStream.eof()) && (ntaxa < Ntaxa)) {
            theStream >> temp;
            SpeciesNames[ntaxa] = temp;
            int nsite = 0;

            char c = ' ';
            do {
                c = theStream.get();
                if ((!theStream.eof()) && (c != ' ') && (c != '\n') && (c != '\t') && (c != 13)) {
                    if (c == '(') {
                        while (c != ')') { theStream >> c; }
                    } else if (c == '{') {
                        while (c != '}') { theStream >> c; }
                    } else {
                        int p = 0;
                        if (DNAcomp != 0) {
                            while ((p < DNAN) && (c != DNAset[p])) { p++; }
                            if (p == DNAN) { DNAcomp = 0; }
                        }
                        p = 0;
                        if (RNAcomp != 0) {
                            while ((p < RNAN) && (c != RNAset[p])) { p++; }
                            if (p == RNAN) { RNAcomp = 0; }
                        }
                        p = 0;
                        if (AAcomp != 0) {
                            while ((p < AAN) && (c != AAset[p])) { p++; }
                            if (p == AAN) { AAcomp = 0; }
                        }
                    }
                    nsite++;
                }
            } while ((!theStream.eof()) && (nsite < Nsite));
            if (theStream.eof()) {
                if (nsite < Nsite) { return 0; }
            }
            ntaxa++;
        }
        if (theStream.eof()) {
            if (ntaxa < Ntaxa) { return 0; }
        }
        if (DNAcomp != 0) {
            statespace = new DNAStateSpace;
        } else if (RNAcomp != 0) {
            statespace = new RNAStateSpace;
        } else if (AAcomp != 0) {
            statespace = new ProteinStateSpace;
        } else {
            return 0;
        }
    } catch (...) { return 0; }
    return 1;
}

void FileSequenceAlignment::ReadPhylipSequential(string filespec) {
    ifstream theStream(filespec.c_str());
    ReadPhylipSequentialFromStream(theStream);
}

void FileSequenceAlignment::ReadPhylipSequentialFromStream(istream &theStream) {
    try {
        string temp;
        theStream >> temp;
        if (IsInt(temp) == 0) {
            cerr << "error when reading data\n";
            cerr << "data should be formatted as follows:\n";
            cerr << "#taxa #sites\n";
            cerr << "name1 seq1.....\n";
            cerr << "name2 seq2.....\n";
            cerr << "...\n";
            cerr << '\n';
            exit(1);
        }
        Ntaxa = Int(temp);
        theStream >> temp;
        if (IsInt(temp) == 0) {
            cerr << "error when reading data\n";
            cerr << "data should be formatted as follows:\n";
            cerr << "#taxa #sites\n";
            cerr << "name1 seq1.....\n";
            cerr << "name2 seq2.....\n";
            cerr << "...\n";
            cerr << '\n';
            exit(1);
        }
        Nsite = Int(temp);

        Data.assign(Ntaxa, std::vector<int>(Nsite, 0));
        std::vector<std::string> SpeciesNames(Ntaxa, "");

        int ntaxa = 0;
        while ((!theStream.eof()) && (ntaxa < Ntaxa)) {
            theStream >> temp;
            SpeciesNames[ntaxa] = temp;
            int nsite = 0;

            char c;
            do {
                c = theStream.get();
                if ((!theStream.eof()) && (c != ' ') && (c != '\n') && (c != '\t') && (c != 13)) {
                    if (c == '(') {
                        Data[ntaxa][nsite] = unknown;
                        while (c != ')') { theStream >> c; }
                    } else if (c == '{') {
                        Data[ntaxa][nsite] = unknown;
                        while (c != '}') { theStream >> c; }
                    } else {
                        ostringstream s;
                        s << c;
                        Data[ntaxa][nsite] = statespace->GetState(s.str());
                    }
                    nsite++;
                }
            } while ((!theStream.eof()) && (nsite < Nsite));
            ntaxa++;
        }
        taxset = new TaxonSet(SpeciesNames);
    } catch (...) {
        cerr << "error while reading data file\n";
        exit(1);
    }
}

int FileSequenceAlignment::TestPhylip(string filespec, int repeattaxa) {
    ifstream theStream(filespec.c_str());
    try {
        string temp;
        theStream >> temp;
        if (IsInt(temp) == 0) {
            cerr << "error when reading data\n";
            cerr << "data should be formatted as follows:\n";
            cerr << "#taxa #sites\n";
            cerr << "name1 seq1.....\n";
            cerr << "name2 seq2.....\n";
            cerr << "...\n";
            cerr << '\n';
            exit(1);
        }
        Ntaxa = Int(temp);
        theStream >> temp;
        if (IsInt(temp) == 0) {
            cerr << "error when reading data\n";
            cerr << "data should be formatted as follows:\n";
            cerr << "#taxa #sites\n";
            cerr << "name1 seq1.....\n";
            cerr << "name2 seq2.....\n";
            cerr << "...\n";
            cerr << '\n';
            exit(1);
        }
        Nsite = Int(temp);

        std::vector<std::string> SpeciesNames(Ntaxa, "");

        int AAcomp = 1;
        int DNAcomp = 1;
        int RNAcomp = 1;

        int l = 0;
        int block = 0;
        while (l < Nsite) {
            block++;
            int m = 0;
            for (int i = 0; i < Ntaxa; i++) {
                if ((l == 0) || (repeattaxa != 0)) {
                    string temp;
                    theStream >> temp;
                    if (l == 0) {
                        SpeciesNames[i] = temp;
                    } else {
                        if (temp != SpeciesNames[i]) { return 0; }
                    }
                }

                unsigned char c;
                int k = l;
                do {
                    c = theStream.get();
                    if ((!theStream.eof()) && (c != ' ') && (c != '\n') && (c != '\t') &&
                        (c != 13)) {
                        if (c == '(') {
                            while (c != ')') { theStream >> c; }
                        } else if (c == '{') {
                            while (c != '}') { theStream >> c; }
                        } else {
                            int p = 0;
                            if (DNAcomp != 0) {
                                while ((p < DNAN) && (c != DNAset[p])) { p++; }
                                if (p == DNAN) { DNAcomp = 0; }
                            }
                            p = 0;
                            if (RNAcomp != 0) {
                                while ((p < RNAN) && (c != RNAset[p])) { p++; }
                                if (p == RNAN) { RNAcomp = 0; }
                            }
                            p = 0;
                            if (AAcomp != 0) {
                                while ((p < AAN) && (c != AAset[p])) { p++; }
                                if (p == AAN) { AAcomp = 0; }
                            }
                        }
                        k++;
                    }
                } while ((!theStream.eof()) && (c != '\n') && (c != 13) && (c != 10));
                if (theStream.eof()) {
                    if (i < Ntaxa - 1) {
                        cerr << "error : found " << i << " taxa instead of " << Ntaxa
                             << " in datafile\n";
                        exit(1);
                    }
                }
                c = theStream.peek();
                while ((!theStream.eof()) && ((c == '\n') || (c == 13))) {
                    c = theStream.get();
                    c = theStream.peek();
                }
                if (m == 0) {
                    m = k;
                } else {
                    if (m != k) {
                        cerr << "in test phylip\n";
                        cerr << "error when reading data non matching number of sequences "
                                "in block "
                                "number "
                             << block << " for taxon " << i + 1 << " " << SpeciesNames[i] << '\n';
                        cerr << "taxa : " << i << '\t' << SpeciesNames[i] << '\n';
                        cerr << "read " << k << " instead of " << m << "characters\n";
                        exit(1);
                    }
                }
            }
            l = m;
        }
        if (l < Nsite) {
            cerr << "error : reached end of stream \n";
            cerr << "data should be formatted as follows:\n";
            cerr << "#taxa #sites\n";
            cerr << "name1 seq1.....\n";
            cerr << "name2 seq2.....\n";
            cerr << "...\n";
            cerr << '\n';
            exit(1);
        }
        if (DNAcomp != 0) {
            statespace = new DNAStateSpace;
            // cerr << "dna sequences\n";
        } else if (RNAcomp != 0) {
            statespace = new DNAStateSpace;
            // cerr << "rna sequences\n";
        } else if (AAcomp != 0) {
            statespace = new ProteinStateSpace;
            // cerr << "protein sequences\n";
        } else {
            // cerr << "format not recognised\n";
            return 0;
        }
    } catch (...) {
        cerr << "error while reading data file\n";
        return 0;
    }
    return 1;
}

void FileSequenceAlignment::ReadPhylip(string filespec, int repeattaxa) {
    ifstream theStream(filespec.c_str());
    try {
        string temp;
        theStream >> temp;
        if (IsInt(temp) == 0) {
            cerr << "error when reading data\n";
            cerr << "data should be formatted as follows:\n";
            cerr << "#taxa #sites\n";
            cerr << "name1 seq1.....\n";
            cerr << "name2 seq2.....\n";
            cerr << "...\n";
            cerr << '\n';
            exit(1);
        }
        Ntaxa = Int(temp);
        theStream >> temp;
        if (IsInt(temp) == 0) {
            cerr << "error when reading data\n";
            cerr << "data should be formatted as follows:\n";
            cerr << "#taxa #sites\n";
            cerr << "name1 seq1.....\n";
            cerr << "name2 seq2.....\n";
            cerr << "...\n";
            cerr << '\n';
            exit(1);
        }
        Nsite = Int(temp);
        // cerr << Ntaxa << '\t' << Nsite << '\n';

        Data.assign(Ntaxa, std::vector<int>(Nsite, 0));
        std::vector<std::string> SpeciesNames(Ntaxa, "");

        int l = 0;
        int block = 0;
        while (l < Nsite) {
            block++;
            int m = 0;
            for (int i = 0; i < Ntaxa; i++) {
                if ((l == 0) || (repeattaxa != 0)) {
                    string temp;
                    theStream >> temp;
                    if (l == 0) {
                        SpeciesNames[i] = temp;
                    } else {
                        if (temp != SpeciesNames[i]) {
                            cerr << "error when reading data: read " << temp << " instead of "
                                 << SpeciesNames[i] << '\n';
                            exit(1);
                        }
                    }
                }

                unsigned char c;
                int k = l;
                do {
                    c = theStream.get();
                    if ((!theStream.eof()) && (c != ' ') && (c != '\n') && (c != '\t') &&
                        (c != 13)) {
                        if (c == '(') {
                            Data[i][k] = unknown;
                            while (c != ')') { theStream >> c; }
                        } else if (c == '{') {
                            Data[i][k] = unknown;
                            while (c != '}') { theStream >> c; }
                        } else {
                            ostringstream s;
                            s << c;
                            Data[i][k] = statespace->GetState(s.str());
                        }
                        k++;
                    }
                } while ((!theStream.eof()) && (c != '\n') && (c != 13));
                if (theStream.eof()) {
                    if (i < Ntaxa - 1) {
                        cerr << "error : found " << i << " taxa instead of " << Ntaxa
                             << " in datafile\n";
                        exit(1);
                    }
                }
                c = theStream.peek();
                while ((!theStream.eof()) && ((c == '\n') || (c == 13))) {
                    c = theStream.get();
                    c = theStream.peek();
                }

                if (m == 0) {
                    m = k;
                } else {
                    if (m != k) {
                        cerr << "error when reading data non matching number of sequences "
                                "in block "
                                "number "
                             << block << " for taxon " << i << " " << SpeciesNames[i] << '\n';
                        cerr << "taxa : " << i << '\t' << SpeciesNames[i] << '\n';
                        cerr << "read " << k << " instead of " << m << "characters\n";
                        exit(1);
                    }
                }
            }
            l = m;
        }
        if (l < Nsite) {
            cerr << "error : reached end of stream \n";
            cerr << "data should be formatted as follows:\n";
            cerr << "#taxa #sites\n";
            cerr << "name1 seq1.....\n";
            cerr << "name2 seq2.....\n";
            cerr << "...\n";
            cerr << '\n';
            exit(1);
        }
        taxset = new TaxonSet(SpeciesNames);
    } catch (...) { cerr << "error while reading data file\n"; }
}
