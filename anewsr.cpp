// cl anewsr.cpp /std:c++17
// anewsr.exe -f sample_data/sample_genome.fa -g sample_data/sample_genome.gff -o sample_data/sample_genome_output.csv --cas9 -v

#include <iostream>
#include <cstring>
#include <fstream>
#include <string>
#include <map>
#include <unordered_map>
#include <regex.h>
#include <sstream>
#include "argparse.hpp"

using namespace std;

unordered_map<string, int> antiswap{
{"A",0},
{"T",1},
{"U",1},
{"C",2},
{"G",3},
};

char complement(char n)
{
    switch (n)
    {
    case 'A':
        return 'U';
    case 'T':
        return 'A';
    case 'G':
        return 'C';
    case 'C':
        return 'G';
    }
    return ' ';
}

double matrix_multiply(double** a, int** b, int c, int d) {

    double trace = 0;
    double** res = new double* [c];

    for (int i = 0; i < c; i++) {
        res[i] = new double[c];

        for (int j = 0; j < c; j++) {
            res[i][j] = 0;
            double thr = 5;
            for (int k = 0; k < d; k++) {
                if(a[i][k] > -thr && a[i][k] < thr && b[k][j] > -thr && b[k][j] < thr)
                res[i][j] += a[i][k] * b[k][j];
                }
              //reject garbage values that might have slipped in        
        }
        trace += res[i][i];
    }
    return trace;
}

vector<int> find_pam_site2(string target, string sequence){
    vector<int> sites;
    if (target.length() > sequence.length())
        return sites;
    else 
        for (int i = 0;i < sequence.length() - target.length();i++) {
            
            string s = sequence.substr(i, target.length());            
            if (s == target)
                sites.push_back(i);
        }
    return sites;
}

string get_reverse_complement(string sequence) {

    unordered_map<int, string> swap{
  {0, "A"},
  {1, "T"},
  {1, "U"},
  {2, "C"},
  {3, "G"},

    };
    string comp = "";
    for (int i = 0;i < sequence.length();i++) {
        string s = { sequence[i] };
       
        comp += swap[antiswap[s]+pow(-1,(antiswap[s]%2))];
      
    }
    return comp;
}



double rs1_score(string sequence) {
    int** x = new int* [sequence.length()];
    for (int i = 0; i < sequence.length(); i++) {
        x[i] = new int[4];  
        string delta = { sequence[i] };
        x[i][antiswap[delta]] = 1;
    }

    int** y = new int* [sequence.length()];
    for (int i = 0; i < sequence.length()-1; i++) {
        y[i] = new int[16];
        string theta = { sequence[i]};
        string thetwo = { sequence[i + 1] };
        y[i][4*antiswap[theta]+antiswap[thetwo]] = 1;    
    } 

    string first_order[] = { "G02", "A03", "C03", "C04", "C05",
        "G05", "A06", "C06", "C07", "G07",
        "A12", "A15", "C15", "A16", "C16",
        "T16", "A17", "G17", "C18", "G18",
        "A19", "C19", "G20", "T20", "G21",
        "T21", "C22", "T22", "T23", "C24",
        "G24", "T24", "A25", "C25", "T25",
        "G28", "T28", "C29", "G30" };
   
    double first_scores[] = { -0.2753771, -0.3238875, 0.17212887, -0.1006662, -0.2018029,
    0.24595663, 0.03644004, 0.09837684, -0.7411813, -0.3932644,
    -0.466099, 0.08537695, -0.013814, 0.27262051, 0.1190226,
    -0.2859442, 0.09745459, -0.1755462, -0.3457955, -0.6780964,
    0.22508903, -0.5077941, -0.4173736, -0.054307, 0.37989937,
    -0.0907126, 0.05782332, -0.5305673, -0.8770074, -0.8762358,
    0.27891626, -0.4031022, -0.0773007, 0.28793562, -0.2216372,
    -0.6890167, 0.11787758, -0.1604453, 0.38634258 };
   
    string second_order[] = { "GT02", "GC05", "AA06", "TA06", "GG07",
        "GG12", "TA12", "TC12", "TT12", "GG13",
        "GA14", "GC14", "TG17", "GG19", "TC19",
        "CC20", "TG20", "AC21", "CG21", "GA21",
        "GG21", "TC22", "CG23", "CT23", "AA24",
        "AG24", "AG25", "CG25", "TG25", "GT27",
        "GG29" };

    double second_scores[] = { -0.6257787, 0.30004332, -0.8348362, 0.76062777, -0.4908167,
    -1.5169074, 0.7092612, 0.49629861, -0.5868739, -0.3345637,
    0.76384993, -0.5370252, -0.7981461, -0.6668087, 0.35318325,
    0.74807209, -0.3672668, 0.56820913, 0.32907207, -0.8364568,
    -0.7822076, -1.029693, 0.85619782, -0.4632077, -0.5794924,
    0.64907554, -0.0773007, 0.28793562, -0.2216372, 0.11787758,
    -0.69774 };   

    double intersect = 0.59763615;
    double low_gc = -0.2026259;
    double high_gc = -0.1665878;

    double** first_matrix = new double* [4];
    for (int i = 0;i < 4;i++) {
        first_matrix[i] = new double[30];
    }

    for (int i = 0;i <39;i++) {
       
        string delta = first_order[i].substr(0, 1);
        int val = stoi(first_order[i].substr(1))-1;
        first_matrix[antiswap[delta]][val] = first_scores[i];
    }

    //first order ends above

    double** second_matrix = new double* [16];
    for (int i = 0;i < 16;i++) {
        second_matrix[i] = new double[29];
    }

    for (int i = 0;i < 31;i++) {

        string delta = second_order[i].substr(0,1) ;
        string deltwo = second_order[i].substr(1,1);
        int val = stoi(second_order[i].substr(2)) - 1;
        second_matrix[4 * antiswap[delta] + antiswap[deltwo]][val] = second_scores[i];
       
    }
    //second order ends above

    double score_first = matrix_multiply(first_matrix, x, 4, 30);
    double score_second = matrix_multiply(second_matrix, y, 16, 29);
    double gc_count = 0;
    for (int l = 0;l < sequence.length();l++) {
        if (sequence[l] == 'G' || sequence[l] == 'C')
            gc_count++;
    }

    if (gc_count < 10)
        gc_count = low_gc;
    else gc_count = high_gc;
    
    double score = ((1 / (1 + exp(-(intersect + gc_count + score_first + score_second)))));      
    return score;

}

void import_fasta_file(int argc, string argv, vector<string>& chromosome_fasta, vector<string>& sequence_fasta) {
    ifstream fasta;
    fasta.open(argv);

    if (!fasta)
        fprintf(stderr, "Error opening FASTA file '%s'\n", argv);
    

    else {
        string line_buf, seq;

        while (getline(fasta, line_buf)) {

            if (line_buf.at(0) == '>') {
                chromosome_fasta.push_back(line_buf.substr(1));
                if (!seq.empty()) {
                    sequence_fasta.push_back(seq);
                    seq.clear();
                }
            }
            else {
                string temp(line_buf);
                seq = seq + temp;
            }
        }
        sequence_fasta.push_back(seq);

        fasta.close();
        printf("Genome file successfully imported and formatted\n");

    }
}


//verbatim as in .c
void import_gff_file(int argc, string argv, vector<string>& chromosome_gff, vector<string>& source_gff,
    vector<string>& feature_gff, vector<string>& start_gff, vector<string>& end_gff,
    vector<string>& score_gff, vector<string>& strand_gff, vector<string>& phase_gff,
    vector<string>& attributes_gff) {

    ifstream gff;
    gff.open(argv);
    string line, word;
    string find("##");
    int i = 0;
    if (!gff)
        fprintf(stderr, "Error opening GFF file. '%s'\n");
    

    else {
        printf("Annotation file successfully imported. \n");
        while (getline(gff, line)) {

            if (line.find(find) == std::string::npos) {

                i = 0;

                stringstream s(line);

                while (getline(s, word, '\t')) {
                    if (i == 0) chromosome_gff.push_back(word);
                    else if (i == 1) source_gff.push_back(word);
                    else if (i == 2) feature_gff.push_back(word);
                    else if (i == 3) start_gff.push_back(word);
                    else if (i == 4) end_gff.push_back(word);
                    else if (i == 5) score_gff.push_back(word);
                    else if (i == 6) strand_gff.push_back(word);
                    else if (i == 7) phase_gff.push_back(word);
                    else if (i == 8) attributes_gff.push_back(word);
                    i += 1;
                }
            }

        }
        gff.close();
        printf("Annotation database successfully generated \n");
    }
}

string correct_substr(string input, int start, int end) {
    string output = "";
    for (int i = start;i < end;i++)
        output += input.at(i);
    return output;
}

string get_id(char type, string s) {
    string out = "";
    out = out + type + correct_substr(s,3, s.length());
        int randlimit = 7;
    static const char alphanum[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    

    for (int i = 0; i < randlimit; i++) {
        out += alphanum[rand() % (sizeof(alphanum) - 1)];
    }
    return out;
    
}

string get_gRNA_sequence(string input_sequence)
{
    std::transform(input_sequence.begin(), input_sequence.end(), input_sequence.begin(), ::toupper);
    transform(
        begin(input_sequence),
        end(input_sequence),
        begin(input_sequence),
        complement);
    reverse(input_sequence.begin(), input_sequence.end());
    return input_sequence;
}

int main(int argc, char** argv) {
   // test();
    vector<string> chromosome_fasta, sequence_fasta, chromosome_gff, source_gff, feature_gff,
        start_gff, end_gff, score_gff, strand_gff, phase_gff, attributes_gff;
    vector<int> pam_sites_postive, pam_sites_negative;
    
    argparse::ArgumentParser parser("");

    parser.add_argument("-f", "--fasta")
        .help(" path to input file in FASTA format")
        .required();

    parser.add_argument("-g", "--gff")
        .help(" path to input file in GFF format");

    parser.add_argument("-p", "--phytozome")
        .help(" path to input annotation info file in TXT format");

    parser.add_argument("-o", "--output")
        .help(" path to output file")
        .default_value("data.csv");

    parser.add_argument("-l", "--length")
        .help(" length of the gRNA sequence. default = 20")
        .default_value(20)
        .scan<'i', int>();

    parser.add_argument("-L", "--flanking")
        .help(" length of flanking region for verification, default = 200")
        .default_value(200)
        .scan<'i', int>();

    parser.add_argument("--cas9")
        .help(" specifies that design will be made for the Cas9 CRISPR system")
        .default_value(false)
        .implicit_value(true);

    parser.add_argument("--cpf1")
        .help(" specifies that design will be made for the Cpf1 CRISPR system")
        .default_value(false)
        .implicit_value(true);

    parser.add_argument("--CUDA")
        .help(" runs processing steps utilizing the GPU instead of CPU where possible")
        .default_value(false)
        .implicit_value(true);

    parser.add_argument("-v", "--verbose")
        .help(" prints visual indicators for each iteration")
        .default_value(false)
        .implicit_value(true);  

    try {
        parser.parse_args(argc, argv);
    }
    catch (const std::runtime_error & err) {
        std::cerr << err.what() << std::endl;
        std::cerr << parser;
        std::exit(1);
    }

    if (parser.get("-o").find(".csv") == string::npos) {
        parser.get("-o") = "data.csv";
    }
    ofstream myfile(parser.get("-o"));
    if (!myfile.is_open())
    {
        cout << "\nThe output file cannot be written to. Please try again.\n";
    }
    else {
        myfile << "crispr_id, crispr_sys, sequence, long_sequence, chromosome, start_pos, end_pos, cutsite, strand, on_site_score, features\n";
        

        import_fasta_file(argc, parser.get("-f"), chromosome_fasta, sequence_fasta);
        import_gff_file(argc, parser.get("-g"), chromosome_gff, source_gff, feature_gff, start_gff,
            end_gff, score_gff, strand_gff, phase_gff, attributes_gff);
        
        auto input = parser.get("-f");
        auto in2 = parser.get("-g");
        string crispr_guide;
        vector<string> complete_dataset;
        for (int i = 0;i < chromosome_fasta.size();i++) {
            for (string sequence : sequence_fasta) {

                pam_sites_postive.clear();
                if (parser["--cas9"] == true) {
                    regex re("(?=.GG)");
                    
                    pam_sites_postive = find_pam_site2("GG", sequence);
                    pam_sites_postive[0]--;
                    pam_sites_postive[1]--;
                    
                    for (int j : pam_sites_postive) {
                       
                        int pam_location[] = { j - (parser.get<int>("-l") + 1),j - 1 };
                        
                        if (pam_location[0] >= 5 && pam_location[0] + 5 <= sequence.length() + 10 && pam_location[1] >= 5 && pam_location[1] <= sequence.length() + 10) {
                            string shortseq = get_gRNA_sequence(correct_substr(sequence, pam_location[0], pam_location[1]));
                            string longseq = get_gRNA_sequence(correct_substr(sequence, pam_location[0]-5, pam_location[1]+5));
                         
                            crispr_guide = get_id('A', chromosome_fasta[i])+", cas9, " + shortseq + ", " + longseq + ", " + chromosome_fasta[i] + ", " + to_string(pam_location[0]) + ", " +
                                to_string(pam_location[1]) + ", " + to_string(pam_location[1] - 3) + ", +, "+to_string( rs1_score(longseq))+", completed";
                            myfile << crispr_guide + "\n";
                          
                        }
                    }
                  
                    pam_sites_postive = find_pam_site2("CC", sequence);
                    pam_sites_postive[0]--;
                    pam_sites_postive[1]--;
                    for (int j : pam_sites_postive) {
                        int pam_location[] = { j + 1,(j + (parser.get<int>("-l") + 1)) };
                        if (pam_location[0] >= 5 && pam_location[0] + 5 <= sequence.length() + 10 && pam_location[1] >= 5 && pam_location[1] <= sequence.length() + 10) {
                            string shortseq = get_gRNA_sequence(correct_substr(sequence, pam_location[0], pam_location[1]));
                            string longseq = get_gRNA_sequence(correct_substr(sequence, pam_location[0] - 5, pam_location[1] + 5));
                            crispr_guide = get_id('A', chromosome_fasta[i])+", cas9, " + shortseq + ", " + longseq + ", " + chromosome_fasta[i] + ", " + to_string(pam_location[0]) + ", " +
                                to_string(pam_location[1]) + ", " + to_string(pam_location[1] - 3) + ", -, " + to_string(rs1_score(longseq)) + ", completed";
                            myfile << crispr_guide + "\n";

                       
                        }
                    }

                }//cas9 ends above, cpf1 starts below
                if (parser["--cpf1"] == true) {
                    regex rg("(?=TTT.)");
                    pam_sites_postive = find_pam_site2("TTT", sequence);
                    pam_sites_postive[0]--;
                    pam_sites_postive[1]--;

                    for (int j : pam_sites_postive) {
                        int pam_location[] = { j + 4,(j + (parser.get<int>("-l") + 4)) };
                        if (pam_location[0] >= 0 && pam_location[0] + 5 <= sequence.length() && pam_location[1] >= 05 && pam_location[1] <= sequence.length()) {
                            string shortseq = get_gRNA_sequence(correct_substr(sequence, pam_location[0], pam_location[1]));
                            string longseq = get_gRNA_sequence(correct_substr(sequence, pam_location[0] - 5, pam_location[1] + 5));
                            crispr_guide = get_id('B', chromosome_fasta[i])+", cpf1, " + shortseq + ", " + longseq + ", " + chromosome_fasta[i] + ", " + to_string(pam_location[0]) + ", " +
                                to_string(pam_location[1]) + ", " + to_string(pam_location[1] - 3) + ", +, " + to_string(rs1_score(longseq)) + ", completed";
                            myfile << crispr_guide + "\n";
                           
                        }
                    }
                   
                    pam_sites_postive = find_pam_site2("AAA", sequence);
                    pam_sites_postive[0]--;
                    pam_sites_postive[1]--;

                    for (int j : pam_sites_postive) {
                        
                        int pam_location[] = { j - (parser.get<int>("-l") + 1),(parser.get<int>("-l") - 1) };
                        if (pam_location[0] >= 0 && pam_location[0] + 5 <= sequence.length() && pam_location[1] >= 05 && pam_location[1] <= sequence.length()) {
                            string shortseq = get_gRNA_sequence(sequence.substr(pam_location[0], pam_location[1]));
                            string longseq = get_gRNA_sequence(sequence.substr(pam_location[0] - 5, pam_location[1] + 5));
                            crispr_guide = get_id('B', chromosome_fasta[i])+", cpf1, " + shortseq + ", " + longseq + ", " + chromosome_fasta[i] + ", " + to_string(pam_location[0]) + ", " +
                                to_string(pam_location[1]) + ", " + to_string(pam_location[1] - 3) + ", -, " + to_string(rs1_score(longseq)) + ", completed";
                            myfile << crispr_guide + "\n";                           
                        }
                    }

                }
            }
        }
        printf("Output file successfully generated. \n");
    } 
	
	return 100;
}

void test() {
    string test = "UCCGUGUGCGUACGUAAAAUCAGUAUACAA";
    cout << "Start Test";
    cout << get_reverse_complement(test) + "\n";
    rs1_score(test);
    cout << get_id('A', "Chr01") + "\n";
    vector<int> fpm = find_pam_site2("AA", test);
    cout << fpm.at(0);
    cout << "End Test";
}