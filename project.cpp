#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <fstream>
#include <unordered_map>
#include <regex>
#include <stdexcept>
#include <utility>
#include <chrono>



#include "spoa/spoa.hpp"

using namespace std;
namespace fs = filesystem;

size_t getCurrentRSS() {
#if defined(__linux__)
    ifstream status_file("/proc/self/status");
    string line;
    while (getline(status_file, line)) {
        if (line.rfind("VmRSS:", 0) == 0) {
            istringstream iss(line);
            string key, value, unit;
            iss >> key >> value >> unit;
            return stoul(value); // KB
        }
    }
#endif
    return 0;
}


struct FastqRead {
  string id;
  string sequence;
  string plus;
  string quality;
};


/*Function for parsing FASTQ files.
  It takes filename and desired length of sequences, parses whole file, and only returns vector<FastqRead> only containing sequences of 
  desired length. 
  EV*/
vector<FastqRead> parseFile(const string& filename, size_t targetLength = 0){
  vector<FastqRead> reads;
  ifstream file(filename);

  if(!file) {
    cerr << "Error opening file: " << filename << endl;
    return reads;
  }


  int n_longer = 0;
  string line;
  while (getline(file, line)){
    FastqRead read;
    read.id = line;
    getline(file, read.sequence);
    getline(file, read.plus);
    getline(file, read.quality);

    if(read.sequence.length() == targetLength or targetLength == 0) {
      reads.push_back(read);
    }
    else {
      n_longer++;
    }
  }
  
  return reads;

}

/*Function for parsing ground truth files. Maps ground_truth_sequence_id to ground_truth_sequence for later use
EV*/
unordered_map<string, string> parseGT_File(const string& filename, size_t targetLength = 249) {
    unordered_map<string, string> id_to_sequence;

    ifstream file(filename);
    if (!file) {
        cerr << " Error opening GT file: " << filename << "\n";
        return id_to_sequence;
    }

    cout << "Parsing GT file: " << filename << "\n";

    string line, current_id;
    while (getline(file, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            current_id = line.substr(1);
            auto space_pos = current_id.find(' ');
            if (space_pos != string::npos) {
                current_id = current_id.substr(0, space_pos);  
            }
        } else {
            if (!current_id.empty() && line.length() == targetLength) {
                id_to_sequence[current_id] = line;
            }
        }
    }

    return id_to_sequence;
}

/*Function for getting filenames of sample reads. It takes path to directory with all fastq samples as parameter and returns vector<path> that
  includes all desired file paths
  EV*/
vector<fs::path> iterateDirectory(const string& dirPath){
    vector<fs::path> files;
    for(const auto& entry : fs::directory_iterator(dirPath)){
      if(entry.path().extension() == ".fastq" || entry.path().extension() == ".fasta" ){//take into consideration only fastq or fasta files
        if(entry.path().filename().string().compare(0, 1, "J") == 0) //take into consideration only files that starts with "J", indicating deer samples (can be changed for different purposes)
          files.push_back(entry.path());
      }
    }
    return files;
  }


//LM
//use spoa library to generate graph from which consensus and MSA are calculated
spoa::Graph generateGraph(const vector<string>& reads){
  
  cout << "generating graph..." << endl;
  auto alignmentEngine = spoa::AlignmentEngine::Create(
      spoa::AlignmentType::kNW, 0, -1, -100);
  spoa::Graph graph{};

  for (const auto& read : reads) {
    auto alignment = alignmentEngine->Align(read, graph);
    graph.AddAlignment(alignment, read);
  }
  cout << "graph generated :D" << endl;
  return graph;
}


//LM
//returns hamming distance between 2 sequences
size_t hammingDistance(const string& a, const string& b) {

    if (a.size() != b.size()) {
        throw invalid_argument("Sequences must have equal length for Hamming distance");
    }
    size_t dist = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) ++dist;
    }
    return dist;
}

//LM
//calculates consensus sequence (used as centroid) of given sequences 
string computeCentroid(const vector<string>& sequences) {

    if (sequences.empty()) return "";

    auto alignment_engine = spoa::AlignmentEngine::Create(
        spoa::AlignmentType::kNW, 3, -1, -5);  // Adjusted parameters
    spoa::Graph graph{};

    for (const auto& seq : sequences) {
        auto alignment = alignment_engine->Align(seq, graph);
        graph.AddAlignment(alignment, seq);
    }
    return graph.GenerateConsensus();
}

//LM
//iterates over readings and generates clusters

//k = number of missmatches tolerated for one cluster
//s = min nuber of cluster members
vector<vector<string> > clusteringAlgorithm(const vector<string> &sequences, size_t k, size_t minClusterSize) {
  
  vector<vector<string> > clusters;

  for (const string &seq: sequences) {
    bool added = false;

    // Try adding to an existing cluster
    for (auto &cluster: clusters) {
      // Compare to the cluster's representative (first sequence)
      if (hammingDistance(seq, cluster[0]) < k) {
        cluster.push_back(seq);
        added = true;
        break;
      }
    }

    // If not similar to any existing cluster, create new cluster
    if (!added) {
      clusters.push_back({seq});
    }
  }

  // Filter out small clusters
  vector<vector<string> > filtered;
  for (const auto& cluster : clusters){
    if (cluster.size() > minClusterSize)
      filtered.push_back(cluster);
  }
  
  return filtered;
}

//LM
//tests appropriate hamming distance used for calculating clusters
void cluster_parameter_test(const vector<string> & readingsList) {

  for (size_t k = 0; k < 50; ++k) {
    vector<vector<string>> clusters = clusteringAlgorithm(readingsList, k, 0);
      cout<<k<<","<<clusters.size()<<endl;
      
  }
}


//LM
//returns map where key is sample ID and value is a vector of corresponding sequence string
unordered_map<string, vector<string>> get_readingsList(const vector<fs::path> & filesToParse, int targetLength = 0) {
  vector<string> chosenSequences;
  unordered_map<std::string, vector<string>> readingsList;

  for (const auto& filePath : filesToParse) {
    string filename = filePath.filename().string();
    string sample_id = filename.substr(0, filename.find("_"));

    vector<FastqRead> reads = parseFile(filePath.string(), targetLength);//filter for deer sequences of length 296
      vector<string> chosenSequences;
      for (const auto& read : reads) {
      chosenSequences.push_back(read.sequence);
    }
    readingsList[sample_id]=chosenSequences;
    chosenSequences.clear();
  }
    cout<<"Files parsed."<<endl;

  return readingsList;
}

//LM
//function which generates corresponding .msa files - replaces every sequence in sample with globally aligned sequence - does this for every sample
void generateMSA_Files(const unordered_map<string, vector<string>> & readingsList) {

    // Go up two levels from current directory to get to project root
    fs::path project_root = fs::current_path().parent_path();
    fs::path output_dir = project_root /"data"/ "msa_files";

    // Create output directory if it doesn't exist
    bool createdNewly = fs::create_directories(output_dir);
    if (!createdNewly && !fs::is_empty(output_dir)) {
        return;
    }
    cout << "Output directory path: " << output_dir << endl;
    cout << "Directory created :) "<<endl;

    for (auto const& reading: readingsList) {
        string sample_id = reading.first;
        auto reads = reading.second;
        auto graph = generateGraph(reads);

        //use spoa to get consensus and msa
        auto consensus = graph.GenerateConsensus();
        auto msa = graph.GenerateMultipleSequenceAlignment();

        fs::path output_path = output_dir / (sample_id + ".msa");

        // Write MSA to file
        ofstream outfile(output_path);
        if (!outfile.is_open()) {
            cerr << "Error opening file: " << output_path << endl;
            continue;
        }
        // Write each sequence in the MSA to separate line
        for (const auto& aligned_seq : msa) {
            outfile << aligned_seq << '\n';
        }
        outfile.close();
        cout << "Created MSA file: " << sample_id << ".msa" << " with " << msa.size() << " sequences\n";
    }
}

//LM
//creates clusters list for each sample (takes data from msa file)
unordered_map<string, vector<vector<string>>>
 create_clusters(size_t maxDist = 10, size_t minClusterSize = 5) {
    unordered_map<string, vector<vector<string>>>
 clustered_samples;


    fs::path project_root = fs::current_path().parent_path();
    fs::path input_dir = project_root /"data"/ "msa_files";
    if (!fs::exists(input_dir)) {
        cerr << "Output directory does not exist.\n";
        return clustered_samples;
    }

    for (const auto& entry : fs::directory_iterator(input_dir)) {//iterates over each directory
        if (entry.path().extension() == ".msa") {
            ifstream infile(entry.path());
            if (!infile.is_open()) {
                cerr << "Failed to open " << entry.path() << "\n";
                continue;
            }

            vector<string> sequences;
            string line;
            while (getline(infile, line)) {
                if (!line.empty())
                    sequences.push_back(line);
            }

            infile.close();

            // Cluster the sequences
            auto clusters = clusteringAlgorithm(sequences, maxDist, minClusterSize);

            // Get sample ID from filename
            string sample_id = entry.path().stem().string(); // e.g., "J30"

            clustered_samples[sample_id] = clusters;
        }
    }

    return clustered_samples;
}


//function for parsing .fasta file that contains centroids of all clusters. Returns vector<string> containing all sequences from given file
//EV
vector<string> parse_fasta_sequences(const string& fasta_path) {
    vector<string> sequences;
    ifstream file(fasta_path);
    if (!file.is_open()) {
        throw runtime_error("Cannot open FASTA file: " + fasta_path);
    }

    string line, current_seq;
    while (getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!current_seq.empty()) {
                sequences.push_back(current_seq);
                current_seq.clear();
            }
        } else {
            current_seq += line;
        }
    }
    if (!current_seq.empty()) {
        sequences.push_back(current_seq);
    }

    return sequences;
}


/*creates directory for saving FASTA files. For each sample file this function calculates consensus sequence for all of the clusters,
 and saves them as FASTA files  
 EV*/
void writeCentroidsToFasta(const unordered_map<string, vector<vector<string>>>& cluster_list) {
    fs::path project_root = fs::current_path().parent_path();
    fs::path fasta_output_dir = project_root / "data" / "consensus_files";
    fs::path clusters_output_dir = project_root / "data" / "clusters_txt";

    //create directories if they dont exist
    bool created = fs::create_directories(fasta_output_dir);
    fs::create_directories(clusters_output_dir);
    if (!created) {
        return;
    }

    for (const auto& [sample_id, clusters] : cluster_list) {
        // open .fasta file to get centroid sequences
        ofstream fasta_out(fasta_output_dir / (sample_id + ".fasta"));
        if (!fasta_out.is_open()) {
            cerr << "Cannot write to FASTA file for sample " << sample_id << endl;
            continue;
        }

        // open .txt file to write all sequences of each cluster
        ofstream cluster_out(clusters_output_dir / (sample_id + ".txt"));
        if (!cluster_out.is_open()) {
            cerr << "Cannot write to cluster TXT file for sample " << sample_id << endl;
            continue;
        }

        for (size_t i = 0; i < clusters.size(); ++i) {
            const auto& cluster = clusters[i];
            string centroid = computeCentroid(cluster);

            // write centroid to .fasta
            fasta_out << ">cluster" << i + 1 << "\n" << centroid << "\n";

            // write all sequences of each cluster to .txt
            cluster_out << "cluster" << i + 1 << "\n";
            for (const auto& seq : cluster) {
                cluster_out << seq << "\n";
            }
            cluster_out << "\n";
        }

        fasta_out.close();
        cluster_out.close();
        cout << "Wrote centroids and cluster sequences for " << sample_id << "\n";
    }
}

//implementation of Smith-Waterman algorithm
//EV
int smithWaterman(const string& seq1, const string& seq2, int match = 2, int mismatch = -1, int gap = -1) {
    int m = seq1.size();
    int n = seq2.size();
    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));

    int maxScore = 0;

    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int score = (seq1[i - 1] == seq2[j - 1]) ? match : mismatch;
            dp[i][j] = max({
                0,
                dp[i - 1][j - 1] + score,
                dp[i - 1][j] + gap,
                dp[i][j - 1] + gap
            });
            maxScore = max(maxScore, dp[i][j]);
        }
    }

    return maxScore;
}

//compare each cluster centroid with every ground truth seqence using SmithWaterman algorithm. Write out best match with score
//EV
void compareCentroidsWithGT(const vector<string>& centroids, const unordered_map<string, string>& gt_sequences) {
    int idx=0;
    for (const auto& centroid : centroids) {
        int bestScore = numeric_limits<int>::min();
        string bestGT_ID;

        for (const auto& [gt_id, gt_seq] : gt_sequences) {
            int score = smithWaterman(centroid, gt_seq);
            if (score > bestScore) {
                bestScore = score;
                bestGT_ID = gt_id;
            }
        }

        cout << "Cluster" << idx++ << " matched best with GT ID: " 
                  << bestGT_ID << " (score: " << bestScore << "/498)\n"; //498 score is maximum that sw algorithm can give in our case
    }
}

//parse ground truth and cluster files, then compare each cluster consensus to ground truth sequences
//EV
void checkGroundTruth(fs::path path, int index){

    auto groundTruthSequences = parseGT_File(path.string());
    fs::path project_root = fs::current_path().parent_path();
    fs::path directory = project_root /"data"/ "consensus_files";

    vector<string> sequencesToCheck;
    if(index == 29){
      cout << "Checking J29.fasta" << endl;
      sequencesToCheck = parse_fasta_sequences((directory / "J29.fasta").string());
    }
    else if(index == 30){
      cout << "checking J30.fasta" << endl;
      sequencesToCheck = parse_fasta_sequences((directory / "J30.fasta").string());
    }
    
    compareCentroidsWithGT(sequencesToCheck, groundTruthSequences);
}

int main(int argc, char **argv) {
    auto start = chrono::high_resolution_clock::now();
    size_t mem_start = getCurrentRSS();

  if (argc < 3) {
    cerr << "Missing folder path/s!" << endl;
    return 1;
  }
  string sampleFolder = argv[1]; 
  vector<fs::path> filesToParse = iterateDirectory(sampleFolder);

  unordered_map<string, vector<string> > readingsList = get_readingsList(filesToParse, 296);

  generateMSA_Files(readingsList);
    //map with sample id as key and cluster list ad value
    unordered_map<string, vector<vector<string>> > cluster_list = create_clusters(10, 0);
  
  //generate .fasta files for every cluster in each sample file  
  writeCentroidsToFasta(cluster_list);

  //check if we get good results
  string GT_Folder = argv[2];
  vector<fs::path> GT_files = iterateDirectory(GT_Folder);

  int index = 29;
  for (const auto& path : GT_files) {
    checkGroundTruth(path, index++);
  }

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    std::cout << "Execution time: " << duration.count() << " ms" << std::endl;

    size_t mem_end = getCurrentRSS();
    std::cout << "Memory cost: " << (mem_end - mem_start) << " KB" << std::endl;

  return 0;

}
