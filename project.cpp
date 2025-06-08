#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <fstream>
#include <unordered_map>
#include <regex>
#include <stdexcept>
#include <utility>
#include "spoa/spoa.hpp"

using namespace std;
namespace fs = filesystem;

struct FastqRead {
  string id;
  string sequence;
  string plus;
  string quality;
};


/*Function for parsing FASTQ files.
  It takes filename and desired length of sequences, parses whole file, and only returns vector<FastqRead> only containing sequences of 
  desired length.*/
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
  
  if (reads.empty()) {
    cerr << "Error parsing file: " << filename <<"for target length:"<<targetLength<<endl;
  }
  return reads;

}

/*Function for parsing ground truth files. Maps ground_truth_sequence_id to ground_truth_sequence for later use*/
std::unordered_map<std::string, std::string> parseGT_File(const std::string& filename, size_t targetLength = 249) {
    std::unordered_map<std::string, std::string> id_to_sequence;

    std::ifstream file(filename);
    if (!file) {
        std::cerr << " Error opening GT file: " << filename << "\n";
        return id_to_sequence;
    }

    std::cout << "Parsing GT file: " << filename << "\n";

    std::string line, current_id;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            current_id = line.substr(1);
            auto space_pos = current_id.find(' ');
            if (space_pos != std::string::npos) {
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
  includes all desired file paths*/
vector<fs::path> iterateDirectory(const string& dirPath){
    vector<fs::path> files;
    for(const auto& entry : fs::directory_iterator(dirPath)){
      if(entry.path().extension() == ".fastq" || entry.path().extension() == ".fasta" ){//or entry.path().extension() == ".fasta"){ //take into consideration only fastq or fasta files
        if(entry.path().filename().string().compare(0, 1, "J") == 0) //take into consideration only files that starts with "J", indicating deer samples (can be changed for different purposes)
          files.push_back(entry.path());
      }
    }
    return files;
  }


//use spoa to generate graph from which consensus and MSA are calculated
spoa::Graph generateGraph(const vector<std::string>& reads){
  std::cout << "number of reads: " << reads.size() << std::endl;
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


size_t hammingDistance(const std::string& a, const std::string& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Sequences must have equal length for Hamming distance");
    }
    size_t dist = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) ++dist;
    }
    return dist;
}

//calculates consensus sequence of given sequences
std::string computeCentroid(const std::vector<std::string>& sequences) {
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

std::vector<std::vector<std::string>> clusteringWithCentroid(
    const std::vector<std::string>& sequences,
    size_t maxDist,
    size_t minClusterSize)
{
    if (sequences.empty()) return {};

    // Verify uniform sequence length
    const size_t ref_len = sequences[0].size();
    for (const auto& seq : sequences) {
        if (seq.size() != ref_len) {
            throw std::invalid_argument("All sequences must be of equal length");
        }
    }

    // Initial clustering using first sequence as centroid
    std::vector<std::vector<std::string>> clusters;
    std::vector<std::string> centroids;

    for (const auto& seq : sequences) {
        bool assigned = false;
        for (size_t i = 0; i < centroids.size(); ++i) {
            if (hammingDistance(seq, centroids[i]) <= maxDist) {
                clusters[i].push_back(seq);
                assigned = true;
                break;
            }
        }
        if (!assigned) {
            clusters.push_back({seq});
            centroids.push_back(seq);
        }
    }

    // Batch centroid recomputation
    for (size_t i = 0; i < clusters.size(); ++i) {
        centroids[i] = computeCentroid(clusters[i]);
    }

    // Refinement iterations (max 5)
    const size_t MAX_ITER = 5;
    for (size_t iter = 0; iter < MAX_ITER; ++iter) {
        std::vector<std::vector<std::string>> new_clusters(centroids.size());

        // Assign sequences to closest centroid
        for (const auto& seq : sequences) {
            size_t min_dist = maxDist + 1;
            int best_index = -1;

            for (size_t i = 0; i < centroids.size(); ++i) {
                size_t dist = hammingDistance(seq, centroids[i]);
                if (dist < min_dist) {
                    min_dist = dist;
                    best_index = i;
                }
            }

            if (min_dist <= maxDist) {
                new_clusters[best_index].push_back(seq);
            } else {
                new_clusters.push_back({seq});
            }
        }

        // Remove empty clusters
        std::vector<std::vector<std::string>> non_empty;
        for (auto& cluster : new_clusters) {
            if (!cluster.empty()) {
                non_empty.push_back(std::move(cluster));
            }
        }

        // Update clusters and centroids
        clusters = non_empty;
        centroids.clear();
        for (auto& cluster : clusters) {
            centroids.push_back(computeCentroid(cluster));
        }
    }

    // Filter small clusters
    std::vector<std::vector<std::string>> result;
    for (auto& cluster : clusters) {
        if (cluster.size() >= minClusterSize) {
            result.push_back(std::move(cluster));
        }
    }

    return result;
}



//iterates over reading and generates clusters
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
  //cout << "Unfiltered clusters count: " << clusters.size() << endl;
  for (const auto& cluster : clusters){
    //cout << "  Cluster size: " << cluster.size() << endl;
    if (cluster.size() > minClusterSize)
      filtered.push_back(cluster);
  }
  
  return filtered;
}

void cluster_parameter_test(const vector<std::string> & readingsList) {
  for (size_t k = 0; k < 50; ++k) {
    vector<vector<string>> clusters = clusteringAlgorithm(readingsList, k, 0);
      cout<<k<<","<<clusters.size()<<endl;
      //cout<<"Cluster number for k="<<k<<" -> "<<clusters.size()<<endl;
  }
}

std::unordered_map<std::string, std::vector<std::string>> get_readingsList(const vector<fs::path> & filesToParse, int targetLength = 0) {
  vector<std::string> chosenSequences;
  std::unordered_map<std::string, std::vector<std::string>> readingsList;
  for (const auto& filePath : filesToParse) {
    string filename = filePath.filename().string();
    string sample_id = filename.substr(0, filename.find("_"));
    //cout << "parsing file " << filePath<<" for sample: "<<sample_id<<endl;
    vector<FastqRead> reads = parseFile(filePath.string(), targetLength);//filter for deer sequences of length 296
      vector<std::string> chosenSequences;
      for (const auto& read : reads) {
      chosenSequences.push_back(read.sequence);
    }
    if ( sample_id=="J30") {
      cout<<"J30 sample size: "<<chosenSequences.size()<<endl;
    }
    readingsList[sample_id]=chosenSequences;
    chosenSequences.clear();
  }
    cout<<"Files parsed."<<endl;

  return readingsList;
}
//function which generates corresponeding .msa files - replaces every sequence in sample with globally alligned sequence - does this for every sample
//
void generateMSA_Files(const unordered_map<std::string, std::vector<std::string>> & readingsList) {
    // Go up two levels from current directory to get to project root
    fs::path project_root = fs::current_path().parent_path();
    fs::path output_dir = project_root /"data"/ "msa_files";

    // Create output directory if it doesn't exist
    bool created = fs::create_directories(output_dir);
    if (!created) {
        return;
    }
    std::cout << "Output directory path: " << output_dir << endl;
    std::cout << "Directory created :) "<<endl;

    for (auto const& reading: readingsList) {
        string sample_id = reading.first;
        auto reads = reading.second;
        auto graph = generateGraph(reads);

        //use spoa to get consensus and msa
        auto consensus = graph.GenerateConsensus();
        auto msa = graph.GenerateMultipleSequenceAlignment();

        fs::path output_path = output_dir / (sample_id + ".msa");

        // Write MSA to file
        std::ofstream outfile(output_path);
        if (!outfile.is_open()) {
            std::cerr << "Error opening file: " << output_path << std::endl;
            continue;
        }
        // Write each sequence in the MSA to separate line
        for (const auto& aligned_seq : msa) {
            outfile << aligned_seq << '\n';
        }
        outfile.close();
        std::cout << "Created MSA file: " << sample_id << ".msa" << " with " << msa.size() << " sequences\n";
    }
}

//creates clusters list for each sample (takes data from msa file)
std::unordered_map<std::string, std::vector<std::vector<std::string>>>
 create_clusters(size_t maxDist = 10, size_t minClusterSize = 5) {
    std::unordered_map<std::string, std::vector<std::vector<std::string>>>
 clustered_samples;


    fs::path project_root = fs::current_path().parent_path();
    fs::path input_dir = project_root /"data"/ "msa_files";
    if (!fs::exists(input_dir)) {
        std::cerr << "Output directory does not exist.\n";
        return clustered_samples;
    }

    for (const auto& entry : fs::directory_iterator(input_dir)) {//iterates over each directory
        if (entry.path().extension() == ".msa") {
            std::ifstream infile(entry.path());
            if (!infile.is_open()) {
                std::cerr << "Failed to open " << entry.path() << "\n";
                continue;
            }

            std::vector<std::string> sequences;
            std::string line;
            while (std::getline(infile, line)) {
                if (!line.empty())
                    sequences.push_back(line);
            }

            infile.close();

            // Cluster the sequences
            auto clusters = clusteringAlgorithm(sequences, maxDist, minClusterSize);

            // Get sample ID from filename
            std::string sample_id = entry.path().stem().string(); // e.g., "J30"

            clustered_samples[sample_id] = clusters;
        }
    }

    return clustered_samples;
}


//function for parsing .fasta file that contains centroids of all clusters. Returns vector<string> containing all sequences from given file
std::vector<std::string> parse_fasta_sequences(const std::string& fasta_path) {
    std::vector<std::string> sequences;
    std::ifstream file(fasta_path);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open FASTA file: " + fasta_path);
    }

    std::string line, current_seq;
    while (std::getline(file, line)) {
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
 and saves them as FASTA files  */
void writeCentroidsToFasta(const std::unordered_map<std::string, std::vector<std::vector<std::string>>>& cluster_list){
  fs::path project_root = fs::current_path().parent_path();
  fs::path output_dir = project_root /"data"/ "consensus_files";

    bool created = fs::create_directories(output_dir);
    if (!created) {
        return;
    }
    
    for (const auto& sample_pair : cluster_list) {
        const std::string& sample_id = sample_pair.first;
        const std::vector<std::vector<std::string>>& clusters = sample_pair.second;

        
        std::string filename = output_dir / (sample_id + ".fasta");
        std::ofstream outfile(filename);

        if (!outfile) {
            std::cerr << "Greška pri otvaranju fajla: " << filename << "\n";
            continue;
        }

        int cluster_idx = 0;
        for (const auto& cluster : clusters) {
            std::string centroid = computeCentroid(cluster);
            outfile << ">" << sample_id << "_cluster" << cluster_idx << "\n";
            outfile << centroid << "\n";
            ++cluster_idx;
        }

        outfile.close();
       
    }
}


int smithWaterman(const std::string& seq1, const std::string& seq2, int match = 2, int mismatch = -1, int gap = -1) {
    int m = seq1.size();
    int n = seq2.size();
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));

    int maxScore = 0;

    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int score = (seq1[i - 1] == seq2[j - 1]) ? match : mismatch;
            dp[i][j] = std::max({
                0,
                dp[i - 1][j - 1] + score,
                dp[i - 1][j] + gap,
                dp[i][j - 1] + gap
            });
            maxScore = std::max(maxScore, dp[i][j]);
        }
    }

    return maxScore;
}

//compare each cluster centroid with every ground truth seqence using SmithWaterman algorithm. Write out best match with score
void compareCentroidsWithGT(const std::vector<std::string>& centroids, const std::unordered_map<std::string, std::string>& gt_sequences) {
    int idx=0;
    for (const auto& centroid : centroids) {
        int bestScore = std::numeric_limits<int>::min();
        std::string bestGT_ID;

        for (const auto& [gt_id, gt_seq] : gt_sequences) {
            int score = smithWaterman(centroid, gt_seq);
            if (score > bestScore) {
                bestScore = score;
                bestGT_ID = gt_id;
            }
        }

        std::cout << "Cluster" << idx++ << " matched best with GT ID: " 
                  << bestGT_ID << " (score: " << bestScore << "/498)\n"; //498 score is maximum that sw algorithm can give in our case
    }
}


//parse ground truth and cluster files, then compare each cluster consensus to ground truth sequences
void checkGroundTruth(fs::path path, int index){

    auto groundTruthSequences = parseGT_File(path.string());
    fs::path project_root = fs::current_path().parent_path();
    fs::path directory = project_root /"data"/ "consensus_files";

    std::vector<std::string> sequencesToCheck;
    if(index == 29){
      cout << "Checking J29.fasta" << endl;
      sequencesToCheck = parse_fasta_sequences(directory / "J29.fasta");
    }
    else if(index == 30){
      cout << "checking J30.fasta" << endl;
      sequencesToCheck = parse_fasta_sequences(directory / "J30.fasta");
    }
    
    compareCentroidsWithGT(sequencesToCheck, groundTruthSequences);

}



int main(int argc, char **argv) {
  if (argc < 3) {
    cerr << "Missing folder path/s!" << endl;
    return 1;
  }
  string sampleFolder = argv[1]; 
  vector<fs::path> filesToParse = iterateDirectory(sampleFolder);


  std::unordered_map<std::string, std::vector<std::string> > readingsList = get_readingsList(filesToParse, 296);

  /*
  cout << "readingsList size: " << readingsList.size() << endl;
  for (const auto& pair : readingsList) {
    cout << "readings size: " << pair.second.size() << endl;
  }
  cout << "readingsList[0]" << endl;


  for (const auto& a: readingsList[0]) {
    cout << a <<"||"<< endl;
  }*/

  //cluster_parameter_test(readingsList["J30"]);

  //vector<vector<std::string>> GT_readingsList;
  /*
  vector<fs::path> GT_Files = iterateDirectory(GT_Folder);
  std::unordered_map<std::string, std::vector<std::string> > GT_readingsList = get_readingsList(GT_Files, 0);*/
  //cout << "GT_readingsList size: " << GT_readingsList.size() << endl;

  generateMSA_Files(readingsList);
    //map with sample id as key and cluster list ad value
    std::unordered_map<std::string, std::vector<std::vector<std::string>> > cluster_list = create_clusters(10, 0);
    cout<<"cluster_list size: "<<cluster_list.size()<<endl;
    cout<<"cluster lists sizes:"<<endl;
    for (const auto& list: cluster_list) {
        cout<<list.second.size()<<endl;
    }

  //generate .fasta files for every cluster in each sample file  
  writeCentroidsToFasta(cluster_list);

  //check if we get good results
  string GT_Folder = argv[2];
  vector<fs::path> GT_files = iterateDirectory(GT_Folder);

  int index = 29;
  for (const auto& path : GT_files) {
    checkGroundTruth(path, index++);
  }

  return 0;

}
