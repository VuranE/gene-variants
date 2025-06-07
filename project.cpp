#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <fstream>
#include <unordered_map>

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
  //cout << "parser" << endl;
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
  //printf("n_diff: %d\n",n_longer);
  if (reads.empty()) {
    cerr << "Error parsing file: " << filename <<"for target length:"<<targetLength<<endl;
  }
  return reads;

}


vector<std::string> parseGT_File(const string& filename, size_t targetLength = 249){
  vector<std::string> reads;
  ifstream file(filename);

  if(!file) {
    cerr << "Error opening GT file: " << filename << endl;
    return reads;
  }
  cout<< "parsing GT file: "<<filename << endl;



  int n = 0;
  string line;
  while (getline(file, line)){
    if (line[0] == '>')
      continue;
    if(line.length() == targetLength) {
      reads.push_back(line);
      n++;
    }
  }
  printf("n: %d\n",n);

  return reads;

}

/*Function for getting filenames of sample reads. It takes path to directory with all fastq samples as parameter and returns vector<path> that
  includes all desired file paths*/
vector<fs::path> iterateDirectory(const string& dirPath){
    vector<fs::path> files;
    for(const auto& entry : fs::directory_iterator(dirPath)){
      if(entry.path().extension() == ".fastq" or entry.path().extension() == ".fasta"){ //take into consideration only fastq or fasta files
        //if(entry.path().filename().string().compare(0, 1, "J") == 0) //take into consideration only files that starts with "J", indicating deer samples (can be changed for different purposes)
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
//
void generateMSA_File(const unordered_map<std::string, std::vector<std::string>> & readingsList) {

}

int main(int argc, char **argv) {
  if (argc < 3) {
    cerr << "Missing folder path/s!" << endl;
    return 1;
  }

  string sampleFolder = argv[1];
  vector<fs::path> filesToParse = iterateDirectory(sampleFolder);

  
  /*combine all acceptable sequences to one vector*/
  //vector<vector<std::string>> readingsList;
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

  string GT_Folder = argv[2];
  vector<std::string> GT_Sequences;
  //vector<vector<std::string>> GT_readingsList;
  vector<fs::path> GT_Files = iterateDirectory(GT_Folder);
  std::unordered_map<std::string, std::vector<std::string> > GT_readingsList = get_readingsList(GT_Files, 0);


  cout << "GT_readingsList size: " << GT_readingsList.size() << endl;
  for (const auto& pair : GT_readingsList) {
    cout << "readings size: " << pair.second.size() << endl;
  }
  //cluster_parameter_test(GT_readingsList[0]);
  //vector<vector<string>> clusters = clusteringAlgorithm(GT_readingsList["J30B"], 12, 0);
  //cout << clusters.size() << endl;

  //cout << "readingsList[J30] size:"<< readingsList["J30"].size() << endl;
  auto graph = generateGraph(readingsList["J30"]);

  generateMSA_File(readingsList);

  //use spoa to get consensus and msa
  auto consensus = graph.GenerateConsensus();
  auto msa = graph.GenerateMultipleSequenceAlignment();

  cout << "msa size:"<< msa.size() << endl;/*
  for (const auto& seq : msa) {
    cout << seq << endl;
  }*/

  //cluster_parameter_test(msa);



  /*vector<vector<string>> clusters = clusteringAlgorithm(msa, 12, 100);
  
  int clusterCount = 1;
  for(vector<string> cluster: clusters){
    cout << "Cluster no." << clusterCount++ << endl;
    cout << "sequences:" << endl;

    int sequenceCount = 1;
    for(string sequence : cluster){
      cout << sequenceCount++ << ") " << sequence << endl;
    }
  }  */

  /*
    POGLAVLJE 4.5
    https://repozitorij.vef.unizg.hr/islandora/object/vef%3A641
  */

  return 0;

}
