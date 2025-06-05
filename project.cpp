#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <fstream>
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
vector<FastqRead> parseFile(const string& filename, size_t targetLength = 296){
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


    if(read.sequence.length() == targetLength) {
      reads.push_back(read);
    }
    else {
      n_longer++;
    }
  }
  //printf("n_longer: %d",n_longer);

  return reads;

}

/*Function for getting filenames of sample reads. It takes path to directory with all fastq samples as parameter and returns vector<path> that
  includes all desired file paths*/
vector<fs::path> iterateDirectory(const string& dirPath){
    vector<fs::path> files;
    for(const auto& entry : fs::directory_iterator(dirPath)){
      if(entry.path().extension() == ".fastq"){ //take into consideration only fastq files
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
      spoa::AlignmentType::kNW, 5, -4, -8);
  cout << "DEBUG 1" << endl;
  spoa::Graph graph{};
  cout << "DEBUG 2" << endl;

  for (const auto& read : reads) {
    auto alignment = alignmentEngine->Align(read, graph);
    graph.AddAlignment(alignment, read);
  }
  return graph;
}


// Hamming distance (number of mismatches)
size_t hammingDistance(const string &a, const string &b) {
  size_t dist = 0;
  for (size_t i = 0; i < a.size(); ++i)
    if (a[i] != b[i]) ++dist;
  return dist;
}

//iterates over reading and generates clusters
//k = number of missmatches tolerated for one cluster
//s = min nuber of cluster members

vector<vector<string> > clusteringAlgorithm(const vector<string> &readings, size_t k, size_t minClusterSize) {
  vector<vector<string> > clusters;

  for (const string &seq: readings) {
    bool added = false;

    // Try adding to an existing cluster
    for (auto &cluster: clusters) {
      // Compare to the cluster's representative (first sequence)
      if (hammingDistance(seq, cluster[0]) <= k) {
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
  cout << "Unfiltered clusters count: " << clusters.size() << endl;
  for (const auto& cluster : clusters){
    cout << "  Cluster size: " << cluster.size() << endl;
    if (cluster.size() > minClusterSize)
      filtered.push_back(cluster);
  }

  return filtered;
}


int main(int argc, char **argv) {
  if (argc < 3) {
    cerr << "Missing folder path/s!" << endl;
    return 1;
  }

  string sampleFolder = argv[1];
  vector<fs::path> filesToParse = iterateDirectory(sampleFolder);

  
  /*combine all acceptable sequences to one vector*/
  vector<std::string> chosenSequences;
  for (const auto& filePath : filesToParse) {
    cout << "parsing file " << filePath << endl;
    vector<FastqRead> reads = parseFile(filePath.string());
    for (const auto& read : reads) {
      chosenSequences.push_back(read.sequence);
    }
  }
  string GT_Folder = argv[2];
  vector<std::string> GT_Sequences;
  vector<fs::path> GT_Files = iterateDirectory(GT_Folder);
  for (const auto& GT_filePath : GT_Files) {
    vector<FastqRead> reads = parseFile(GT_filePath.string());
    for (const auto& read : reads) {
      GT_Sequences.push_back(read.sequence);
    }
  }

  vector<vector<string>> clusters = clusteringAlgorithm(GT_Sequences, 100, 0);
  cout << clusters.size() << endl;


  /*
  for (auto& cluster : clusters) {

  }
  auto graph = generateGraph(chosenSequences);

  //use spoa to get consensus and msa
  auto consensus = graph.GenerateConsensus();
  auto msa = graph.GenerateMultipleSequenceAlignment();

  cout << "msa:" << endl;
  for (const auto& seq : msa) {
    cout << seq << endl;
  }

  */

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