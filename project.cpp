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
  cout << "parser" << endl;
  vector<FastqRead> reads;
  ifstream file(filename);

  if(!file) {
    cerr << "Error opening file: " << filename << endl;
    return reads;
  }

  
  string line;
  while (getline(file, line)){
    FastqRead read;
    read.id = line;
    getline(file, read.sequence);
    getline(file, read.plus);
    getline(file, read.quality);

    if(read.sequence.length() == targetLength)
      reads.push_back(read);
  }

  return reads;

}

/*Function for getting filenames of sample reads. It takes path to directory with all fastq samples as parameter and returns vector<path> that
  includes all desired file paths*/
vector<fs::path> iterateDirectory(const string& dirPath){
    vector<fs::path> files;
    for(const auto& entry : fs::directory_iterator(dirPath)){
      if(entry.path().extension() == ".fastq"){ //take into consideration only fastq files
        if(entry.path().filename().string().compare(0, 1, "J") == 0) //take into consideration only files that starts with "J", indicating deer samples (can be changed for different purpouses)
          files.push_back(entry.path());
          
      }
    }
    return files;
  }


  //use spoa to generate consensus and msa of given sequences
pair<string, vector<string>> consensusAndMSA(const vector<FastqRead>& reads){
  cout << "calculating consensus and msa" << endl;
  auto alignmentEngine = spoa::AlignmentEngine::Create(
      spoa::AlignmentType::kNW, 0, -1, -100); 

  spoa::Graph graph{};
 
  
  for (const auto& read : reads) {
    auto alignment = alignmentEngine->Align(read.sequence, graph);
    graph.AddAlignment(alignment, read.sequence);
  }
  
  auto consensus = graph.GenerateConsensus();
  
  auto msa = graph.GenerateMultipleSequenceAlignment();
  
  return {consensus, msa};

}

//
vector<vector<string>> clusteringAlgorithm(const vector<string>& msa, size_t k, size_t s) {
    vector<vector<string>> clusters;

    for (const string& seq : msa) {
        bool added = false;

        for (auto& cluster : clusters) {
            bool isSimilar = true;
            for (const string& member : cluster) {
                size_t diff = 0;
                for (size_t i = 0; i < seq.size(); ++i) {
                    if (seq[i] != member[i]) ++diff;
                    if (diff >= k) {
                        isSimilar = false;
                        break;
                    }
                }
                if (!isSimilar) break;
            }

            if (isSimilar) {
                cluster.push_back(seq);
                added = true;
                break;
            }
        }

        if (!added) {
            clusters.push_back({seq}); // generate new cluster
        }
    }

    // filter only clusters bigger than "s" parameter
    vector<vector<string>> filteredClusters;
    for (const auto& cluster : clusters) {
        if (cluster.size() > s) {
            filteredClusters.push_back(cluster);
        }
    }

    return filteredClusters;
}


int main(int argc, char** argv) {

  if(argc != 2){
    cerr << "Missing samples folder path!" << endl;
    return 1;
  }

  string sampleFolder = argv[1];
  vector<fs::path> filesToParse = iterateDirectory(sampleFolder);

  
  /*combine all acceptable sequences to one vector*/
  vector<FastqRead> chosenSequences;
  for (const auto& filePath : filesToParse) {
    cout << "parsing file " << filePath << endl;
    vector<FastqRead> reads = parseFile(filePath);
    chosenSequences.insert(chosenSequences.end(), reads.begin(), reads.end());
  }

  //use spoa to get consensus and msa
  auto [consensus, msa] = consensusAndMSA(chosenSequences);

  cout << consensus << endl;

  vector<vector<string>> clusters = clusteringAlgorithm(msa, 12, 100);
  
  int clusterCount = 1;
  for(vector<string> cluster: clusters){
    cout << "Cluster no." << clusterCount++ << endl;
    cout << "sequences:" << endl;

    int sequenceCount = 1;
    for(string sequence : cluster){
      cout << sequenceCount++ << ") " << sequence << endl;
    }
  }

  /*
    POGLAVLJE 4.5
    https://repozitorij.vef.unizg.hr/islandora/object/vef%3A641
  */

  return 0;

}