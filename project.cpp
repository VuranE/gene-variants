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

int main(int argc, char** argv) {

  if(argc != 2){
    cerr << "Missing samples folder path!" << endl;
    return 1;
  }

  string sampleFolder = argv[1];
  vector<fs::path> filesToParse = iterateDirectory(sampleFolder);

  vector<FastqRead> chosenSequences; //vector containing all sequences in all files that are of desired length
  for (const auto& filePath : filesToParse) {
    vector<FastqRead> reads = parseFile(filePath);
    chosenSequences.insert(chosenSequences.end(), reads.begin(), reads.end());
  }

  /*
    GLOBALNO PORAVNANJE NEEDLEMANWUNSCH
    vidi: https://github.com/francescoborando/Sequence_Alignment/blob/main/NeedlemanWunsch.cpp
  */


  


  return 0;


    /*
    spoa example:

  std::vector<std::string> sequences = {
      "CATAAAAGAACGTAGGTCGCCCGTCCGTAACCTGTCGGATCACCGGAAAGGACCCGTAAAGTGATAATGAT",
      "ATAAAGGCAGTCGCTCTGTAAGCTGTCGATTCACCGGAAAGATGGCGTTACCACGTAAAGTGATAATGATTAT",
      "ATCAAAGAACGTGTAGCCTGTCCGTAATCTAGCGCATTTCACACGAGACCCGCGTAATGGG",
      "CGTAAATAGGTAATGATTATCATTACATATCACAACTAGGGCCGTATTAATCATGATATCATCA",
      "GTCGCTAGAGGCATCGTGAGTCGCTTCCGTACCGCAAGGATGACGAGTCACTTAAAGTGATAAT",
      "CCGTAACCTTCATCGGATCACCGGAAAGGACCCGTAAATAGACCTGATTATCATCTACAT"
  };

  auto alignment_engine = spoa::AlignmentEngine::Create(
      spoa::AlignmentType::kNW, 3, -5, -3);  // linear gaps

  spoa::Graph graph{};

  for (const auto& it : sequences) {
    auto alignment = alignment_engine->Align(it, graph);
    graph.AddAlignment(alignment, it);
  }

  auto consensus = graph.GenerateConsensus();

  std::cerr << ">Consensus LN:i:" << consensus.size() << std::endl
            << consensus << std::endl;

  auto msa = graph.GenerateMultipleSequenceAlignment();

  for (const auto& it : msa) {
    std::cerr << it << std::endl;
  }

  return 0;
  */
}