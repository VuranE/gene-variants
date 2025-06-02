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


  //use spoa to generate consensus and msa of given sequences
pair<string, vector<string>> consensusAndMSA(const vector<FastqRead>& reads){
  auto alignmentEngine = spoa::AlignmentEngine::Create(
      spoa::AlignmentType::kNW, 3, -5, -3);  // linear gaps

  spoa::Graph graph{};

  for (const auto& read : reads) {
    auto alignment = alignment_engine->Align(read.sequence, graph);
    graph.AddAlignment(alignment, read.sequence);
  }

  auto consensus = graph.GenerateConsensus();

  auto msa = graph.GenerateMultipleSequenceAlignment();

  return {consensus, msa};

}

//clustering of sequences in each file
vector<string> runClustering(const vector<FastqRead>& reads){
  auto [consensus, msa] = consensusAndMSA(reads);

}

int main(int argc, char** argv) {

  if(argc != 2){
    cerr << "Missing samples folder path!" << endl;
    return 1;
  }

  string sampleFolder = argv[1];
  vector<fs::path> filesToParse = iterateDirectory(sampleFolder);

  for (const auto& filePath : filesToParse) {
    vector<FastqRead> reads = parseFile(filePath);
    
    vector<string> consensusVariants = runClustering(reads);

    for (size_t i = 0; i < consensusVariants.size() && i < 4; ++i) {
      cout << ">consensus_" << file.filename().string() << "_v" << i + 1 << endl;
      cout << consensusVariants[i] << endl;
    }
  }

  /*
    GLOBALNO PORAVNANJE NEEDLEMANWUNSCH
    vidi: https://github.com/francescoborando/Sequence_Alignment/blob/main/NeedlemanWunsch.cpp

    POGLAVLJE 4.5
    https://repozitorij.vef.unizg.hr/islandora/object/vef%3A641
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