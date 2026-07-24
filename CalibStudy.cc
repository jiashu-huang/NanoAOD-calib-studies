#include <TBranch.h>
#include <TClass.h>
#include <TDataType.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TObjArray.h>
#include <TTree.h>

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>

namespace {

std::string trim(const std::string& text) {
  const auto first = text.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return "";
  const auto last = text.find_last_not_of(" \t\r\n");
  return text.substr(first, last - first + 1);
}

std::string inputRootFile(const std::string& inputList) {
  std::ifstream stream(inputList);
  if (!stream) throw std::runtime_error("Cannot open input address file: " + inputList);

  std::string line;
  while (std::getline(stream, line)) {
    line = trim(line);
    if (line.empty() || line.front() == '#') continue;
    if (line.size() < 5 || line.substr(line.size() - 5) != ".root") {
      throw std::runtime_error("First non-comment input line is not a .root file: " + line);
    }
    return line;
  }
  throw std::runtime_error("No .root file found in input address file: " + inputList);
}

bool validProductionTag(const std::string& tag) {
  return !tag.empty() && tag != "." && tag != ".." &&
         tag.find('/') == std::string::npos && tag.find('\\') == std::string::npos &&
         tag.find('\0') == std::string::npos;
}

std::string branchType(TBranch* branch) {
  const char* className = branch->GetClassName();
  if (className && *className) return className;

  TClass* expectedClass = nullptr;
  EDataType expectedType = kOther_t;
  if (branch->GetExpectedType(expectedClass, expectedType) == 0) {
    if (expectedClass) return expectedClass->GetName();
    if (expectedType != kOther_t) return TDataType::GetTypeName(expectedType);
  }

  auto* leaves = branch->GetListOfLeaves();
  if (leaves && leaves->GetEntriesFast() > 0) {
    if (leaves->GetEntriesFast() == 1) {
      return static_cast<TLeaf*>(leaves->At(0))->GetTypeName();
    }
    std::ostringstream types;
    for (int i = 0; i < leaves->GetEntriesFast(); ++i) {
      if (i) types << ", ";
      auto* leaf = static_cast<TLeaf*>(leaves->At(i));
      types << leaf->GetName() << ':' << leaf->GetTypeName();
    }
    return types.str();
  }
  return "unknown";
}

void writeFields(TTree& tree, const std::filesystem::path& report) {
  std::ofstream output(report);
  if (!output) throw std::runtime_error("Cannot write report: " + report.string());
  output << "# Branch schema for TTree '" << tree.GetName() << "'\n";
  output << "# branch_name\tROOT/C++ data_type\n";
  auto* branches = tree.GetListOfBranches();
  for (int i = 0; i < branches->GetEntriesFast(); ++i) {
    auto* branch = static_cast<TBranch*>(branches->At(i));
    output << branch->GetName() << '\t' << branchType(branch) << '\n';
  }
}

void writeFirstEvent(TTree& tree, const std::filesystem::path& report) {
  std::ofstream output(report);
  if (!output) throw std::runtime_error("Cannot write report: " + report.string());
  output << "# Entry 0 dump for TTree '" << tree.GetName() << "'\n";
  if (tree.GetEntries() == 0) {
    output << "# The tree contains no entries.\n";
    return;
  }
  output.close();

  const int originalStdout = dup(fileno(stdout));
  if (originalStdout == -1) throw std::runtime_error("Cannot save stdout: " + std::string(std::strerror(errno)));
  if (!std::freopen(report.c_str(), "a", stdout)) {
    close(originalStdout);
    throw std::runtime_error("Cannot redirect stdout to: " + report.string());
  }
  tree.Show(0);
  std::fflush(stdout);
  if (dup2(originalStdout, fileno(stdout)) == -1) {
    close(originalStdout);
    throw std::runtime_error("Cannot restore stdout: " + std::string(std::strerror(errno)));
  }
  close(originalStdout);
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc != 3 && argc != 4) {
    std::cerr << "Usage: " << argv[0] << " <input-address-txt> <production-tag> [tree-name]\n";
    return 2;
  }

  try {
    const std::string inputList = argv[1];
    const std::string productionTag = argv[2];
    const std::string treeName = argc == 4 ? argv[3] : "Events";
    if (!validProductionTag(productionTag)) {
      throw std::runtime_error("production-tag must be one safe path component (no '/' or '\\\\')");
    }

    const auto rootFile = inputRootFile(inputList);
    const std::filesystem::path outputDirectory =
        std::filesystem::path("outputs") / productionTag;
    std::filesystem::create_directories(outputDirectory);

    std::unique_ptr<TFile> file(TFile::Open(rootFile.c_str(), "READ"));
    if (!file || file->IsZombie()) throw std::runtime_error("Cannot open ROOT file: " + rootFile);
    auto* tree = dynamic_cast<TTree*>(file->Get(treeName.c_str()));
    if (!tree) throw std::runtime_error("TTree '" + treeName + "' not found in: " + rootFile);

    writeFields(*tree, outputDirectory / "fields.txt");
    writeFirstEvent(*tree, outputDirectory / "0th-event.txt");
    std::cout << "Wrote reports to " << outputDirectory << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "calib-study: " << error.what() << '\n';
    return 1;
  }
}
