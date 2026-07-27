#include <TBranch.h>
#include <TCanvas.h>
#include <TClass.h>
#include <TDataType.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLeaf.h>
#include <TObjArray.h>
#include <TTree.h>

#include <algorithm>
#include <cerrno>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

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

int jetIndex(const std::string& branchName, const std::string& prefix) {
  if (branchName.rfind(prefix, 0) != 0 || branchName.size() == prefix.size()) return -1;

  int index = 0;
  for (std::size_t position = prefix.size(); position < branchName.size(); ++position) {
    const char character = branchName[position];
    if (character < '0' || character > '9') return -1;
    index = index * 10 + (character - '0');
  }
  return index;
}

std::map<int, TLeaf*> indexedLeaves(TTree& tree, const std::string& prefix) {
  std::map<int, TLeaf*> leaves;
  auto* branches = tree.GetListOfBranches();
  for (int branchNumber = 0; branchNumber < branches->GetEntriesFast(); ++branchNumber) {
    auto* branch = static_cast<TBranch*>(branches->At(branchNumber));
    const int index = jetIndex(branch->GetName(), prefix);
    if (index < 0) continue;

    auto* leaf = tree.GetLeaf(branch->GetName());
    if (!leaf) {
      throw std::runtime_error("Cannot read scalar leaf for branch: " +
                               std::string(branch->GetName()));
    }
    leaves.emplace(index, leaf);
  }
  if (leaves.empty()) throw std::runtime_error("No branches found with prefix: " + prefix);
  return leaves;
}

struct CutflowBranches {
  std::map<int, TLeaf*> pt;
  std::map<int, TLeaf*> eta;
  std::array<std::map<int, TLeaf*>, 4> matched;
};

CutflowBranches cutflowBranches(TTree& tree) {
  CutflowBranches branches{
      indexedLeaves(tree, "ak4JetPt"),
      indexedLeaves(tree, "ak4JetEta"),
      {indexedLeaves(tree, "ak4MatchedHadB_"), indexedLeaves(tree, "ak4MatchedLepB_"),
       indexedLeaves(tree, "ak4MatchedHadQ1_"), indexedLeaves(tree, "ak4MatchedHadQ2_")}};

  for (const auto& [index, _] : branches.pt) {
    if (!branches.eta.count(index)) {
      throw std::runtime_error("Missing ak4JetEta branch for jet index " + std::to_string(index));
    }
    for (const auto& matched : branches.matched) {
      if (!matched.count(index)) {
        throw std::runtime_error("Missing matched-jet branch for jet index " + std::to_string(index));
      }
    }
  }
  return branches;
}

std::vector<int> matchedJetIndices(const std::map<int, TLeaf*>& matched) {
  std::vector<int> indices;
  for (const auto& [index, flag] : matched) {
    if (flag->GetValue() == 1.0) indices.push_back(index);
  }
  return indices;
}

double eventWeight(TLeaf* finalWeight) {
  return finalWeight ? finalWeight->GetValue() : 1.0;
}

struct CutflowStep {
  std::string label;
  Long64_t events;
  double weightedEvents;
};

void writeCutflow(TTree& tree, const std::filesystem::path& report) {
  const CutflowBranches branches = cutflowBranches(tree);
  auto* finalWeight = tree.GetLeaf("finalWeight");
  const Long64_t initial = tree.GetEntries();
  Long64_t fourSelectedJets = 0;
  Long64_t allQuarksMatched = 0;
  Long64_t distinctMatchedJets = 0;
  double initialWeight = 0.0;
  double fourSelectedJetsWeight = 0.0;
  double allQuarksMatchedWeight = 0.0;
  double distinctMatchedJetsWeight = 0.0;

  for (Long64_t entry = 0; entry < initial; ++entry) {
    tree.GetEntry(entry);
    const double weight = eventWeight(finalWeight);
    initialWeight += weight;

    int selectedJets = 0;
    for (const auto& [index, pt] : branches.pt) {
      const double eta = branches.eta.at(index)->GetValue();
      if (pt->GetValue() >= 20.0 && std::abs(eta) < 2.5) ++selectedJets;
    }
    if (selectedJets < 4) continue;
    ++fourSelectedJets;
    fourSelectedJetsWeight += weight;

    std::array<std::vector<int>, 4> quarkMatchedJetIndices;
    bool everyQuarkMatched = true;
    for (std::size_t quark = 0; quark < branches.matched.size(); ++quark) {
      quarkMatchedJetIndices[quark] = matchedJetIndices(branches.matched[quark]);
      if (quarkMatchedJetIndices[quark].empty()) everyQuarkMatched = false;
    }
    if (!everyQuarkMatched) continue;
    ++allQuarksMatched;
    allQuarksMatchedWeight += weight;

    std::set<int> uniqueMatchedJets;
    std::size_t numberOfMatches = 0;
    for (const auto& quarkMatches : quarkMatchedJetIndices) {
      numberOfMatches += quarkMatches.size();
      uniqueMatchedJets.insert(quarkMatches.begin(), quarkMatches.end());
    }
    // A valid one-to-one matching has one jet for every quark and four unique jets.
    if (numberOfMatches == 4 && uniqueMatchedJets.size() == 4) {
      ++distinctMatchedJets;
      distinctMatchedJetsWeight += weight;
    }
  }

  std::ofstream output(report);
  if (!output) throw std::runtime_error("Cannot write report: " + report.string());
  output << "# Cutflow for TTree '" << tree.GetName() << "'\n";
  output << "# Jet requirement: PT >= 20 GeV and abs(eta) < 2.5\n";
  output << "# Matched branches: ak4MatchedHadB_*, ak4MatchedLepB_*, "
            "ak4MatchedHadQ1_*, ak4MatchedHadQ2_*\n";
  output << "# The final step requires exactly one matched jet per quark and four unique jet indices.\n\n";
  output << "# Efficiencies use " << (finalWeight ? "finalWeight" : "unit event weights") << ".\n\n";
  output << std::left << std::setw(54) << "cut" << std::right << std::setw(12) << "events"
         << std::setw(18) << "weighted yield" << std::setw(18) << "vs previous"
         << std::setw(18) << "vs initial" << '\n';
  output << std::string(120, '-') << '\n';

  const std::array<CutflowStep, 4> steps = {{
      {"Initial input", initial, initialWeight},
      {"At least four jets passing PT and eta requirements", fourSelectedJets, fourSelectedJetsWeight},
      {"All four semileptonic-top quarks matched", allQuarksMatched, allQuarksMatchedWeight},
      {"Four matched jets are distinct", distinctMatchedJets, distinctMatchedJetsWeight},
  }};
  output << std::fixed << std::setprecision(2);
  double previousWeight = initialWeight;
  for (const auto& step : steps) {
    const double previousEfficiency =
        previousWeight == 0.0 ? 0.0 : 100.0 * step.weightedEvents / previousWeight;
    const double initialEfficiency =
        initialWeight == 0.0 ? 0.0 : 100.0 * step.weightedEvents / initialWeight;
    output << std::left << std::setw(54) << step.label << std::right << std::setw(12) << step.events
           << std::setw(18) << step.weightedEvents << std::setw(17) << previousEfficiency << '%'
           << std::setw(17) << initialEfficiency
           << "%\n";
    previousWeight = step.weightedEvents;
  }
}

void requireKinematicLeaves(const CutflowBranches& branches,
                            const std::map<int, TLeaf*>& phi,
                            const std::map<int, TLeaf*>& mass) {
  for (const auto& [index, _] : branches.pt) {
    if (!phi.count(index) || !mass.count(index)) {
      throw std::runtime_error("Missing ak4JetPhi or ak4JetMass branch for jet index " +
                               std::to_string(index));
    }
  }
}

double invariantMass(const CutflowBranches& branches,
                     const std::map<int, TLeaf*>& phi,
                     const std::map<int, TLeaf*>& mass,
                     const std::vector<int>& jetIndices) {
  double energy = 0.0;
  double momentumX = 0.0;
  double momentumY = 0.0;
  double momentumZ = 0.0;
  for (const int index : jetIndices) {
    const double pt = branches.pt.at(index)->GetValue();
    const double eta = branches.eta.at(index)->GetValue();
    const double jetPhi = phi.at(index)->GetValue();
    const double jetMass = mass.at(index)->GetValue();
    const double px = pt * std::cos(jetPhi);
    const double py = pt * std::sin(jetPhi);
    const double pz = pt * std::sinh(eta);
    energy += std::sqrt(px * px + py * py + pz * pz + jetMass * jetMass);
    momentumX += px;
    momentumY += py;
    momentumZ += pz;
  }
  const double massSquared = energy * energy - momentumX * momentumX - momentumY * momentumY -
                             momentumZ * momentumZ;
  return std::sqrt(std::max(0.0, massSquared));
}

void writeMassPlot(TTree& tree,
                   const std::string& productionTag,
                   const std::filesystem::path& report,
                   const std::string& histogramName,
                   const std::string& branchName,
                   const std::vector<std::size_t>& quarks,
                   int bins,
                   double minimum,
                   double maximum) {
  const CutflowBranches branches = cutflowBranches(tree);
  const auto phi = indexedLeaves(tree, "ak4JetPhi");
  const auto mass = indexedLeaves(tree, "ak4JetMass");
  requireKinematicLeaves(branches, phi, mass);
  auto* finalWeight = tree.GetLeaf("finalWeight");

  TH1D histogram(histogramName.c_str(),
                 (branchName + " in " + productionTag + ", post CMSSW_CHARGE").c_str(),
                 bins,
                 minimum,
                 maximum);
  histogram.SetDirectory(nullptr);
  histogram.GetXaxis()->SetTitle("Invariant mass [GeV]");
  histogram.GetYaxis()->SetTitle(finalWeight ? "Weighted events" : "Events");
  histogram.Sumw2();

  for (Long64_t entry = 0; entry < tree.GetEntries(); ++entry) {
    tree.GetEntry(entry);
    std::vector<int> jetIndices;
    bool hasOneDistinctMatchPerQuark = true;
    for (const std::size_t quark : quarks) {
      const std::vector<int> matches = matchedJetIndices(branches.matched.at(quark));
      if (matches.size() != 1) {
        hasOneDistinctMatchPerQuark = false;
        break;
      }
      jetIndices.push_back(matches.front());
    }
    std::set<int> uniqueJetIndices(jetIndices.begin(), jetIndices.end());
    if (!hasOneDistinctMatchPerQuark || uniqueJetIndices.size() != jetIndices.size()) continue;
    histogram.Fill(invariantMass(branches, phi, mass, jetIndices), eventWeight(finalWeight));
  }

  TCanvas canvas((histogramName + "_canvas").c_str(), "", 800, 600);
  histogram.Draw("HIST");
  canvas.SaveAs(report.c_str());
  if (!std::filesystem::exists(report)) {
    throw std::runtime_error("Cannot write plot: " + report.string());
  }
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
    writeCutflow(*tree, outputDirectory / "cutflow.txt");
    writeMassPlot(*tree,
                  productionTag,
                  outputDirectory / "hadronic-w-mass.pdf",
                  "hadronic_w_mass",
                  "ak4MatchedHadQ1_ak4MatchedHadQ2_mass",
                  {2, 3},
                  10,
                  30.0,
                  130.0);
    writeMassPlot(*tree,
                  productionTag,
                  outputDirectory / "hadronic-top-mass.pdf",
                  "hadronic_top_mass",
                  "ak4MatchedHadQ1_ak4MatchedHadQ2_ak4MatchedHadB_mass",
                  {2, 3, 0},
                  14,
                  100.0,
                  240.0);
    std::cout << "Wrote reports to " << outputDirectory << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "calib-study: " << error.what() << '\n';
    return 1;
  }
}
