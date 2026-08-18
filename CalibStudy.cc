#include <TBranch.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TClass.h>
#include <TDataType.h>
#include <TFile.h>
#include <TH1D.h>
#include <TKey.h>
#include <TLeaf.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
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
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <limits>
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

std::vector<std::string> inputRootFiles(const std::string& inputList) {
  std::ifstream stream(inputList);
  if (!stream) throw std::runtime_error("Cannot open input address file: " + inputList);

  std::vector<std::string> rootFiles;
  std::string line;
  std::size_t lineNumber = 0;
  while (std::getline(stream, line)) {
    ++lineNumber;
    line = trim(line);
    if (line.empty() || line.front() == '#') continue;
    if (line.size() < 5 || line.substr(line.size() - 5) != ".root") {
      throw std::runtime_error("Input line " + std::to_string(lineNumber) +
                               " is not a .root file: " + line);
    }
    rootFiles.push_back(line);
  }
  if (rootFiles.empty()) {
    throw std::runtime_error("No .root files found in input address file: " + inputList);
  }
  return rootFiles;
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
  if (tree.LoadTree(0) < 0) {
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

std::map<int, std::string> indexedLeaves(TTree& tree, const std::string& prefix) {
  std::map<int, std::string> leaves;
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
    leaves.emplace(index, branch->GetName());
  }
  if (leaves.empty()) throw std::runtime_error("No branches found with prefix: " + prefix);
  return leaves;
}

double leafValue(TTree& tree, const std::string& branchName) {
  auto* leaf = tree.GetLeaf(branchName.c_str());
  if (!leaf) throw std::runtime_error("Cannot read scalar leaf for branch: " + branchName);
  return leaf->GetValue();
}

// The four partons of the semileptonic top decay, in the order their assigned
// jet slots are read, stored and reported.
const std::array<const char*, 4>& genJetIdxBranches() {
  static const std::array<const char*, 4> names = {"GenHadBJetIdx", "GenLepBJetIdx",
                                                   "GenHadQ1JetIdx", "GenHadQ2JetIdx"};
  return names;
}

struct CutflowBranches {
  std::map<int, std::string> pt;
  std::map<int, std::string> eta;
  std::array<std::string, 4> genJetIdx;
};

CutflowBranches cutflowBranches(TTree& tree) {
  CutflowBranches branches{indexedLeaves(tree, "ak4JetPt"), indexedLeaves(tree, "ak4JetEta"), {}};

  for (std::size_t parton = 0; parton < branches.genJetIdx.size(); ++parton) {
    const char* branchName = genJetIdxBranches()[parton];
    if (!tree.GetLeaf(branchName)) {
      throw std::runtime_error(
          "Missing gen-assignment branch: " + std::string(branchName) +
          ". Skims written before the one-to-one jet assignment carry ak4Matched*_N "
          "flags instead and are not supported.");
    }
    branches.genJetIdx[parton] = branchName;
  }

  for (const auto& [index, _] : branches.pt) {
    if (!branches.eta.count(index)) {
      throw std::runtime_error("Missing ak4JetEta branch for jet index " + std::to_string(index));
    }
  }
  return branches;
}

// The saved jet slot each parton is assigned to, -1 where it is unmatched.
std::array<int, 4> assignedJetIndices(TTree& tree, const CutflowBranches& branches) {
  std::array<int, 4> indices{};
  for (std::size_t parton = 0; parton < indices.size(); ++parton) {
    indices[parton] = static_cast<int>(std::llround(leafValue(tree, branches.genJetIdx[parton])));
  }
  return indices;
}

double eventWeight(TTree& tree) {
  auto* finalWeight = tree.GetLeaf("finalWeight");
  return finalWeight ? finalWeight->GetValue() : 1.0;
}

double invariantMass(TTree& tree,
                     const CutflowBranches& branches,
                     const std::map<int, std::string>& phi,
                     const std::map<int, std::string>& mass,
                     const std::vector<int>& jetIndices) {
  double energy = 0.0;
  double momentumX = 0.0;
  double momentumY = 0.0;
  double momentumZ = 0.0;
  for (const int index : jetIndices) {
    const double pt = leafValue(tree, branches.pt.at(index));
    const double eta = leafValue(tree, branches.eta.at(index));
    const double jetPhi = leafValue(tree, phi.at(index));
    const double jetMass = leafValue(tree, mass.at(index));
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

struct PtRegion {
  double minimum;
  double maximum;
  std::string label;
  std::string filename;
};

const std::array<PtRegion, 6>& ptRegions() {
  static const std::array<PtRegion, 6> regions = {{
      {20.0, 30.0, "[20, 30) GeV", "20to30"},
      {30.0, 50.0, "[30, 50) GeV", "30to50"},
      {50.0, 80.0, "[50, 80) GeV", "50to80"},
      {80.0, 110.0, "[80, 110) GeV", "80to110"},
      {110.0, 150.0, "[110, 150) GeV", "110to150"},
      {150.0, std::numeric_limits<double>::infinity(), "[150, inf) GeV", "150plus"},
  }};
  return regions;
}

int ptRegionIndex(double pt) {
  const auto& regions = ptRegions();
  for (std::size_t region = 0; region < regions.size(); ++region) {
    if (pt >= regions[region].minimum && pt < regions[region].maximum) {
      return static_cast<int>(region);
    }
  }
  return -1;
}

// The jets that the studies below score. One index is resolved per event and
// then shared by every study that plots the same jet.
enum class JetRole { LepB = 0, HadB = 1, C = 2, WB = 3 };
constexpr std::size_t kJetRoles = 4;

std::size_t roleIndex(JetRole role) { return static_cast<std::size_t>(role); }

struct PosvsNegGroup {
  std::string stem;
  std::string jetLabel;
  std::string directoryName;
  JetRole role;
};

struct PNetBGroup {
  std::string jetLabel;
  std::string directoryName;
  JetRole role;
};

std::vector<PosvsNegGroup> posvsNegGroups(bool runVcbStudies) {
  std::vector<PosvsNegGroup> groups = {{"lepb", "LepB", "LepB_jet_PosvsNeg", JetRole::LepB}};
  if (runVcbStudies) groups.push_back({"c", "C", "C_jet_PosvsNeg", JetRole::C});
  return groups;
}

std::vector<PNetBGroup> pNetBGroups(bool runVcbStudies) {
  std::vector<PNetBGroup> groups = {{"LepB", "LepB_jet_PNetB", JetRole::LepB},
                                    {"HadB", "HadB_jet_PNetB", JetRole::HadB}};
  if (runVcbStudies) {
    groups.push_back({"C", "C_jet_PNetB", JetRole::C});
    groups.push_back({"WB", "WB_jet_PNetB", JetRole::WB});
  }
  return groups;
}

std::string posvsNegHistogramName(const PosvsNegGroup& group,
                                  const PtRegion& region,
                                  bool positiveLepton) {
  return group.stem + "_part_posvsneg_" + region.filename +
         (positiveLepton ? "_positive" : "_negative");
}

std::string pNetBHistogramName(const PNetBGroup& group, const PtRegion& region) {
  return "pnetb_" + group.jetLabel + '_' + region.filename;
}

const std::array<double, 7>& pNetBBinEdges() {
  static const std::array<double, 7> edges = {0.0, 0.047, 0.245, 0.6734, 0.7862, 0.961, 1.0};
  return edges;
}

const std::array<const char*, 3>& cutflowLabels() {
  static const std::array<const char*, 3> labels = {
      "Initial input",
      "At least four jets passing PT and eta requirements",
      "All four semileptonic-top quarks matched"};
  return labels;
}

// Every study result lives here as a named histogram, so a worker can write its
// share to a partial file and the merge step can add them without knowing which
// study produced what.
class StudyHistograms {
 public:
  TH1D& adopt(std::unique_ptr<TH1D> histogram) {
    histogram->SetDirectory(nullptr);
    const std::string name = histogram->GetName();
    auto [entry, inserted] = histograms_.emplace(name, std::move(histogram));
    if (!inserted) throw std::runtime_error("Duplicate histogram name: " + name);
    return *entry->second;
  }

  TH1D& book(std::unique_ptr<TH1D> histogram) {
    TH1D& booked = adopt(std::move(histogram));
    booked.Sumw2();
    return booked;
  }

  TH1D& at(const std::string& name) {
    auto entry = histograms_.find(name);
    if (entry == histograms_.end()) throw std::runtime_error("Missing histogram: " + name);
    return *entry->second;
  }

  const TH1D& at(const std::string& name) const {
    auto entry = histograms_.find(name);
    if (entry == histograms_.end()) throw std::runtime_error("Missing histogram: " + name);
    return *entry->second;
  }

  bool has(const std::string& name) const { return histograms_.count(name) != 0; }
  std::size_t size() const { return histograms_.size(); }

  const std::map<std::string, std::unique_ptr<TH1D>>& all() const { return histograms_; }

  void add(const StudyHistograms& other, const std::string& label) {
    if (other.size() != size()) {
      throw std::runtime_error("Partial " + label + " holds " + std::to_string(other.size()) +
                               " histograms but " + std::to_string(size()) + " were expected. " +
                               "Were the partials produced with the same --Vcb setting?");
    }
    for (auto& [name, histogram] : histograms_) {
      if (!other.has(name)) {
        throw std::runtime_error("Partial " + label + " is missing histogram: " + name);
      }
      const TH1D& addend = other.at(name);
      if (addend.GetNbinsX() != histogram->GetNbinsX()) {
        throw std::runtime_error("Partial " + label + " has incompatible binning for: " + name);
      }
      histogram->Add(&addend);
    }
  }

 private:
  std::map<std::string, std::unique_ptr<TH1D>> histograms_;
};

StudyHistograms bookStudyHistograms(const std::string& productionTag,
                                    bool runVcbStudies,
                                    bool hasFinalWeight) {
  StudyHistograms histograms;
  const char* yieldTitle = hasFinalWeight ? "Weighted events" : "Events";
  const char* normalizedYieldTitle =
      hasFinalWeight ? "Normalized weighted events" : "Normalized events";

  for (const char* name : {"cutflow_events", "cutflow_weight"}) {
    TH1D& cutflow = histograms.book(std::make_unique<TH1D>(name, "Cutflow", 3, 0.0, 3.0));
    for (int bin = 1; bin <= 3; ++bin) {
      cutflow.GetXaxis()->SetBinLabel(bin, cutflowLabels()[bin - 1]);
    }
  }

  // Carries the finalWeight availability through a partial file so that the
  // merge step can label the cutflow report without reopening any input.
  TH1D& finalWeightFlag = histograms.book(
      std::make_unique<TH1D>("meta_has_final_weight", "finalWeight availability", 1, 0.0, 1.0));
  if (hasFinalWeight) finalWeightFlag.SetBinContent(1, 1.0);

  // The assignment is one-to-one by construction upstream, so this must stay at
  // zero; it is counted rather than asserted so one bad event cannot kill a job.
  histograms.book(std::make_unique<TH1D>("meta_nonunique_assignment",
                                         "Events whose four assigned jets are not distinct", 1,
                                         0.0, 1.0));

  const std::string titleSuffix = " in " + productionTag + ", post CMSSW_CHARGE";
  struct MassPlot {
    const char* name;
    const char* title;
    int bins;
    double minimum;
    double maximum;
  };
  for (const MassPlot& plot :
       {MassPlot{"hadronic_w_mass", "m(HadQ1 jet, HadQ2 jet)", 16, 0.0, 160.0},
        MassPlot{"hadronic_top_mass", "m(HadQ1 jet, HadQ2 jet, HadB jet)", 18,
                 80.0, 260.0}}) {
    TH1D& histogram = histograms.book(std::make_unique<TH1D>(
        plot.name, (plot.title + titleSuffix).c_str(), plot.bins, plot.minimum, plot.maximum));
    histogram.GetXaxis()->SetTitle("Invariant mass [GeV]");
    histogram.GetYaxis()->SetTitle(yieldTitle);
  }

  for (const PosvsNegGroup& group : posvsNegGroups(runVcbStudies)) {
    for (const PtRegion& region : ptRegions()) {
      TH1D& positive = histograms.book(std::make_unique<TH1D>(
          posvsNegHistogramName(group, region, true).c_str(),
          ("ParTPosvsNeg score of " + group.jetLabel + " jet, " + region.label).c_str(), 100, 0.0,
          1.0));
      positive.GetXaxis()->SetTitle("ak4JetParTPosvsNeg score");
      positive.GetYaxis()->SetTitle(yieldTitle);
      histograms.book(std::make_unique<TH1D>(
          posvsNegHistogramName(group, region, false).c_str(), "", 100, 0.0, 1.0));
    }
  }

  for (const PNetBGroup& group : pNetBGroups(runVcbStudies)) {
    for (const PtRegion& region : ptRegions()) {
      TH1D& histogram = histograms.book(std::make_unique<TH1D>(
          pNetBHistogramName(group, region).c_str(),
          ("PNetB score of " + group.jetLabel + " jet, " + region.label).c_str(),
          static_cast<int>(pNetBBinEdges().size() - 1), pNetBBinEdges().data()));
      histogram.GetXaxis()->SetTitle("ak4JetbtagPNetB score");
      histogram.GetXaxis()->SetLabelSize(0.0);
      histogram.GetYaxis()->SetTitle(normalizedYieldTitle);
    }
  }
  return histograms;
}

struct AnalysisBranches {
  CutflowBranches cutflow;
  std::map<int, std::string> phi;
  std::map<int, std::string> mass;
  std::map<int, std::string> parTPosvsNeg;
  std::map<int, std::string> pNetB;
  bool hasFinalWeight;
};

AnalysisBranches analysisBranches(TTree& tree, bool runVcbStudies) {
  AnalysisBranches branches{cutflowBranches(tree),
                            indexedLeaves(tree, "ak4JetPhi"),
                            indexedLeaves(tree, "ak4JetMass"),
                            indexedLeaves(tree, "ak4JetParTPosvsNeg"),
                            indexedLeaves(tree, "ak4JetbtagPNetB"),
                            tree.GetLeaf("finalWeight") != nullptr};

  for (const auto& [index, _] : branches.cutflow.pt) {
    if (!branches.phi.count(index) || !branches.mass.count(index)) {
      throw std::runtime_error("Missing ak4JetPhi or ak4JetMass branch for jet index " +
                               std::to_string(index));
    }
    if (!branches.parTPosvsNeg.count(index)) {
      throw std::runtime_error("Missing ak4JetParTPosvsNeg branch for jet index " +
                               std::to_string(index));
    }
    if (!branches.pNetB.count(index)) {
      throw std::runtime_error("Missing ak4JetbtagPNetB branch for jet index " +
                               std::to_string(index));
    }
  }
  if (!tree.GetLeaf("TriggerLeptonCharge")) {
    throw std::runtime_error("Missing required branch: TriggerLeptonCharge");
  }
  if (runVcbStudies) {
    for (const char* branchName : {"GenQ1PdgId", "GenQ2PdgId"}) {
      if (!tree.GetLeaf(branchName)) {
        throw std::runtime_error("Missing required branch: " + std::string(branchName));
      }
    }
  }
  return branches;
}

// Reading all 350-odd branches per entry dominates the run time, so keep only
// the ones the studies actually touch enabled.
void enableRequiredBranches(TTree& tree, const AnalysisBranches& branches, bool runVcbStudies) {
  std::set<std::string> required;
  const auto require = [&required](const std::map<int, std::string>& leaves) {
    for (const auto& [index, branchName] : leaves) required.insert(branchName);
  };
  require(branches.cutflow.pt);
  require(branches.cutflow.eta);
  for (const auto& branchName : branches.cutflow.genJetIdx) required.insert(branchName);
  require(branches.phi);
  require(branches.mass);
  require(branches.parTPosvsNeg);
  require(branches.pNetB);
  required.insert("TriggerLeptonCharge");
  if (branches.hasFinalWeight) required.insert("finalWeight");
  if (runVcbStudies) {
    required.insert("GenQ1PdgId");
    required.insert("GenQ2PdgId");
  }

  tree.SetBranchStatus("*", false);
  for (const std::string& branchName : required) {
    tree.SetBranchStatus(branchName.c_str(), true);
  }
}

// Returns the hadronic-W jet whose gen quark carries |pdgId|, or -1 when no
// quark or both quarks carry it.
int hadronicWJetByAbsPdgId(TTree& tree, int hadQ1Jet, int hadQ2Jet, int absPdgId) {
  const std::array<std::pair<const char*, int>, 2> quarks = {
      {{"GenQ1PdgId", hadQ1Jet}, {"GenQ2PdgId", hadQ2Jet}}};
  int selectedJet = -1;
  for (const auto& [pdgIdBranch, jet] : quarks) {
    const int pdgId = static_cast<int>(std::lround(leafValue(tree, pdgIdBranch)));
    if (std::abs(pdgId) != absPdgId) continue;
    if (selectedJet >= 0) return -1;
    selectedJet = jet;
  }
  return selectedJet;
}

// One pass over the tree fills every study, so adding a file costs one read of
// the enabled branches rather than one read per plot.
void fillStudyHistograms(TTree& tree, StudyHistograms& histograms, bool runVcbStudies) {
  const AnalysisBranches branches = analysisBranches(tree, runVcbStudies);
  enableRequiredBranches(tree, branches, runVcbStudies);

  const std::vector<PosvsNegGroup> posvsNeg = posvsNegGroups(runVcbStudies);
  const std::vector<PNetBGroup> pNetB = pNetBGroups(runVcbStudies);
  TH1D& cutflowEvents = histograms.at("cutflow_events");
  TH1D& cutflowWeight = histograms.at("cutflow_weight");
  TH1D& nonUniqueAssignment = histograms.at("meta_nonunique_assignment");
  TH1D& hadronicWMass = histograms.at("hadronic_w_mass");
  TH1D& hadronicTopMass = histograms.at("hadronic_top_mass");

  const Long64_t entries = tree.GetEntries();
  for (Long64_t entry = 0; entry < entries; ++entry) {
    tree.GetEntry(entry);
    const double weight = eventWeight(tree);
    cutflowEvents.Fill(0.5);
    cutflowWeight.Fill(0.5, weight);

    int selectedJets = 0;
    for (const auto& [index, ptBranch] : branches.cutflow.pt) {
      const double eta = leafValue(tree, branches.cutflow.eta.at(index));
      if (leafValue(tree, ptBranch) >= 20.0 && std::abs(eta) < 2.5) ++selectedJets;
    }
    if (selectedJets < 4) continue;
    cutflowEvents.Fill(1.5);
    cutflowWeight.Fill(1.5, weight);

    const std::array<int, 4> assignedJets = assignedJetIndices(tree, branches.cutflow);
    bool everyQuarkMatched = true;
    for (const int jet : assignedJets) {
      if (jet < 0 || !branches.cutflow.pt.count(jet)) everyQuarkMatched = false;
    }
    if (!everyQuarkMatched) continue;
    cutflowEvents.Fill(2.5);
    cutflowWeight.Fill(2.5, weight);

    // Distinctness is guaranteed by the upstream one-to-one assignment, so it is
    // no longer a cutflow step -- only a counter that must stay at zero.
    const std::set<int> uniqueAssignedJets(assignedJets.begin(), assignedJets.end());
    if (uniqueAssignedJets.size() != assignedJets.size()) nonUniqueAssignment.Fill(0.5);

    const int hadBJet = assignedJets[0];
    const int lepBJet = assignedJets[1];
    const int hadQ1Jet = assignedJets[2];
    const int hadQ2Jet = assignedJets[3];

    hadronicWMass.Fill(
        invariantMass(tree, branches.cutflow, branches.phi, branches.mass, {hadQ1Jet, hadQ2Jet}),
        weight);
    hadronicTopMass.Fill(invariantMass(tree, branches.cutflow, branches.phi, branches.mass,
                                       {hadQ1Jet, hadQ2Jet, hadBJet}),
                         weight);

    std::array<int, kJetRoles> jetForRole = {lepBJet, hadBJet, -1, -1};
    if (runVcbStudies) {
      jetForRole[roleIndex(JetRole::C)] = hadronicWJetByAbsPdgId(tree, hadQ1Jet, hadQ2Jet, 4);
      jetForRole[roleIndex(JetRole::WB)] = hadronicWJetByAbsPdgId(tree, hadQ1Jet, hadQ2Jet, 5);
    }

    const double leptonCharge = leafValue(tree, "TriggerLeptonCharge");
    for (const PosvsNegGroup& group : posvsNeg) {
      const int jet = jetForRole[roleIndex(group.role)];
      if (jet < 0) continue;
      const int region = ptRegionIndex(leafValue(tree, branches.cutflow.pt.at(jet)));
      if (region < 0) continue;
      if (leptonCharge == 0.0) continue;
      const double score = leafValue(tree, branches.parTPosvsNeg.at(jet));
      histograms
          .at(posvsNegHistogramName(group, ptRegions()[region], leptonCharge > 0.0))
          .Fill(score, weight);
    }

    for (const PNetBGroup& group : pNetB) {
      const int jet = jetForRole[roleIndex(group.role)];
      if (jet < 0) continue;
      const int region = ptRegionIndex(leafValue(tree, branches.cutflow.pt.at(jet)));
      if (region < 0) continue;
      histograms.at(pNetBHistogramName(group, ptRegions()[region]))
          .Fill(leafValue(tree, branches.pNetB.at(jet)), weight);
    }
  }
}

void writeStudyHistograms(const StudyHistograms& histograms, const std::filesystem::path& partial) {
  if (partial.has_parent_path() && !partial.parent_path().empty()) {
    std::filesystem::create_directories(partial.parent_path());
  }
  TFile output(partial.c_str(), "RECREATE");
  if (output.IsZombie()) throw std::runtime_error("Cannot write partial file: " + partial.string());
  for (const auto& [name, histogram] : histograms.all()) {
    output.WriteTObject(histogram.get(), name.c_str());
  }
  output.Close();
}

StudyHistograms readStudyHistograms(const std::filesystem::path& partial) {
  TFile input(partial.c_str(), "READ");
  if (input.IsZombie()) throw std::runtime_error("Cannot read partial file: " + partial.string());

  StudyHistograms histograms;
  TIter nextKey(input.GetListOfKeys());
  while (auto* key = static_cast<TKey*>(nextKey())) {
    std::unique_ptr<TObject> object(key->ReadObj());
    auto* histogram = dynamic_cast<TH1D*>(object.get());
    if (!histogram) continue;
    if (histograms.has(histogram->GetName())) continue;  // skip older key cycles
    object.release();
    histograms.adopt(std::unique_ptr<TH1D>(histogram));
  }
  if (histograms.size() == 0) {
    throw std::runtime_error("Partial file holds no histograms: " + partial.string());
  }
  return histograms;
}

std::vector<std::filesystem::path> partialFiles(const std::filesystem::path& directory) {
  if (!std::filesystem::is_directory(directory)) {
    throw std::runtime_error("Not a directory of partial files: " + directory.string());
  }
  std::vector<std::filesystem::path> partials;
  for (const auto& entry : std::filesystem::directory_iterator(directory)) {
    if (entry.is_regular_file() && entry.path().extension() == ".root") {
      partials.push_back(entry.path());
    }
  }
  if (partials.empty()) {
    throw std::runtime_error("No .root partial files found in: " + directory.string());
  }
  std::sort(partials.begin(), partials.end());
  return partials;
}

StudyHistograms mergeStudyHistograms(const std::vector<std::filesystem::path>& partials) {
  StudyHistograms total = readStudyHistograms(partials.front());
  for (std::size_t partial = 1; partial < partials.size(); ++partial) {
    total.add(readStudyHistograms(partials[partial]), partials[partial].filename().string());
  }

  const double filesWithFinalWeight = total.at("meta_has_final_weight").GetBinContent(1);
  if (filesWithFinalWeight > 0.0 &&
      filesWithFinalWeight < static_cast<double>(partials.size())) {
    std::cerr << "calib-study: warning: finalWeight was present in only "
              << filesWithFinalWeight << " of " << partials.size()
              << " partial files; the merged yields mix weighted and unweighted events.\n";
  }
  return total;
}

struct CutflowStep {
  std::string label;
  Long64_t events;
  double weightedEvents;
};

void renderCutflow(const StudyHistograms& histograms,
                   const std::string& treeName,
                   const std::filesystem::path& report) {
  const TH1D& cutflowEvents = histograms.at("cutflow_events");
  const TH1D& cutflowWeight = histograms.at("cutflow_weight");
  const bool hasFinalWeight = histograms.at("meta_has_final_weight").GetBinContent(1) > 0.0;

  std::ofstream output(report);
  if (!output) throw std::runtime_error("Cannot write report: " + report.string());
  output << "# Cutflow for TTree '" << treeName << "'\n";
  output << "# Jet requirement: PT >= 20 GeV and abs(eta) < 2.5\n";
  output << "# Assignment branches: GenHadBJetIdx, GenLepBJetIdx, GenHadQ1JetIdx, "
            "GenHadQ2JetIdx\n";
  output << "# The final step requires all four to point at a saved jet slot (index >= 0).\n";
  output << "# The assignment is one-to-one upstream, so the four jets are distinct by\n";
  output << "# construction; events violating that: "
         << std::llround(histograms.at("meta_nonunique_assignment").GetBinContent(1)) << "\n\n";
  output << "# Efficiencies use " << (hasFinalWeight ? "finalWeight" : "unit event weights") << ".\n\n";
  output << std::left << std::setw(54) << "cut" << std::right << std::setw(12) << "events"
         << std::setw(18) << "weighted yield" << std::setw(18) << "vs previous"
         << std::setw(18) << "vs initial" << '\n';
  output << std::string(120, '-') << '\n';

  std::array<CutflowStep, 3> steps;
  for (std::size_t step = 0; step < steps.size(); ++step) {
    steps[step] = {cutflowLabels()[step],
                   std::llround(cutflowEvents.GetBinContent(static_cast<int>(step) + 1)),
                   cutflowWeight.GetBinContent(static_cast<int>(step) + 1)};
  }

  output << std::fixed << std::setprecision(2);
  const double initialWeight = steps.front().weightedEvents;
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

void saveCanvas(TCanvas& canvas, const std::filesystem::path& report) {
  canvas.SaveAs(report.c_str());
  if (!std::filesystem::exists(report)) {
    throw std::runtime_error("Cannot write plot: " + report.string());
  }
}

void renderMassPlot(TH1D& histogram, const std::filesystem::path& report) {
  TCanvas canvas((std::string(histogram.GetName()) + "_canvas").c_str(), "", 800, 600);
  histogram.Draw("HIST");
  saveCanvas(canvas, report);
}

// Charge-tag working point that the purity and contamination are quoted at.
constexpr double kPosvsNegCut = 0.5;

struct CutYields {
  double negative;  // "b~ l- v~", drawn in red
  double positive;  // "b l+ v", drawn in blue
  double total;
};

// Weighted yields on the high side of the cut. 0.5 is a bin edge of the 100-bin
// [0, 1] score axis, so the cut bin opens the range; the overflow bin is kept so
// that a score of exactly 1 is not dropped.
CutYields yieldsAboveCut(const TH1D& positive, const TH1D& negative) {
  const int firstBin = positive.FindFixBin(kPosvsNegCut);
  const int overflowBin = positive.GetNbinsX() + 1;
  CutYields yields{negative.Integral(firstBin, overflowBin),
                   positive.Integral(firstBin, overflowBin), 0.0};
  yields.total = yields.negative + yields.positive;
  return yields;
}

void renderParTPosvsNegScorePlots(StudyHistograms& histograms,
                                  const PosvsNegGroup& group,
                                  const std::filesystem::path& outputDirectory) {
  const std::filesystem::path plotDirectory = outputDirectory / group.directoryName;
  std::filesystem::create_directories(plotDirectory);

  for (const PtRegion& region : ptRegions()) {
    TH1D& positive = histograms.at(posvsNegHistogramName(group, region, true));
    TH1D& negative = histograms.at(posvsNegHistogramName(group, region, false));
    positive.SetStats(0);
    negative.SetStats(0);
    positive.SetFillColorAlpha(kBlue, 0.35);
    positive.SetLineColor(kBlue + 1);
    negative.SetFillColorAlpha(kRed, 0.35);
    negative.SetLineColor(kRed + 1);

    positive.SetMinimum(0.0);
    const double tallestBin = std::max(positive.GetMaximum(), negative.GetMaximum());
    const double axisMaximum = tallestBin > 0.0 ? (4.0 / 3.0) * tallestBin : 1.0;
    positive.SetMaximum(axisMaximum);
    const CutYields aboveCut = yieldsAboveCut(positive, negative);

    TCanvas canvas((std::string(positive.GetName()) + "_canvas").c_str(), "", 800, 600);
    canvas.SetBottomMargin(0.16);
    positive.Draw("HIST");
    negative.Draw("HIST SAME");

    TLine cutLine(kPosvsNegCut, 0.0, kPosvsNegCut, axisMaximum);
    cutLine.SetLineStyle(2);
    cutLine.SetLineColor(kGray + 2);
    cutLine.Draw();

    TLegend legend(0.54, 0.73, 0.88, 0.88);
    std::ostringstream positiveLabel;
    positiveLabel << "b l+ v (weight = " << std::fixed << std::setprecision(2)
                  << positive.GetSumOfWeights() << ')';
    std::ostringstream negativeLabel;
    negativeLabel << "b~ l- v~ (weight = " << std::fixed << std::setprecision(2)
                  << negative.GetSumOfWeights() << ')';
    legend.AddEntry(&positive, positiveLabel.str().c_str(), "f");
    legend.AddEntry(&negative, negativeLabel.str().c_str(), "f");
    legend.Draw();

    // Drawn under the axis so that it never lands on top of the distributions.
    std::ostringstream cutSummary;
    cutSummary << "score > " << std::fixed << std::setprecision(1) << kPosvsNegCut << ":  ";
    if (aboveCut.total > 0.0) {
      cutSummary << std::setprecision(2) << "purity (b~ l- v~) = "
                 << 100.0 * aboveCut.negative / aboveCut.total << "%,   "
                 << "contamination (b l+ v) = "
                 << 100.0 * aboveCut.positive / aboveCut.total << '%';
    } else {
      cutSummary << "no events above the cut";
    }
    TLatex summaryText;
    summaryText.SetNDC(true);
    summaryText.SetTextSize(0.024);
    summaryText.SetTextAlign(21);
    summaryText.DrawLatex(0.5, 0.04, cutSummary.str().c_str());

    saveCanvas(canvas, plotDirectory / (region.filename + ".pdf"));
  }
}

void drawPNetBWorkingPointLabels(TCanvas& canvas) {
  const std::array<double, 5> boundaries = {0.047, 0.245, 0.6734, 0.7862, 0.961};
  const std::array<const char*, 5> labels = {"L", "M", "T", "XT", "XXT"};
  TLatex text;
  text.SetNDC(true);
  text.SetTextSize(0.03);
  const double left = canvas.GetLeftMargin();
  const double width = 1.0 - left - canvas.GetRightMargin();
  const double y = canvas.GetBottomMargin() - 0.035;
  text.SetTextAlign(21);
  for (std::size_t boundary = 0; boundary < boundaries.size(); ++boundary) {
    text.DrawLatex(left + width * boundaries[boundary], y, labels[boundary]);
  }
}

void renderPNetBPlots(StudyHistograms& histograms,
                      const PNetBGroup& group,
                      const std::filesystem::path& outputDirectory) {
  const std::filesystem::path plotDirectory = outputDirectory / group.directoryName;
  std::filesystem::create_directories(plotDirectory);

  for (const PtRegion& region : ptRegions()) {
    TH1D& histogram = histograms.at(pNetBHistogramName(group, region));
    histogram.SetStats(0);
    // Read the yields before normalizing, which rescales the bin contents.
    // GetEntries() counts every fill, including those landing outside the axis
    // range; those are neither drawn nor part of the normalization, so they are
    // reported separately instead of being folded into the quoted yield.
    const double entries = histogram.GetEntries();
    const double outsideRange =
        histogram.GetBinContent(0) + histogram.GetBinContent(histogram.GetNbinsX() + 1);
    const double integral = histogram.Integral();
    if (integral > 0.0) histogram.Scale(1.0 / integral);
    histogram.SetMinimum(0.0);
    TCanvas canvas((std::string(histogram.GetName()) + "_canvas").c_str(), "", 800, 600);
    canvas.SetBottomMargin(0.18);
    histogram.Draw("HIST");
    drawPNetBWorkingPointLabels(canvas);

    // Drawn below the working-point labels. The plot is normalized, so this is
    // the only place the absolute scale appears.
    std::ostringstream yieldSummary;
    yieldSummary << std::fixed << std::setprecision(2);
    yieldSummary << "entries = " << std::llround(entries)
                 << ",   weighted yield = " << integral;
    if (outsideRange > 0.0) {
      yieldSummary << "   (" << outsideRange << " outside the axis range, not shown)";
    }
    TLatex summaryText;
    summaryText.SetNDC(true);
    summaryText.SetTextSize(0.022);
    summaryText.SetTextAlign(21);
    summaryText.DrawLatex(0.5, 0.04, yieldSummary.str().c_str());

    saveCanvas(canvas, plotDirectory / (region.filename + ".pdf"));
  }
}

void renderStudies(StudyHistograms& histograms,
                   const std::filesystem::path& outputDirectory,
                   const std::string& treeName,
                   bool runVcbStudies) {
  renderCutflow(histograms, treeName, outputDirectory / "cutflow.txt");
  renderMassPlot(histograms.at("hadronic_w_mass"), outputDirectory / "hadronic-w-mass.pdf");
  renderMassPlot(histograms.at("hadronic_top_mass"), outputDirectory / "hadronic-top-mass.pdf");
  for (const PosvsNegGroup& group : posvsNegGroups(runVcbStudies)) {
    renderParTPosvsNegScorePlots(histograms, group, outputDirectory);
  }
  for (const PNetBGroup& group : pNetBGroups(runVcbStudies)) {
    renderPNetBPlots(histograms, group, outputDirectory);
  }
}

void addRootFiles(TChain& chain,
                  const std::vector<std::string>& rootFiles,
                  const std::string& treeName) {
  for (const auto& rootFile : rootFiles) {
    if (chain.Add(rootFile.c_str()) == 0) {
      throw std::runtime_error("Cannot add TTree '" + treeName + "' from: " + rootFile);
    }
  }
}

enum class Mode { All, Partial, Merge };

void usage(const char* program) {
  std::cerr << "Usage: " << program
            << " <input-address-txt> <production-tag> [tree-name] [--Vcb]\n"
            << "                          [--partial <file.root> | --merge <partials-dir>]\n\n"
            << "  (no mode flag)          run every study and write the reports and plots\n"
            << "  --partial <file.root>   fill the study histograms only and write them to\n"
            << "                          <file.root>; no reports or plots are produced\n"
            << "  --merge <partials-dir>  add every .root partial in <partials-dir> and write\n"
            << "                          the reports and plots for the summed histograms\n";
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc < 3) {
    usage(argv[0]);
    return 2;
  }

  try {
    const std::string inputList = argv[1];
    const std::string productionTag = argv[2];
    std::string treeName = "Events";
    bool treeNameSet = false;
    bool runVcbStudies = false;
    Mode mode = Mode::All;
    std::filesystem::path modePath;
    for (int argument = 3; argument < argc; ++argument) {
      const std::string value = argv[argument];
      if (value == "--Vcb") {
        if (runVcbStudies) throw std::runtime_error("--Vcb may be specified only once");
        runVcbStudies = true;
      } else if (value == "--partial" || value == "--merge") {
        if (mode != Mode::All) {
          throw std::runtime_error("--partial and --merge may not be combined or repeated");
        }
        if (argument + 1 >= argc) throw std::runtime_error("Missing path after " + value);
        mode = value == "--partial" ? Mode::Partial : Mode::Merge;
        modePath = argv[++argument];
      } else if (value.rfind("--", 0) == 0) {
        throw std::runtime_error("Unknown option: " + value);
      } else if (treeNameSet) {
        throw std::runtime_error("Only one optional tree-name may be specified");
      } else {
        treeName = value;
        treeNameSet = true;
      }
    }
    if (!validProductionTag(productionTag)) {
      throw std::runtime_error("production-tag must be one safe path component (no '/' or '\\\\')");
    }

    const auto rootFiles = inputRootFiles(inputList);

    if (mode == Mode::Partial) {
      // --merge only collects .root files, so a partial written under any other
      // extension would be skipped there without a word.
      if (modePath.extension() != ".root") {
        throw std::runtime_error("--partial must name a .root file, not: " + modePath.string());
      }
      TChain tree(treeName.c_str());
      addRootFiles(tree, rootFiles, treeName);
      if (tree.LoadTree(0) < 0) {
        throw std::runtime_error("Cannot load the first entry from the input chain");
      }
      StudyHistograms histograms = bookStudyHistograms(
          productionTag, runVcbStudies, tree.GetLeaf("finalWeight") != nullptr);
      fillStudyHistograms(tree, histograms, runVcbStudies);
      writeStudyHistograms(histograms, modePath);
      std::cout << "Filled study histograms from " << rootFiles.size()
                << " ROOT files and wrote the partial to " << modePath << '\n';
      return 0;
    }

    const std::filesystem::path outputDirectory =
        std::filesystem::path("outputs") / productionTag;
    std::filesystem::create_directories(outputDirectory);

    // Studies 1 and 2 describe the schema and entry 0, which every input file in
    // a production shares, so only the first file needs to be opened for them.
    TChain schemaChain(treeName.c_str());
    addRootFiles(schemaChain, {rootFiles.front()}, treeName);
    if (schemaChain.LoadTree(0) < 0) {
      throw std::runtime_error("Cannot load the first entry from the input chain");
    }
    writeFields(schemaChain, outputDirectory / "fields.txt");
    writeFirstEvent(schemaChain, outputDirectory / "0th-event.txt");

    std::size_t mergedPartials = 0;
    std::unique_ptr<StudyHistograms> histograms;
    if (mode == Mode::Merge) {
      const auto partials = partialFiles(modePath);
      mergedPartials = partials.size();
      histograms = std::make_unique<StudyHistograms>(mergeStudyHistograms(partials));
      // The opposite mismatch stops the run with a missing-histogram error, but
      // this one would quietly drop studies 9-11 from the output.
      if (!runVcbStudies &&
          histograms->has(posvsNegHistogramName(posvsNegGroups(true).back(), ptRegions().front(),
                                                true))) {
        std::cerr << "calib-study: warning: the partial files hold the --Vcb studies but --Vcb "
                     "was not passed to the merge, so studies 9-11 will not be plotted.\n";
      }
    } else {
      TChain tree(treeName.c_str());
      addRootFiles(tree, rootFiles, treeName);
      if (tree.LoadTree(0) < 0) {
        throw std::runtime_error("Cannot load the first entry from the input chain");
      }
      histograms = std::make_unique<StudyHistograms>(bookStudyHistograms(
          productionTag, runVcbStudies, tree.GetLeaf("finalWeight") != nullptr));
      fillStudyHistograms(tree, *histograms, runVcbStudies);
    }

    renderStudies(*histograms, outputDirectory, treeName, runVcbStudies);
    if (mode == Mode::Merge) {
      std::cout << "Merged " << mergedPartials << " partial files and wrote reports to "
                << outputDirectory << '\n';
    } else {
      std::cout << "Processed " << rootFiles.size() << " ROOT files and wrote reports to "
                << outputDirectory << '\n';
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "calib-study: " << error.what() << '\n';
    return 1;
  }
}
