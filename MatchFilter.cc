// Keep the MC events whose four ttbar decay partons each match a distinct jet,
// write them to a slim ROOT file, and report how the matching went.
//
//   ./match-filter <input.root | input-list.txt> <output.root>
//
// Input: Calib_ChargeTagger_JH skims. Upstream already requires HLT_IsoMu24 or
// HLT_Ele30_WPTight_Gsf with a trigger-matched lepton (muon pT >= 26, |eta| < 2.4;
// electron pT >= 32, |eta| < 2.5), and saves up to ten Tight-ID jets with pT > 15
// and |eta| < 4.7, cleaned against that lepton. This program then requires:
//   1. at least four saved jets with pT > 25 GeV and |eta| < 2.4 (analysis jets);
//   2. hadB, lepB, q1 and q2 each greedily matched to a distinct analysis jet.
// Events passing both are written to <output.root>; <output>-report.txt breaks
// down the matching of every event that passes step 1.

#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TTree.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef GIT_COMMIT
#define GIT_COMMIT "unknown"
#endif

namespace {

constexpr double kJetPtMin = 25.0;
constexpr double kJetAbsEtaMax = 2.4;
constexpr double kMaxDeltaR = 0.4;
constexpr double kPad = -99999.0;  // upstream's value for an empty slot or a missing object
constexpr double kPi = 3.14159265358979323846;
constexpr int kJetSlots = 10;

// The trigger paths, each named by the |pdgId| of the lepton it selects, so that
// |TriggerLepton_pdgId| is the path that selected the event.
constexpr int kIsoMu24 = 13;  // HLT_IsoMu24
constexpr int kEle30 = 11;    // HLT_Ele30_WPTight_Gsf

// Parton order, which is also the matching tie-break order. Upstream stores the
// hadronic-W daughters by type: GenQ1 is always down-type, GenQ2 always up-type.
constexpr int kPartons = 4;
constexpr unsigned kAllMatched = (1u << kPartons) - 1;
const std::array<std::string, kPartons> kPartonLabels = {"hadB", "lepB", "q1", "q2"};
const std::array<std::string, kPartons> kPartonInputs = {"GenHadTopB", "GenLepTopB", "GenQ1",
                                                         "GenQ2"};
const std::array<std::string, kPartons> kPartonOutputs = {"GenHadB", "GenLepB", "GenQ1",
                                                          "GenQ2"};

// Upstream's six one-hot W-decay flags, stored as one code in CKM order (V_ud ... V_cb):
// 10 * |up-type pdgId| + |down-type pdgId|, i.e. 10 * |GenQ2PdgId| + |GenQ1PdgId|.
// W -> cb, upstream's GenWtoBC, is 45.
constexpr int kWDecays = 6;
const std::array<std::string, kWDecays> kWDecayFlags = {"GenWtoUD", "GenWtoCD", "GenWtoUS",
                                                        "GenWtoCS", "GenWtoBC", "GenWtoUB"};
constexpr std::array<int, kWDecays> kWDecayCodes = {21, 41, 23, 43, 45, 25};

// Four-vector components, in the order the arrays below store them.
enum { kPt, kEta, kPhi, kMass };
const std::array<std::string, 4> kInputComponents = {"Pt", "Eta", "Phi", "Mass"};
const std::array<std::string, 4> kOutputComponents = {"pt", "eta", "phi", "mass"};

// Jet fields: the four-vector, then the b-tag and charge-tag scores.
const std::array<std::string, 6> kJetInputs = {"ak4JetPt", "ak4JetEta", "ak4JetPhi",
                                               "ak4JetMass", "ak4JetbtagPNetB",
                                               "ak4JetParTPosvsNeg"};
const std::array<std::string, 6> kJetOutputs = {"Jet_pt", "Jet_eta", "Jet_phi",
                                                "Jet_mass", "Jet_btagPNetB", "Jet_ParTPosvsNeg"};

// The upstream branches this program reads, with their stored types.
struct InputEvent {
  UInt_t run = 0, luminosityBlock = 0;
  ULong64_t event = 0;
  Double_t weight = 0;
  Long64_t isoMu24 = 0, ele30 = 0, leptonFlavor = 0, leptonCharge = 0;
  std::array<Double_t, 3> lepton{};  // the trigger lepton's pt, eta, phi
  Double_t metPt = 0, metPhi = 0;
  std::array<std::array<Double_t, kJetSlots>, 6> jet{};
  std::array<std::array<Float_t, kPartons>, 4> parton{};
  Int_t q1PdgId = 0, q2PdgId = 0;
  std::array<Bool_t, kWDecays> wDecay{};
  std::array<Float_t, 3> genLepton{}, neutrino{};  // pt, eta, phi
  Int_t genLeptonCharge = 0;
  UInt_t genLeptonFlavor = 0;
};

// One row of the output tree.
struct SlimEvent {
  UInt_t run = 0, luminosityBlock = 0;
  ULong64_t event = 0;
  Double_t weight = 0;
  std::array<Float_t, 3> lepton{};  // the trigger lepton's pt, eta, phi
  Int_t leptonPdgId = 0;
  Float_t metPt = 0, metPhi = 0;
  Int_t nJet = 0;
  std::array<std::array<Float_t, kJetSlots>, 6> jet{};
  std::array<std::array<Float_t, kPartons>, 4> parton{};
  std::array<Int_t, kPartons> partonJet{};
  Int_t q1PdgId = 0, q2PdgId = 0, wDecay = 0;
  std::array<Float_t, 3> genLepton{}, neutrino{};
  Int_t genLeptonPdgId = 0;
};

std::vector<std::string> inputFiles(const std::string& argument) {
  if (std::filesystem::path(argument).extension() == ".root") return {argument};
  std::ifstream list(argument);
  if (!list) throw std::runtime_error("Cannot open input list: " + argument);
  std::vector<std::string> files;
  std::set<std::string> seen;
  for (std::string line; std::getline(list, line);) {
    line.erase(0, line.find_first_not_of(" \t\r"));
    line.erase(line.find_last_not_of(" \t\r") + 1);
    if (line.empty() || line.front() == '#') continue;
    if (!seen.insert(line).second) throw std::runtime_error("Duplicate input: " + line);
    files.push_back(line);
  }
  if (files.empty()) throw std::runtime_error("No input files in " + argument);
  return files;
}

// Opens every input read-only, chains it, and returns its finalWeight
// normalization record, which must be the same for every input.
std::string addInputs(TChain& chain, const std::vector<std::string>& files) {
  std::string normalization;
  for (const auto& file : files) {
    TFile input(file.c_str(), "READ");
    TTree* tree = nullptr;
    if (!input.IsZombie()) input.GetObject("Events", tree);
    if (!tree) throw std::runtime_error("Cannot read the Events tree from " + file);
    const auto* record = input.Get<TNamed>("finalWeight_provenance");
    if (!record) throw std::runtime_error("Missing finalWeight_provenance in " + file);
    if (normalization.empty()) normalization = record->GetTitle();
    if (normalization != record->GetTitle()) {
      throw std::runtime_error("Inputs were normalized separately: " + file);
    }
    if (chain.Add(file.c_str(), tree->GetEntries()) == 0) {
      throw std::runtime_error("Cannot chain " + file);
    }
  }
  return normalization;
}

template <typename T>
void bindBranch(TChain& chain, const std::string& name, T* address) {
  chain.SetBranchStatus(name.c_str(), true);
  if (chain.SetBranchAddress(name.c_str(), address) < 0) {
    throw std::runtime_error("Cannot read " + name + " with the expected type");
  }
}

void bindInput(TChain& chain, InputEvent& in) {
  chain.SetBranchStatus("*", false);
  bindBranch(chain, "run", &in.run);
  bindBranch(chain, "luminosityBlock", &in.luminosityBlock);
  bindBranch(chain, "event", &in.event);
  bindBranch(chain, "finalWeight", &in.weight);
  bindBranch(chain, "HLT_IsoMu24", &in.isoMu24);
  bindBranch(chain, "HLT_Ele30_WPTight_Gsf", &in.ele30);
  bindBranch(chain, "TriggerLeptonFlav", &in.leptonFlavor);
  bindBranch(chain, "TriggerLeptonCharge", &in.leptonCharge);
  bindBranch(chain, "METPt", &in.metPt);
  bindBranch(chain, "METPhi", &in.metPhi);
  bindBranch(chain, "GenQ1PdgId", &in.q1PdgId);
  bindBranch(chain, "GenQ2PdgId", &in.q2PdgId);
  for (int mode = 0; mode < kWDecays; ++mode) bindBranch(chain, kWDecayFlags[mode], &in.wDecay[mode]);
  bindBranch(chain, "GenWLeptonCharge", &in.genLeptonCharge);
  bindBranch(chain, "GenWLeptonFlavor", &in.genLeptonFlavor);
  // Their masses are not read: the lepton's is the PDG value for its flavour, and
  // upstream sets the neutrino's to zero.
  for (int c = 0; c < 3; ++c) {
    const std::string& component = kInputComponents[c];
    bindBranch(chain, "TriggerLepton" + component, &in.lepton[c]);
    bindBranch(chain, "GenWLepton" + component, &in.genLepton[c]);
    bindBranch(chain, "GenWNeutrino" + component, &in.neutrino[c]);
  }
  for (int c = 0; c < 4; ++c) {
    const std::string& component = kInputComponents[c];
    for (int p = 0; p < kPartons; ++p) bindBranch(chain, kPartonInputs[p] + component, &in.parton[c][p]);
  }
  for (std::size_t field = 0; field < kJetInputs.size(); ++field) {
    for (int slot = 0; slot < kJetSlots; ++slot) {
      bindBranch(chain, kJetInputs[field] + std::to_string(slot), &in.jet[field][slot]);
    }
  }
}

void bookOutput(TTree& tree, SlimEvent& out) {
  const auto scalarBranch = [&tree](const std::string& name, void* address, const char* type) {
    tree.Branch(name.c_str(), address, (name + '/' + type).c_str());
  };
  const auto arrayBranch = [&tree](const std::string& name, void* address, const std::string& size,
                             const char* type) {
    tree.Branch(name.c_str(), address, (name + '[' + size + "]/" + type).c_str());
  };
  scalarBranch("run", &out.run, "i");
  scalarBranch("luminosityBlock", &out.luminosityBlock, "i");
  scalarBranch("event", &out.event, "l");
  scalarBranch("finalWeight", &out.weight, "D");
  for (int c = 0; c < 3; ++c) {
    scalarBranch("TriggerLepton_" + kOutputComponents[c], &out.lepton[c], "F");
  }
  scalarBranch("TriggerLepton_pdgId", &out.leptonPdgId, "I");
  scalarBranch("MET_pt", &out.metPt, "F");
  scalarBranch("MET_phi", &out.metPhi, "F");
  scalarBranch("nJet", &out.nJet, "I");
  for (std::size_t field = 0; field < kJetOutputs.size(); ++field) {
    arrayBranch(kJetOutputs[field], out.jet[field].data(), "nJet", "F");
  }
  for (int p = 0; p < kPartons; ++p) {
    for (int c = 0; c < 4; ++c) {
      scalarBranch(kPartonOutputs[p] + '_' + kOutputComponents[c], &out.parton[c][p], "F");
    }
    scalarBranch(kPartonOutputs[p] + "_jetIdx", &out.partonJet[p], "I");
  }
  scalarBranch("GenQ1_pdgId", &out.q1PdgId, "I");
  scalarBranch("GenQ2_pdgId", &out.q2PdgId, "I");
  scalarBranch("GenW_decay", &out.wDecay, "I");
  for (int c = 0; c < 3; ++c) {
    scalarBranch("GenLepton_" + kOutputComponents[c], &out.genLepton[c], "F");
    scalarBranch("GenNeutrino_" + kOutputComponents[c], &out.neutrino[c], "F");
  }
  scalarBranch("GenLepton_pdgId", &out.genLeptonPdgId, "I");
}

void requireFinite(double value, const char* what, Long64_t entry) {
  if (!std::isfinite(value)) {
    throw std::runtime_error(std::string("Non-finite ") + what + " at entry " +
                             std::to_string(entry));
  }
}

// Upstream fills slots from the front, so a filled slot after an empty one means
// the input is not the layout this program expects.
template <std::size_t N>
int filledSlots(const std::array<Double_t, N>& pt, const char* collection, Long64_t entry) {
  int filled = 0;
  while (filled < static_cast<int>(N) && pt[filled] != kPad) ++filled;
  for (std::size_t slot = filled; slot < N; ++slot) {
    if (pt[slot] != kPad) {
      throw std::runtime_error(std::string("Gap in the saved ") + collection + " slots at entry " +
                               std::to_string(entry));
    }
  }
  return filled;
}

// (a - b + pi) mod 2pi - pi, with the non-negative modulus numpy's % gives, so the
// result is bit-identical to upstream's matching.
double deltaPhi(double a, double b) {
  double shifted = std::fmod(a - b + kPi, 2.0 * kPi);
  if (shifted < 0.0) shifted += 2.0 * kPi;
  return shifted - kPi;
}

// Greedy one-to-one assignment, the rule of Calib_ChargeTagger_JH's
// _greedy_jet_parton_assignment: take the smallest DeltaR < kMaxDeltaR among the
// remaining parton-jet pairs, retire both, and repeat. Exact ties go to the earlier
// parton, then the harder jet. Returns, per parton, a position in `jets`, or -1.
std::array<int, kPartons> greedyMatch(const InputEvent& in, const std::vector<int>& jets) {
  constexpr double kTaken = std::numeric_limits<double>::infinity();
  const int nJets = static_cast<int>(jets.size());
  std::vector<double> deltaR(kPartons * nJets, kTaken);
  for (int p = 0; p < kPartons; ++p) {
    if (in.parton[kEta][p] == static_cast<Float_t>(kPad)) continue;  // parton absent
    for (int k = 0; k < nJets; ++k) {
      const double dEta = static_cast<double>(in.parton[kEta][p]) - in.jet[kEta][jets[k]];
      const double dPhi = deltaPhi(in.parton[kPhi][p], in.jet[kPhi][jets[k]]);
      const double distance = std::sqrt(dEta * dEta + dPhi * dPhi);
      if (distance < kMaxDeltaR) deltaR[p * nJets + k] = distance;
    }
  }

  std::array<int, kPartons> match;
  match.fill(-1);
  for (int round = 0; round < kPartons; ++round) {
    int best = -1;
    for (int cell = 0; cell < kPartons * nJets; ++cell) {
      if (deltaR[cell] < kTaken && (best < 0 || deltaR[cell] < deltaR[best])) best = cell;
    }
    if (best < 0) break;
    const int parton = best / nJets;
    const int jet = best % nJets;
    match[parton] = jet;
    for (int k = 0; k < nJets; ++k) deltaR[parton * nJets + k] = kTaken;
    for (int p = 0; p < kPartons; ++p) deltaR[p * nJets + jet] = kTaken;
  }
  return match;
}

struct Tally {
  Long64_t events = 0;
  double weight = 0.0;
  void add(double w) {
    ++events;
    weight += w;
  }
};

std::string partonSet(unsigned mask) {
  std::string set = "{";
  for (int p = 0; p < kPartons; ++p) {
    if (!(mask & (1u << p))) continue;
    if (set.size() > 1) set += ", ";
    set += kPartonLabels[p];
  }
  return set + '}';
}

void writeReport(const std::filesystem::path& report, const std::string& inputArgument,
                 std::size_t nInputs, const std::array<Tally, 3>& cutflow,
                 const std::array<Tally, kAllMatched + 1>& matched) {
  std::ofstream out(report);
  if (!out) throw std::runtime_error("Cannot write report: " + report.string());
  const auto percent = [](double part, double whole) { return whole == 0 ? 0.0 : 100 * part / whole; };
  const auto row = [&out, &percent](const std::string& label, const Tally& tally,
                                    const Tally& whole) {
    out << std::left << std::setw(48) << label << std::right << std::setw(10) << tally.events
        << std::setw(10) << std::setprecision(2) << percent(tally.events, whole.events)
        << std::setw(16) << std::setprecision(4) << tally.weight << std::setw(10)
        << std::setprecision(2) << percent(tally.weight, whole.weight) << '\n';
  };
  const auto header = [&out](const std::string& first, const std::string& share) {
    out << std::left << std::setw(48) << first << std::right << std::setw(10) << "events"
        << std::setw(10) << share << std::setw(16) << "weighted" << std::setw(10) << share << '\n'
        << std::string(94, '-') << '\n';
  };

  out << std::fixed;
  out << "# match-filter report (commit " << GIT_COMMIT << ")\n"
      << "# Input: " << inputArgument << " (" << nInputs << " file(s))\n"
      << "# Upstream already applied the trigger, the trigger-lepton cuts, the Tight\n"
      << "#   jet ID, and jet cleaning against the trigger lepton.\n"
      << "# Analysis jets: saved jets with pT > 25 GeV and |eta| < 2.4.\n"
      << "# Matching: greedy one-to-one among analysis jets, smallest DeltaR first,\n"
      << "#   DeltaR < 0.4.\n"
      << "# Partons: hadB, lepB, q1 = down-type W daughter (d, s, b), q2 = up-type (u, c).\n"
      << "# Weighted columns sum finalWeight, which is signed.\n\n";

  header("Cutflow", "%");
  row("Input", cutflow[0], cutflow[0]);
  row("At least 4 analysis jets", cutflow[1], cutflow[0]);
  row("All four partons matched (written)", cutflow[2], cutflow[0]);

  out << "\nMatching among events with at least 4 analysis jets; % of those events\n";
  header("Matched partons", "%");
  for (int n = kPartons; n >= 0; --n) {
    Tally group;
    for (unsigned mask = 0; mask <= kAllMatched; ++mask) {
      if (__builtin_popcount(mask) != n) continue;
      group.events += matched[mask].events;
      group.weight += matched[mask].weight;
    }
    row(std::to_string(n) + "/4 matched", group, cutflow[1]);
    if (n == 0 || n == kPartons) continue;
    for (int mask = kAllMatched; mask >= 0; --mask) {
      if (__builtin_popcount(mask) == n) row("    " + partonSet(mask), matched[mask], cutflow[1]);
    }
  }
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <input.root | input-list.txt> <output.root>\n";
    return 2;
  }

  try {
    const std::vector<std::string> inputs = inputFiles(argv[1]);
    const std::filesystem::path outputPath = argv[2];
    if (outputPath.extension() != ".root") {
      throw std::runtime_error("The output must be a .root file: " + outputPath.string());
    }
    for (const auto& input : inputs) {
      if (input.find("://") == std::string::npos && std::filesystem::exists(outputPath) &&
          std::filesystem::equivalent(input, outputPath)) {
        throw std::runtime_error("Refusing to overwrite an input: " + input);
      }
    }
    std::filesystem::path reportPath = outputPath;
    reportPath.replace_filename(outputPath.stem().string() + "-report.txt");

    TChain chain("Events");
    const std::string normalization = addInputs(chain, inputs);
    if (chain.LoadTree(0) < 0) throw std::runtime_error("The inputs contain no events");
    InputEvent in;
    bindInput(chain, in);

    if (outputPath.has_parent_path()) std::filesystem::create_directories(outputPath.parent_path());
    TFile output(outputPath.c_str(), "RECREATE");
    if (output.IsZombie()) throw std::runtime_error("Cannot write " + outputPath.string());
    auto* tree = new TTree("Events", "Fully matched events; see match_filter_provenance");
    SlimEvent out;
    bookOutput(*tree, out);

    std::array<Tally, 3> cutflow;
    std::array<Tally, kAllMatched + 1> matched;
    const Long64_t entries = chain.GetEntries();
    for (Long64_t entry = 0; entry < entries; ++entry) {
      if (chain.GetEntry(entry) <= 0) {
        throw std::runtime_error("Cannot read entry " + std::to_string(entry));
      }
      requireFinite(in.weight, "finalWeight", entry);
      cutflow[0].add(in.weight);

      // Upstream already selected on the trigger. TriggerLepton_pdgId stands for the
      // path, so stop if the trigger lepton does not belong to one that fired.
      const bool pathFired = (in.leptonFlavor == kIsoMu24 && in.isoMu24 == 1) ||
                             (in.leptonFlavor == kEle30 && in.ele30 == 1);
      if (!pathFired || std::abs(in.leptonCharge) != 1) {
        throw std::runtime_error("Unexpected trigger lepton at entry " + std::to_string(entry));
      }

      std::vector<int> jets;
      const int savedJets = filledSlots(in.jet[kPt], "jet", entry);
      for (int slot = 0; slot < savedJets; ++slot) {
        requireFinite(in.jet[kPt][slot], "jet pT", entry);
        requireFinite(in.jet[kEta][slot], "jet eta", entry);
        requireFinite(in.jet[kPhi][slot], "jet phi", entry);
        if (in.jet[kPt][slot] > kJetPtMin && std::abs(in.jet[kEta][slot]) < kJetAbsEtaMax) {
          jets.push_back(slot);
        }
      }
      if (jets.size() < 4) continue;
      cutflow[1].add(in.weight);

      for (int p = 0; p < kPartons; ++p) {
        requireFinite(in.parton[kEta][p], "parton eta", entry);
        requireFinite(in.parton[kPhi][p], "parton phi", entry);
      }
      const std::array<int, kPartons> match = greedyMatch(in, jets);
      unsigned mask = 0;
      for (int p = 0; p < kPartons; ++p) {
        if (match[p] >= 0) mask |= 1u << p;
      }
      matched[mask].add(in.weight);
      if (mask != kAllMatched) continue;

      const int q1 = std::abs(in.q1PdgId);
      const int q2 = std::abs(in.q2PdgId);
      if (!(q1 == 1 || q1 == 3 || q1 == 5) || !(q2 == 2 || q2 == 4)) {
        throw std::runtime_error("GenQ1 is not down-type or GenQ2 is not up-type at entry " +
                                 std::to_string(entry));
      }
      int wDecay = 0;
      for (int mode = 0; mode < kWDecays; ++mode) {
        if (!in.wDecay[mode]) continue;
        if (wDecay != 0) {
          throw std::runtime_error("Several W-decay flags set at entry " + std::to_string(entry));
        }
        wDecay = kWDecayCodes[mode];
      }
      if (wDecay != 10 * q2 + q1) {
        throw std::runtime_error("The W-decay flags disagree with GenQ1PdgId and GenQ2PdgId at entry " +
                                 std::to_string(entry));
      }
      cutflow[2].add(in.weight);

      out.run = in.run;
      out.luminosityBlock = in.luminosityBlock;
      out.event = in.event;
      out.weight = in.weight;
      // The trigger lepton only; its flavour gives its mass. The PDG ID is negative
      // for a positive lepton (mu+ is -13), as for GenLepton_pdgId.
      for (int c = 0; c < 3; ++c) out.lepton[c] = static_cast<Float_t>(in.lepton[c]);
      out.leptonPdgId = static_cast<Int_t>(-in.leptonCharge * in.leptonFlavor);
      out.metPt = static_cast<Float_t>(in.metPt);
      out.metPhi = static_cast<Float_t>(in.metPhi);
      out.nJet = static_cast<Int_t>(jets.size());
      for (std::size_t field = 0; field < kJetInputs.size(); ++field) {
        for (std::size_t k = 0; k < jets.size(); ++k) {
          out.jet[field][k] = static_cast<Float_t>(in.jet[field][jets[k]]);
        }
      }
      out.parton = in.parton;
      out.partonJet = match;
      out.q1PdgId = in.q1PdgId;
      out.q2PdgId = in.q2PdgId;
      out.wDecay = wDecay;
      out.genLepton = in.genLepton;
      out.neutrino = in.neutrino;
      // Upstream stores the W's charged lepton as |pdgId| and charge; this is its
      // PDG ID, negative for a positive lepton (e+ is -11).
      const int flavor = static_cast<int>(in.genLeptonFlavor);
      if (!(flavor == 11 || flavor == 13 || flavor == 15) || std::abs(in.genLeptonCharge) != 1) {
        throw std::runtime_error("Unexpected W-decay lepton at entry " + std::to_string(entry));
      }
      out.genLeptonPdgId = -in.genLeptonCharge * flavor;
      tree->Fill();
    }

    auto* cutflowEvents = new TH1D("cutflow_events", "Cutflow (events)", 3, 0, 3);
    auto* cutflowWeight = new TH1D("cutflow_weight", "Cutflow (sum of finalWeight)", 3, 0, 3);
    const std::array<const char*, 3> steps = {"Input", "At least 4 analysis jets",
                                              "All four partons matched"};
    for (int step = 0; step < 3; ++step) {
      for (TH1D* histogram : {cutflowEvents, cutflowWeight}) {
        histogram->GetXaxis()->SetBinLabel(step + 1, steps[step]);
      }
      cutflowEvents->SetBinContent(step + 1, static_cast<double>(cutflow[step].events));
      cutflowWeight->SetBinContent(step + 1, cutflow[step].weight);
    }

    std::ostringstream provenance;
    provenance << "match-filter commit " << GIT_COMMIT << "\n"
               << "analysis jets: saved jets with pT > " << kJetPtMin << " GeV, |eta| < "
               << kJetAbsEtaMax << "; events need at least 4\n"
               << "matching: greedy one-to-one, smallest DeltaR first, DeltaR < " << kMaxDeltaR
               << ", among analysis jets; order hadB, lepB, q1, q2\n"
               << "q1 = down-type W daughter, q2 = up-type\n"
               << "TriggerLepton_pdgId: PDG ID of the trigger lepton (mu+ = -13), which gives "
                  "its mass; |value| is the path that selected the event: 13 = HLT_IsoMu24, "
                  "11 = HLT_Ele30_WPTight_Gsf\n"
               << "GenLepton_pdgId: PDG ID of the W's charged lepton (e- = 11, e+ = -11); its "
                  "mass is the PDG mass; GenNeutrino is massless\n"
               << "GenW_decay: 10 * |up-type pdgId| + |down-type pdgId| (CKM order): "
                  "21 = ud, 23 = us, 25 = ub, 41 = cd, 43 = cs, 45 = cb (upstream GenWtoBC)\n"
               << "inputs:\n";
    for (const auto& input : inputs) provenance << input << '\n';

    output.cd();
    tree->Write();
    cutflowEvents->Write();
    cutflowWeight->Write();
    TNamed("finalWeight_provenance", normalization.c_str()).Write();
    TNamed("match_filter_provenance", provenance.str().c_str()).Write();
    output.Close();

    writeReport(reportPath, argv[1], inputs.size(), cutflow, matched);
    std::cout << "Wrote " << cutflow[2].events << " of " << cutflow[0].events << " events to "
              << outputPath << " and the report to " << reportPath << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "match-filter: " << error.what() << '\n';
    return 1;
  }
}
