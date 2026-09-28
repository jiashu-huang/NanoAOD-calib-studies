"""End-to-end checks with small ROOT fixtures. Run in the CMSSW ROOT environment."""

from array import array
from pathlib import Path
import subprocess
import tempfile
import unittest

import ROOT

ROOT.gROOT.SetBatch(True)
EXECUTABLE = Path(__file__).resolve().parents[1] / "calib-study"
W_DECAYS = {"UD": (1, 2), "CD": (1, 4), "US": (2, 3),
            "CS": (3, 4), "BC": (4, 5), "UB": (2, 5)}


class CalibStudyTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)

    def fixture(self, name, events, weighted=True, missing=None):
        values = {
            "GenHadBJetIdx": 0, "GenLepBJetIdx": 1,
            "GenHadQ1JetIdx": 2, "GenHadQ2JetIdx": 3,
            "GenQ1PdgId": 5, "GenQ2PdgId": -4, "TriggerLeptonCharge": 1,
        }
        values.update({"GenWto" + decay: int(decay == "BC") for decay in W_DECAYS})
        for slot, pt in enumerate([40, 45, 50, 55, 17, 21]):
            for stem, value in {"Pt": pt, "Eta": .2, "Phi": .5 * slot,
                                "Mass": 5, "ParTPosvsNeg": .5, "btagPNetB": .5}.items():
                values[f"ak4Jet{stem}{slot}"] = value
        if weighted:
            values["finalWeight"] = 1
        if missing:
            del values[missing]
        path = self.directory / name
        output = ROOT.TFile(str(path), "RECREATE")
        tree = ROOT.TTree("Events", "Events")
        buffers = {key: array("d", [value]) for key, value in values.items()}
        for key, buffer in buffers.items():
            tree.Branch(key, buffer, key + "/D")
        for event in events:
            daughters = tuple(sorted(abs(event.get(key, values[key]))
                                     for key in ["GenQ1PdgId", "GenQ2PdgId"]))
            flags = {"GenWto" + decay: int(daughters == pair)
                     for decay, pair in W_DECAYS.items()}
            for key, buffer in buffers.items():
                buffer[0] = event.get(key, flags.get(key, values[key]))
            tree.Fill()
        tree.Write()
        output.Close()
        return path

    def run_study(self, files, tag="test", extra=(), success=True):
        addresses = self.directory / "inputs.txt"
        addresses.write_text("".join(str(path) + "\n" for path in files))
        result = subprocess.run(
            [str(EXECUTABLE), str(addresses), tag, "--Vcb", *map(str, extra)],
            cwd=self.directory, capture_output=True, text=True,
        )
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def partial(self, files, name="one.root", success=True):
        path = self.directory / "partials" / name
        result = self.run_study(files, extra=["--partial", path], success=success)
        return path, result

    def test_cutflow_endpoints_and_signed_weights(self):
        path = self.fixture("events.root", [
            {"finalWeight": 2, "ak4JetbtagPNetB1": 1, "ak4JetParTPosvsNeg1": 1},
            {"finalWeight": -.5, "TriggerLeptonCharge": -1},
            {"finalWeight": 1, "ak4JetPt2": 17, "ak4JetPt3": 17},
            {"finalWeight": 3, "GenHadQ2JetIdx": -1},
        ])
        partial, _ = self.partial([path])
        with ROOT.TFile(str(partial)) as output:
            self.assertEqual(list(output.Get("cutflow_events"))[1:4], [4, 3, 2])
            self.assertEqual(list(output.Get("cutflow_weight"))[1:4], [5.5, 4.5, 1.5])
            pnet = output.Get("pnetb_LepB_30to50")
            self.assertEqual(pnet.GetBinContent(6), 2)
            self.assertEqual(pnet.GetBinContent(7), 0)
            part = output.Get("lepb_part_posvsneg_30to50_positive")
            self.assertEqual(part.GetBinContent(100), 2)
            self.assertEqual(part.GetBinContent(101), 0)

    def test_mass_selection_keeps_matched_jet_below_20(self):
        path = self.fixture("events.root", [{"GenHadQ2JetIdx": 4}])
        partial, _ = self.partial([path])
        with ROOT.TFile(str(partial)) as output:
            self.assertEqual(output.Get("hadronic_w_mass").GetEntries(), 1)
            self.assertEqual(sum(output.Get("pnetb_UpType_" + region).GetEntries()
                                 for region in ["20to30", "30to50", "50to80", "80to110",
                                                "110to150", "150plus"]), 0)

    def test_up_type_includes_u_and_c(self):
        path = self.fixture("events.root", [{"GenQ2PdgId": -2}, {"GenQ2PdgId": -4}])
        partial, _ = self.partial([path])
        with ROOT.TFile(str(partial)) as output:
            self.assertEqual(output.Get("pnetb_UpType_50to80").GetEntries(), 2)

    def test_w_decay_split_preserves_charge_totals(self):
        events = []
        for down, up in [(1, 2), (1, 4), (3, 2), (3, 4), (5, 4), (5, 2)]:
            for charge in [1, -1]:
                for score, weight in [(0.25, 2), (0.75, -.5), (1, 3), (-1, 1), (1.1, 4)]:
                    events.append({"GenQ1PdgId": down, "GenQ2PdgId": -up,
                                   "TriggerLeptonCharge": charge, "finalWeight": weight,
                                   "ak4JetParTPosvsNeg3": score})
        partial, _ = self.partial([self.fixture("decays.root", events)])
        with ROOT.TFile(str(partial)) as output:
            for charge in ["positive", "negative"]:
                stem = "uptype_part_posvsneg_50to80_" + charge
                total = output.Get(stem)
                components = [output.Get(stem + "_Wto" + decay)
                              for decay in ["UD", "CD", "US", "CS", "CB", "UB"]]
                for component in components:
                    self.assertEqual(component.GetEntries(), 5)
                    for bin_index, weight in [(0, 1), (26, 2), (76, -.5), (100, 3), (101, 4)]:
                        self.assertEqual(component.GetBinContent(bin_index), weight)
                        self.assertEqual(component.GetSumw2()[bin_index], weight ** 2)
                for bin_index in range(total.GetNbinsX() + 2):
                    self.assertEqual(total.GetBinContent(bin_index),
                                     sum(h.GetBinContent(bin_index) for h in components))
                    self.assertEqual(total.GetSumw2()[bin_index],
                                     sum(h.GetSumw2()[bin_index] for h in components))
                self.assertEqual(total.GetEntries(), sum(h.GetEntries() for h in components))

    def test_w_decay_uses_saved_flags(self):
        # The split follows the saved flags, not the dataset label or PDG IDs.
        path = self.fixture("flag.root", [{"GenWtoBC": 0, "GenWtoUD": 1}])
        partial, _ = self.partial([path])
        with ROOT.TFile(str(partial)) as output:
            stem = "uptype_part_posvsneg_50to80_positive_Wto"
            self.assertEqual(output.Get(stem + "UD").GetEntries(), 1)
            self.assertEqual(output.Get(stem + "CB").GetEntries(), 0)

    def test_invalid_w_decay_flags_fail(self):
        for event, message in [({"GenWtoBC": 0}, "No W daughter flag"),
                               ({"GenWtoUD": 1}, "Multiple W daughter flags"),
                               ({"GenWtoBC": .5}, "must be 0 or 1")]:
            with self.subTest(event=event):
                path = self.fixture("flags.root", [event])
                _, result = self.partial([path], success=False)
                self.assertIn(message, result.stderr)
        path = self.fixture("missing-flag.root", [{}], missing="GenWtoCS")
        _, result = self.partial([path], success=False)
        self.assertIn("Missing required W daughter flag: GenWtoCS", result.stderr)

    def test_sentinel_remains_underflow(self):
        path = self.fixture("events.root", [{"ak4JetbtagPNetB1": -1}])
        partial, _ = self.partial([path])
        with ROOT.TFile(str(partial)) as output:
            self.assertEqual(output.Get("pnetb_LepB_30to50").GetBinContent(0), 1)

    def test_bad_inputs_fail(self):
        good = self.fixture("good.root", [{}])
        cases = [([good, good], "Duplicate"),
                 ([good, self.directory / "absent.root"], "Cannot read TTree"),
                 ([good, self.fixture("unit.root", [{}], weighted=False)], "schema differs"),
                 ([good, self.fixture("missing.root", [{}], missing="ak4JetPhi2")], "schema differs")]
        for files, message in cases:
            with self.subTest(message=message):
                _, result = self.partial(files, success=False)
                self.assertIn(message, result.stderr)

    def test_invalid_values_fail(self):
        for event, message in [({"finalWeight": float("nan")}, "Non-finite"),
                               ({"GenHadBJetIdx": .1}, "Invalid jet assignment"),
                               ({"GenHadBJetIdx": 10}, "outside saved slots"),
                               ({"GenHadBJetIdx": 1}, "Non-distinct"),
                               ({"TriggerLeptonCharge": 0}, "TriggerLeptonCharge")]:
            with self.subTest(event=event):
                path = self.fixture("invalid.root", [event])
                _, result = self.partial([path], success=False)
                self.assertIn(message, result.stderr)

    def test_merge_matches_direct(self):
        first = self.fixture("first.root", [{"finalWeight": 2},
                                            {"GenQ1PdgId": 1, "GenQ2PdgId": -2}])
        second = self.fixture("second.root", [{"finalWeight": -.5, "GenQ2PdgId": -2},
                                              {"GenQ1PdgId": 3, "GenQ2PdgId": -4,
                                               "TriggerLeptonCharge": -1}])
        self.partial([first], "first.root")
        self.partial([second], "second.root")
        self.run_study([first, second], "direct")
        self.run_study([first, second], "merged", ["--merge", self.directory / "partials"])
        with ROOT.TFile(str(self.directory / "outputs/direct/histograms.root")) as direct, \
             ROOT.TFile(str(self.directory / "outputs/merged/histograms.root")) as merged:
            for key in direct.GetListOfKeys():
                if key.GetName() == "meta_has_final_weight":
                    continue  # number of partial groups, not number of events
                left, right = direct.Get(key.GetName()), merged.Get(key.GetName())
                self.assertEqual(list(left), list(right), key.GetName())
                self.assertEqual(list(left.GetSumw2()), list(right.GetSumw2()), key.GetName())

    def test_unsplit_vcb_partial_requires_regeneration(self):
        source = self.fixture("source.root", [{}])
        partial, _ = self.partial([source])
        with ROOT.TFile(str(partial), "UPDATE") as output:
            for key in list(output.GetListOfKeys()):
                if "_Wto" in key.GetName():
                    output.Delete(key.GetName() + ";*")
        result = self.run_study([source], extra=["--merge", partial.parent], success=False)
        self.assertIn("Partial lacks the UpType W daughter split; rerun the map step", result.stderr)

    def test_mixed_weight_partials_fail(self):
        first = self.fixture("first.root", [{}])
        second = self.fixture("second.root", [{}], weighted=False)
        self.partial([first], "first.root")
        self.partial([second], "second.root")
        result = self.run_study([first], extra=["--merge", self.directory / "partials"], success=False)
        self.assertIn("Cannot mix weighted and unweighted", result.stderr)

    def test_mismatched_axis_fails(self):
        first = self.fixture("first.root", [{}])
        self.partial([first], "first.root")
        second, _ = self.partial([first], "second.root")
        with ROOT.TFile(str(second), "UPDATE") as output:
            hist = output.Get("hadronic_w_mass")
            hist.SetBins(16, 0, 320)
            hist.Write("hadronic_w_mass", ROOT.TObject.kOverwrite)
        result = self.run_study([first], extra=["--merge", self.directory / "partials"], success=False)
        self.assertIn("incompatible axis", result.stderr)


if __name__ == "__main__":
    unittest.main()
