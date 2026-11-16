"""Tests for the normalized provenance DAG ledger (issue #231)."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.cli import build_parser, main
from candidate_qtl.provenance_ledger import provenance_ledger, render_ledger_report, shared_derivation_paths

# Two route outputs (routeA, routeB) both derive from the same PMID through one
# shared database_annotation node: the exact structure issue #231 asks to
# surface as a shared derivation path.
SHARED_ENTRIES = (
    "node_id\tnode_kind\tsource_id\tlabel\tparent_node\tsource_ref\n"
    "exp1\tprimary_experiment\tPMID:1234\toriginal experiment\t\tpaper DOI\n"
    "ann1\tdatabase_annotation\tPMID:1234\tdatabase annotation entry\texp1\tdb version 1\n"
    "routeA\troute_output\trouteA\troute A output\tann1\t\n"
    "routeB\troute_output\trouteB\troute B output\tann1\t\n"
)

# Two fully independent derivations: no shared source anywhere.
INDEPENDENT_ENTRIES = (
    "node_id\tnode_kind\tsource_id\tlabel\tparent_node\n"
    "exp1\tprimary_experiment\tPMID:111\tpaper one\t\n"
    "routeA\troute_output\trouteA\troute A output\texp1\n"
    "exp2\tprimary_experiment\tPMID:222\tpaper two\t\n"
    "routeB\troute_output\trouteB\troute B output\texp2\n"
)

# The full chain named in issue #231: primary_experiment -> database_annotation
# -> review -> route_output. Both route outputs are only reachable through the
# same review article, so a convergence between them is shared-source, not
# independent replication, even though neither route cites the review directly.
REVIEW_CHAIN_ENTRIES = (
    "node_id\tnode_kind\tsource_id\tlabel\tparent_node\tsource_ref\n"
    "exp1\tprimary_experiment\tPMID:5000\toriginal experiment\t\tpaper DOI\n"
    "ann1\tdatabase_annotation\tPMID:5000\tdatabase annotation entry\texp1\tdb version 1\n"
    "rev1\treview\tPMID:5000\treview article citing the annotation\tann1\t\n"
    "routeA\troute_output\trouteA\troute A output\trev1\t\n"
    "routeB\troute_output\trouteB\troute B output\trev1\t\n"
)

# Same shape, but the mediating node is a literature citation rather than a
# review article: primary_experiment -> database_annotation -> literature ->
# route_output.
LITERATURE_CHAIN_ENTRIES = (
    "node_id\tnode_kind\tsource_id\tlabel\tparent_node\tsource_ref\n"
    "exp2\tprimary_experiment\tPMID:6000\toriginal experiment\t\tpaper DOI\n"
    "ann2\tdatabase_annotation\tPMID:6000\tdatabase annotation entry\texp2\tdb version 1\n"
    "lit2\tliterature\tPMID:6000\tliterature citation of the annotation\tann2\t\n"
    "routeC\troute_output\trouteC\troute C output\tlit2\t\n"
    "routeD\troute_output\trouteD\troute D output\tlit2\t\n"
)


class ProvenanceLedgerTests(unittest.TestCase):
    def _build(self, root: Path, entries_text: str, case_id: str = "") -> tuple[Path, Path, dict]:
        entries = root / "entries.tsv"
        entries.write_text(entries_text, encoding="utf-8")
        ledger_dir, out_dir = root / "ledger", root / "out"
        result = provenance_ledger(str(entries), str(ledger_dir), str(out_dir), case_id=case_id)
        return ledger_dir, out_dir, result

    def test_shared_derivation_pair_is_detected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            ledger_dir, _, _ = self._build(Path(tmp), SHARED_ENTRIES)
            shared = shared_derivation_paths(str(ledger_dir))
            self.assertEqual(1, len(shared))
            entry = shared[0]
            self.assertEqual("PMID:1234", entry["source_id"])
            self.assertEqual(["routeA", "routeB"], entry["route_outputs"])
            self.assertEqual([["routeA", "routeB"]], entry["pairs"])

    def test_independent_sources_have_no_shared_paths(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            ledger_dir, _, _ = self._build(Path(tmp), INDEPENDENT_ENTRIES)
            self.assertEqual([], shared_derivation_paths(str(ledger_dir)))

    def test_review_mediated_chain_resolves_to_shared_root(self) -> None:
        # primary_experiment -> database_annotation -> review -> route_output
        # (both routes). shared_derivation_paths must still walk through the
        # review hop and resolve both outputs to the same primary-experiment
        # root, so a convergence reached only via a review article is flagged
        # as shared-source rather than independent.
        with tempfile.TemporaryDirectory() as tmp:
            ledger_dir, _, _ = self._build(Path(tmp), REVIEW_CHAIN_ENTRIES)
            edges = pd.read_csv(ledger_dir / "edges.tsv", sep="\t", dtype=str).fillna("").set_index("child_node")
            self.assertEqual("reviewed_in", edges.loc["rev1", "edge_kind"])
            shared = shared_derivation_paths(str(ledger_dir))
            self.assertEqual(1, len(shared))
            entry = shared[0]
            self.assertEqual("PMID:5000", entry["source_id"])
            self.assertEqual(["routeA", "routeB"], entry["route_outputs"])
            self.assertEqual([["routeA", "routeB"]], entry["pairs"])

    def test_literature_mediated_chain_resolves_to_shared_root(self) -> None:
        # primary_experiment -> database_annotation -> literature -> route_output
        # (both routes), mirroring the review case above for the literature
        # node kind.
        with tempfile.TemporaryDirectory() as tmp:
            ledger_dir, _, _ = self._build(Path(tmp), LITERATURE_CHAIN_ENTRIES)
            edges = pd.read_csv(ledger_dir / "edges.tsv", sep="\t", dtype=str).fillna("").set_index("child_node")
            self.assertEqual("reviewed_in", edges.loc["lit2", "edge_kind"])
            shared = shared_derivation_paths(str(ledger_dir))
            self.assertEqual(1, len(shared))
            entry = shared[0]
            self.assertEqual("PMID:6000", entry["source_id"])
            self.assertEqual(["routeC", "routeD"], entry["route_outputs"])
            self.assertEqual([["routeC", "routeD"]], entry["pairs"])

    def test_parentless_source_becomes_primary_experiment_root(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            ledger_dir, _, _ = self._build(Path(tmp), SHARED_ENTRIES)
            nodes = pd.read_csv(ledger_dir / "nodes.tsv", sep="\t", dtype=str).fillna("").set_index("node_id")
            self.assertEqual("primary_experiment", nodes.loc["exp1", "node_kind"])
            self.assertEqual("PMID:1234", nodes.loc["exp1", "source_id"])
            self.assertEqual("original experiment", nodes.loc["exp1", "label"])

    def test_auto_generated_node_ids_for_children(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            entries = root / "entries.tsv"
            entries.write_text(
                "node_kind\tsource_id\tlabel\tparent_node\n"
                "database_annotation\tDB:999\tannotation\tPMID:1\n"
                "primary_experiment\tPMID:1\tpaper\t\n",
                encoding="utf-8",
            )
            ledger_dir, _, _ = self._build(root, entries.read_text(encoding="utf-8"))
            nodes = pd.read_csv(ledger_dir / "nodes.tsv", sep="\t", dtype=str).fillna("").set_index("node_id")
            edges = pd.read_csv(ledger_dir / "edges.tsv", sep="\t", dtype=str).fillna("")
            self.assertIn("database_annotation:DB:999", nodes.index)
            self.assertEqual("primary_experiment", nodes.loc["PMID:1", "node_kind"])
            self.assertEqual("annotates", edges.set_index("child_node").loc["database_annotation:DB:999", "edge_kind"])

    def test_duplicate_node_id_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            entries = root / "entries.tsv"
            entries.write_text(
                "node_id\tnode_kind\tsource_id\tlabel\tparent_node\n"
                "dup\tprimary_experiment\ts1\ta\t\n"
                "dup\troute_output\ts2\tb\t\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "duplicate node_id"):
                self._build(root, entries.read_text(encoding="utf-8"))

    def test_edge_to_missing_parent_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            entries = root / "entries.tsv"
            entries.write_text(
                "node_kind\tsource_id\tlabel\tparent_node\n"
                "route_output\trouteA\toutput\tNO_SUCH_NODE\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "unknown parent"):
                self._build(root, entries.read_text(encoding="utf-8"))

    def test_child_equals_parent_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            entries = root / "entries.tsv"
            entries.write_text(
                "node_id\tnode_kind\tsource_id\tlabel\tparent_node\n"
                "loop\troute_output\ts1\toutput\tloop\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "self-loop"):
                self._build(root, entries.read_text(encoding="utf-8"))

    def test_all_four_output_files_are_written(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            entries = root / "entries.tsv"
            entries.write_text(SHARED_ENTRIES, encoding="utf-8")
            ledger_dir, out_dir = root / "ledger", root / "out"
            self.assertEqual(0, main(["provenance-ledger", "--entries", str(entries),
                                      "--ledger-dir", str(ledger_dir), "--out-dir", str(out_dir),
                                      "--case-id", "case-001"]))
            for path in (ledger_dir / "nodes.tsv", ledger_dir / "edges.tsv",
                         out_dir / "provenance_graph.tsv", out_dir / "provenance_ledger_manifest.json"):
                self.assertTrue(path.exists(), path)
            nodes = pd.read_csv(ledger_dir / "nodes.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(4, len(nodes))
            edges = pd.read_csv(ledger_dir / "edges.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(3, len(edges))
            self.assertEqual(["annotates", "derived_from", "derived_from"], sorted(edges["edge_kind"].tolist()))
            graph = pd.read_csv(out_dir / "provenance_graph.tsv", sep="\t", dtype=str).fillna("")
            self.assertEqual(4, len(graph))
            self.assertIn("parent_node", graph.columns)
            manifest = json.loads((out_dir / "provenance_ledger_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual("case-001", manifest["case_id"])
            self.assertEqual(4, manifest["node_count"])
            self.assertEqual(1, manifest["shared_derivation_source_count"])
            self.assertEqual(1, manifest["shared_derivation_pair_count"])
            self.assertTrue(any("shared by multiple routes is counted once" in line for line in manifest["scientific_limits"]))

    def test_report_names_shared_source_and_limit(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            ledger_dir, _, _ = self._build(Path(tmp), SHARED_ENTRIES)
            report = render_ledger_report(str(ledger_dir))
            self.assertIn("PMID:1234", report)
            self.assertIn("counted once", report)


class ProvenanceLedgerCliTests(unittest.TestCase):
    def test_provenance_ledger_command_is_registered(self) -> None:
        parser = build_parser()
        subparser_action = next(action for action in parser._subparsers._group_actions if hasattr(action, "choices"))
        self.assertIn("provenance-ledger", subparser_action.choices)
        subparser = subparser_action.choices["provenance-ledger"]
        optionals = [action.dest for action in subparser._actions]
        for flag in ("entries", "ledger_dir", "out_dir", "case_id", "json"):
            self.assertIn(flag, optionals, flag)

    def test_cli_json_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            entries = root / "entries.tsv"
            entries.write_text(SHARED_ENTRIES, encoding="utf-8")
            ledger_dir, out_dir = root / "ledger", root / "out"
            result = main(["provenance-ledger", "--entries", str(entries),
                           "--ledger-dir", str(ledger_dir), "--out-dir", str(out_dir), "--json"])
            self.assertEqual(0, result)


if __name__ == "__main__":
    unittest.main()
