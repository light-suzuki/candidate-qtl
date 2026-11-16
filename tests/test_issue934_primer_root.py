"""Directory inputs retain their type through case-run digest and resume."""
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import yaml

from candidate_qtl.cli import main
from candidate_qtl import workflow_extensions as workflow
from tests import test_issue_255 as fixture


class PrimerRootTests(unittest.TestCase):
    def inputs(self, root, ready=False):
        case, resources, _, _ = fixture.Issue255ContentDigestTests()._evidence_case(root)
        package = root / "primerblast-oss" / "primerblast_oss"
        package.mkdir(parents=True)
        (package / "__init__.py").write_text("# package\n")
        (package / "__main__.py").write_text("# entry 1\n")
        case_data = yaml.safe_load(case.read_text())
        resource_data = yaml.safe_load(resources.read_text())
        case_data["assay"]["primerblast_root_resource"] = "primer_root"
        resource_data["resources"]["primer_root"] = {"path": str(package.parent)}
        if ready:
            for name, text in {"genome": ">chr1\n" + "ACGT" * 100 + "\n",
                               "gff": "##gff-version 3\nchr1\tfixture\tgene\t1\t40\t.\t+\t.\tID=GeneA\n",
                               "genomic": ">GeneA\nACGTACGT\n>GeneB\nACGTACGT\n", "db": "fixture DB"}.items():
                path = root / (name + ".txt")
                path.write_text(text)
                resource_data["resources"][name] = {"path": str(path)}
            case_data["target"] = {"fasta_resource": "genome", "gff3_resource": "gff", "nucleotide_db_resource": "db"}
            case_data["inputs"]["candidate_genomic_resource"] = "genomic"
        case.write_text(yaml.safe_dump(case_data))
        resources.write_text(yaml.safe_dump(resource_data))
        return case, resources, package

    def run_case(self, root, case, resources, *flags):
        self.assertEqual(0, main(["case-run", "--case", str(case), "--resources", str(resources),
                                 "--out-dir", str(root / "out"), *flags]))
        return json.loads((root / "out/case_state.json").read_text())

    def test_optional_root_finishes_without_assay_and_resume_reacts_to_code(self):
        for strict in (False, True):
            with self.subTest(strict=strict), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                case, resources, package = self.inputs(root)
                flags = ["--strict-fingerprints"] if strict else []
                first = self.run_case(root, case, resources, *flags)
                self.assertEqual("REPORT_READY", first["state"])
                self.assertFalse((root / "out/primer/plan").exists())
                snapshot = (root / "out/rank/ranked_candidates.tsv").read_bytes()
                self.assertEqual(first["steps"]["primer"], self.run_case(root, case, resources, "--resume", *flags)["steps"]["primer"])
                entry = package / "__main__.py"
                stat = entry.stat()
                entry.write_text("# entry 2\n")
                os.utime(entry, ns=(stat.st_atime_ns, stat.st_mtime_ns))
                changed = self.run_case(root, case, resources, "--resume", *flags)
                self.assertNotEqual(first["steps"]["primer"]["input_hash"], changed["steps"]["primer"]["input_hash"])
                self.assertEqual(snapshot, (root / "out/rank/ranked_candidates.tsv").read_bytes())
                fp = changed["steps"]["primer"]["primerblast_root_fingerprint"]
                self.assertEqual("directory", fp["resource_type"])
                self.assertIn("primerblast_oss/__main__.py", fp["files"])
                (package.parent / "huge-db.bin").write_bytes(b"unrelated database")
                self.assertEqual(fp, workflow._primerblast_root_fingerprint(package.parent))

    def test_ready_inputs_reach_real_plan_and_root_errors_remain_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            case, resources, package = self.inputs(root, ready=True)
            with patch.object(workflow, "primer_plan", wraps=workflow.primer_plan) as plan:
                self.run_case(root, case, resources)
                self.assertEqual(1, plan.call_count)
            self.assertTrue((root / "out/primer/plan/run_direct_sequencing.sh").is_file())
            (package / "__main__.py").unlink()
            with self.assertRaisesRegex(ValueError, "lacks primerblast_oss"):
                self.run_case(root, case, resources, "--resume")
            with self.assertRaisesRegex(ValueError, "must be a directory"):
                workflow._primerblast_root_fingerprint(package / "__init__.py")
            with self.assertRaisesRegex(ValueError, "must be a directory"):
                workflow._primerblast_root_fingerprint(root / "missing")
