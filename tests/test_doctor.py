from __future__ import annotations

import unittest

from candidate_qtl.core import doctor


class DoctorVersionTests(unittest.TestCase):
    def test_records_python_and_package_versions(self) -> None:
        report = doctor()
        self.assertRegex(report["python_version"], r"^\d+\.\d+\.\d+")
        packages = report["package_versions"]
        # Installed runtime dependencies report a concrete version, not None.
        for dist in ("pandas", "pysam", "biopython", "PyYAML"):
            self.assertIn(dist, packages)
            self.assertIsNotNone(packages[dist])
            self.assertRegex(str(packages[dist]), r"\d")
        # tool_versions only lists tools that are actually installed.
        self.assertIsInstance(report["tool_versions"], dict)


if __name__ == "__main__":
    unittest.main()
