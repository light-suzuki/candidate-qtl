"""Unit tests for the deterministic blast-panel freeze layer (issue #374 / #355).

Covers the row-level verdicts (V-ACC/V-SEQ/V-AA/V-LEN/V-GENERIC/V-SHA/V-PROV),
the fail-closed column rules (V-COL/V-FORBIDDEN/V-DUP/V-RANK), the freeze /
verify roundtrip including drift detection (NC5), determinism, manifest
content, the pre-interval firewall and min_panel_size handling.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from candidate_qtl.blast_panel import (
    DEFAULT_THRESHOLDS,
    PANEL_COLUMNS,
    TRUTH_ABSENCE_ASSERTION,
    _forbidden_columns,
    _generic_motif,
    _has_local_tandem_array,
    _row_sha256,
    freeze_blast_panel,
    validate_panel_row,
    verify_blast_panel_freeze,
)

FIXTURE = Path(__file__).resolve().parent / "fixtures" / "route_a"

S = "CDLYQWLICSGRWNGGYHCLKFRVGRYVDSYIMRQYFNMFLEVAWCGNRTDNRMNHMNNQDKTTEQPNHVNKQKNLGSRNSEYNGGSQYKGNLRRCKLTIYEEPFSWSDAPQCWSCIFET"
PA = "MYTCNTWCWQTEEFKCHFAIGTHRPYTYCTNKWFDNQDNCCEEVEEAAPRMLMKADWYYRIDDDAPRQVSDLGPMAKINYNRGSIDRINHLCVWDLHQHYYWWCPFYKAFREYFAKEKQR"
SHORT20 = "MKTLFVAEGSINPDRQHCYW"
CYCLIC20 = "CDEFGHIKLMNPQRSTVWY"
PROMPT_FP = "0123456789abcdef" * 4

# Public UniProt sequences (retrieved via the issue #374 r2 UniProt cache;
# byte-identical to the frozen bait panels of runs-issue374-r2). The three
# LysM receptor kinases were false positives of the old global dipeptide
# count; OsSGR/TT8/LDOX also fired it. All six must pass the local tandem rule.
NFR1_Q70KR8 = "MKLKTGLLLFFILLLGHVCFHVESNCLKGCDLALASYYILPGVFILQNITTFMQSEIVSSNDAITSYNKDKILNDINIQSFQRLNIPFPCDCIGGEFLGHVFEYSASKGDTYETIANLYYANLTTVDLLKRFNSYDPKNIPVNAKVNVTVNCSCGNSQVSKDYGLFITYPIRPGDTLQDIANQSSLDAGLIQSFNPSVNFSKDSGIAFIPGRYKNGVYVPLYHRTAGLASGAAVGISIAGTFVLLLLAFCMYVRYQKKEEEKAKLPTDISMALSTQDASSSAEYETSGSSGPGTASATGLTSIMVAKSMEFSYQELAKATNNFSLDNKIGQGGFGAVYYAELRGKKTAIKKMDVQASTEFLCELKVLTHVHHLNLVRLIGYCVEGSLFLVYEHIDNGNLGQYLHGSGKEPLPWSSRVQIALDAARGLEYIHEHTVPVYIHRDVKSANILIDKNLRGKVADFGLTKLIEVGNSTLQTRLVGTFGYMPPEYAQYGDISPKIDVYAFGVVLFELISAKNAVLKTGELVAESKGLVALFEEALNKSDPCDALRKLVDPRLGENYPIDSVLKIAQLGRACTRDNPLLRPSMRSLVVALMTLSSLTEDCDDESSYESQTLINLLSVR"
LYK3_Q6UD73 = "MNLKNGLLLFILFLDCVFFKVESKCVKGCDVALASYYIIPSIQLRNISNFMQSKIVLTNSFDVIMSYNRDVVFDKSGLISYTRINVPFPCECIGGEFLGHVFEYTTKEGDDYDLIANTYYASLTTVELLKKFNSYDPNHIPVKAKINVTVICSCGNSQISKDYGLFVTYPLRSDDTLAKIATKAGLDEGLIQNFNQDANFSIGSGIVFIPGRDQNGHFFPLYSRTGIAKGSAVGIAMAGIFGLLLFVIYIYAKYFQKKEEEKTKLPQTSRAFSTQDASGSAEYETSGSSGHATGSAAGLTGIMVAKSTEFTYQELAKATNNFSLDNKIGQGGFGAVYYAELRGEKTAIKKMDVQASSEFLCELKVLTHVHHLNLVRLIGYCVEGSLFLVYEHIDNGNLGQYLHGIGTEPLPWSSRVQIALDSARGLEYIHEHTVPVYIHRDVKSANILIDKNLRGKVADFGLTKLIEVGNSTLHTRLVGTFGYMPPEYAQYGDVSPKIDVYAFGVVLYELITAKNAVLKTGESVAESKGLVQLFEEALHRMDPLEGLRKLVDPRLKENYPIDSVLKMAQLGRACTRDNPLLRPSMRSIVVALMTLSSPTEDCDDDSSYENQSLINLLSTR"
NFP_Q0GXS4 = "MSAFFLPSSSHALFLVLMLFFLTNISAQPLYISETNFTCPVDSPPSCETYVAYRAQSPNFLSLSNISDIFNLSPLRIAKASNIEAEDKKLIPDQLLLVPVTCGCTKNHSFANITYSIKQGDNFFILSITSYQNLTNYLEFKNFNPNLSPTLLPLDTKVSVPLFCKCPSKNQLNKGIKYLITYVWQDNDNVTLVSSKFGASQVEMLAENNHNFTASTNRSVLIPVTSLPKLDQPSSNGRKSSSQNLALIIGISLGSAFFILVLTLSLVYVYCLKMKRLNRSTSSSETADKLLSGVSGYVSKPTMYEIDAIMEGTTNLSDNCKIGESVYKANIDGRVLAVKKIKKDASEELKILQKVNHGNLVKLMGVSSDNDGNCFLVYEYAENGSLEEWLFSESSKTSNSVVSLTWSQRITIAMDVAIGLQYMHEHTYPRIIHRDITTSNILLGSNFKAKIANFGMARTSTNSMMPKIDVFAFGVVLIELLTGKKAMTTKENGEVVILWKDFWKIFDLEGNREERLRKWMDPKLESFYPIDNALSLASLAVNCTADKSLSRPTIAEIVLCLSLLNQPSSEPMLERSLTSGLDAEATHVVTSVIAR"
OSSGR_Q652K1 = "MAAATSTMSLLPPITQQQRWHAADSLVVLASRCHNSRRRRRCRYVVPRARLFGPAIFEASKLKVLFLGVDEEKHQHPGKLPRTYTLTHSDVTARLTLAVSHTINRAQLQGWYNKLQRDEVVAEWKKVQGHMSLHVHCHISGGHVLLDLIAGLRYYIFRKELPVVLKAFVHGDGNLFSRHPELEEATVWVYFHSNLPRFNRVECWGPLRDAGAPPEEDDAVAAAAAEEAAAEQMPAAGEWPRRCPGQCDCCFPPYSLIPWPHQHDVAAADGQPQQ"
TT8_Q9FT81 = "MDESSIIPAEKVAGAEKKELQGLLKTAVQSVDWTYSVFWQFCPQQRVLVWGNGYYNGAIKTRKTTQPAEVTAEEAALERSQQLRELYETLLAGESTSEARACTALSPEDLTETEWFYLMCVSFSFPPPSGMPGKAYARRKHVWLSGANEVDSKTFSRAILAKSAKIQTVVCIPMLDGVVELGTTKKVREDVEFVELTKSFFYDHCKTNPKPALSEHSTYEVHEEAEDEEEVEEEMTMSEEMRLGSPDDEDVSNQNLHSDLHIESTHTLDTHMDMMNLMEEGGNYSQTVTTLLMSHPTSLLSDSVSTSSYIQSSFATWRVENGKEHQQVKTAPSSQWVLKQMIFRVPFLHDNTKDKRLPREDLSHVVAERRRREKLNEKFITLRSMVPFVTKMDKVSILGDTIAYVNHLRKRVHELENTHHEQQHKRTRTCKRKTSEEVEVSIIENDVLLEMRCEYRDGLLLDILQVLHELGIETTAVHTSVNDHDFEAEIRAKVRGKKASIAEVKRAIHQVIIHDTNL"
LDOX_Q96323 = "MVAVERVESLAKSGIISIPKEYIRPKEELESINDVFLEEKKEDGPQVPTIDLKNIESDDEKIRENCIEELKKASLDWGVMHLINHGIPADLMERVKKAGEEFFSLSVEEKEKYANDQATGKIQGYGSKLANNASGQLEWEDYFFHLAYPEEKRDLSIWPKTPSDYIEATSEYAKCLRLLATKVFKALSVGLGLEPDRLEKEVGGLEELLLQMKINYYPKCPQPELALGVEAHTDVSALTFILHNMVPGLQLFYEGKWVTAKCVPDSIVMHIGDTLEILSNGKYKSILHRGLVNKEKVRISWAVFCEPPKDKIVLKPLPEMVSVESPAKFPPRTFAQHIEHKLFGKEQEELVSEKND"


def _row(**overrides: str) -> dict[str, str]:
    row = {
        "panel_id": "row1",
        "mechanism": "mechanism_x",
        "family": "family_x",
        "synonyms": "syn1",
        "reference_gene_name": "GENE1",
        "protein_accession": "Q9TES1",
        "taxon": "Testis syntheticus",
        "source": "local_record_cache",
        "source_url": "file://fixture",
        "retrieved_at": "2026-08-12T00:00:00+00:00",
        "verification_status": "verified",
        "verification_accession": "Q9TES1",
        "verification_sequence_sha256": _row_sha256(S),
        "verification_source_url": "file://fixture/verified-record.json",
        "verification_retrieved_at": "2026-08-12T00:00:00+00:00",
        "trait_direction": "increase",
        "query_sequence": S,
        "query_sequence_sha256": _row_sha256(S),
        "provider": "test_fixture",
        "model": "none",
        "prompt_fingerprint": PROMPT_FP,
        "review_status": "verified",
        "review_reason": "",
        "rank_contribution": "none",
    }
    row.update(overrides)
    if "protein_accession" in overrides and "verification_accession" not in overrides:
        row["verification_accession"] = row["protein_accession"]
    if "query_sequence" in overrides and "verification_sequence_sha256" not in overrides:
        row["verification_sequence_sha256"] = _row_sha256(row["query_sequence"])
    return row


def _panel_tsv(root: Path, rows: list[dict[str, str]]) -> Path:
    path = root / "panel.tsv"
    path.write_text(
        "\t".join(PANEL_COLUMNS) + "\n"
        + "".join("\t".join(str(row.get(column, "")) for column in PANEL_COLUMNS) + "\n" for row in rows),
        encoding="utf-8",
    )
    return path


def _freeze(root: Path, rows: list[dict[str, str]], **kwargs: object) -> dict[str, object]:
    panel = _panel_tsv(root, rows)
    return freeze_blast_panel(str(panel), str(root / "out"), **kwargs)


class RowValidationTests(unittest.TestCase):
    def test_verified_row_passes(self) -> None:
        verdict = validate_panel_row(_row())
        self.assertEqual("verified", verdict["review_status"])
        self.assertEqual("", verdict["review_reason"])
        self.assertEqual(len(S), verdict["query_length"])
        self.assertEqual(S, verdict["query_sequence"])
        self.assertFalse(verdict["is_generic_motif"])

    def test_sequence_is_upper_cased_and_sha_recomputed(self) -> None:
        verdict = validate_panel_row(_row(query_sequence=S.lower(), query_sequence_sha256="wrong"))
        self.assertEqual(S, verdict["query_sequence"])
        self.assertEqual(_row_sha256(S), verdict["query_sequence_sha256"])

    def test_bad_accession_format_is_first_verdict(self) -> None:
        # V-ACC fires before V-SEQ/V-AA/V-LEN: even a broken sequence cannot
        # mask a broken accession.
        for accession in ("P0FAKE", "abc123", "", "UNIPROT1", "Q0TEST", "Q0TEST1"):
            verdict = validate_panel_row(_row(protein_accession=accession, query_sequence=""))
            self.assertEqual("review_required", verdict["review_status"], accession)
            self.assertEqual("bad_accession_format", verdict["review_reason"], accession)

    def test_valid_accession_formats_pass(self) -> None:
        for accession in ("Q9TES1", "P93003", "Q9LZS3", "A0ACM8R380", "F8T7T1", "NP_123456", "XP_001234.5", "YP_987654", "WP_000000001", "ENSP00000269305"):
            verdict = validate_panel_row(_row(protein_accession=accession, verification_accession=accession))
            self.assertNotEqual("bad_accession_format", verdict["review_reason"], accession)

    def test_empty_sequence(self) -> None:
        verdict = validate_panel_row(_row(query_sequence="", query_sequence_sha256=""))
        self.assertEqual("review_required", verdict["review_status"])
        self.assertEqual("empty_sequence", verdict["review_reason"])

    def test_non_aa_and_high_unknown_fraction(self) -> None:
        high_unknown = "X" * 31 + CYCLIC20 + "ACDEFGHIK"  # 60 aa, X fraction > 0.5
        verdict = validate_panel_row(_row(query_sequence=high_unknown, query_sequence_sha256=_row_sha256(high_unknown)))
        self.assertEqual("invalid_aa_or_high_unknown", verdict["review_reason"])

        has_b_residue = "B" + "A" * 59
        verdict = validate_panel_row(_row(query_sequence=has_b_residue, query_sequence_sha256=_row_sha256(has_b_residue)))
        self.assertEqual("invalid_aa_or_high_unknown", verdict["review_reason"])

        low_unknown = "X" * 14 + CYCLIC20 * 2 + "ACDEFGHIKLMN"  # 66 aa, X fraction 0.21
        verdict = validate_panel_row(_row(query_sequence=low_unknown, query_sequence_sha256=_row_sha256(low_unknown)))
        self.assertNotEqual("invalid_aa_or_high_unknown", verdict["review_reason"])

    def test_short_query(self) -> None:
        verdict = validate_panel_row(_row(query_sequence=SHORT20, query_sequence_sha256=_row_sha256(SHORT20)))
        self.assertEqual("short_query", verdict["review_reason"])
        self.assertEqual(20, verdict["query_length"])

    def test_generic_motif_single_aa_fraction(self) -> None:
        seq = "A" * 20 + CYCLIC20 * 2  # 60 aa, A fraction 0.5, run 20
        triggered, detail = _generic_motif(seq)
        self.assertTrue(triggered)
        self.assertEqual("single_amino_acid_fraction", detail)
        verdict = validate_panel_row(_row(query_sequence=seq, query_sequence_sha256=_row_sha256(seq)))
        self.assertEqual("generic_motif", verdict["review_reason"])
        self.assertTrue(verdict["is_generic_motif"])

    def test_generic_motif_homopolymer_run(self) -> None:
        seq = "A" * 12 + (CYCLIC20 * 5)[:108]  # 120 aa, run 12, A fraction 0.1
        triggered, detail = _generic_motif(seq)
        self.assertTrue(triggered)
        self.assertEqual("homopolymer_run", detail)
        verdict = validate_panel_row(_row(query_sequence=seq, query_sequence_sha256=_row_sha256(seq)))
        self.assertEqual("generic_motif", verdict["review_reason"])

    def test_generic_motif_dipeptide_repeat(self) -> None:
        seq = "ACDEG" * 20  # 100 aa: each AA fraction 0.2, run 1, perfect period-5 tandem
        self.assertTrue(_has_local_tandem_array(seq))
        triggered, detail = _generic_motif(seq)
        self.assertTrue(triggered)
        self.assertEqual("tandem_repeat", detail)
        verdict = validate_panel_row(_row(query_sequence=seq, query_sequence_sha256=_row_sha256(seq)))
        self.assertEqual("generic_motif", verdict["review_reason"])

    def test_generic_motif_tandem_repeat_boundary_values(self) -> None:
        # min_len=12 (SEG local window): a perfect period-2 tandem span of
        # exactly 12 residues fires; 10 residues stays below the boundary.
        self.assertTrue(_has_local_tandem_array("KL" * 6))
        self.assertFalse(_has_local_tandem_array("KL" * 5))
        self.assertTrue(_has_local_tandem_array("AC" * 6))
        self.assertFalse(_has_local_tandem_array("KL" * 5 + "K"))

    def test_generic_motif_short_homopolymer_does_not_trip_tandem_rule(self) -> None:
        # A 5-residue "AAAAA" stretch inside a longer sequence has no >= 12 aa
        # tandem span: neither the helper nor the generic-motif gate may fire.
        seq = "M" + "AAAAA" + CYCLIC20 * 3
        self.assertFalse(_has_local_tandem_array(seq))
        triggered, detail = _generic_motif(seq)
        self.assertFalse(triggered)
        self.assertEqual("", detail)

    def test_natural_lysm_receptor_kinases_do_not_trip_tandem_repeat(self) -> None:
        # NFR1 / LYK3 / NFP were false positives of the old global dipeptide
        # count (LL/GL/SS scattered pairs). The local tandem rule must not fire.
        for accession, seq in (("Q70KR8", NFR1_Q70KR8), ("Q6UD73", LYK3_Q6UD73), ("Q0GXS4", NFP_Q0GXS4)):
            with self.subTest(accession=accession):
                self.assertFalse(_has_local_tandem_array(seq))
                triggered, detail = _generic_motif(seq)
                self.assertFalse(triggered)
                self.assertEqual("", detail)

    def test_other_natural_bait_proteins_do_not_trip_tandem_repeat(self) -> None:
        # OsSGR / TT8 / LDOX also fired the old dipeptide count; none of them
        # holds a >= 12 aa perfect tandem span, so the new rule passes them.
        for accession, seq in (("Q652K1", OSSGR_Q652K1), ("Q9FT81", TT8_Q9FT81), ("Q96323", LDOX_Q96323)):
            with self.subTest(accession=accession):
                self.assertFalse(_has_local_tandem_array(seq))
                triggered, detail = _generic_motif(seq)
                self.assertFalse(triggered)
                self.assertEqual("", detail)

    def test_sequence_fingerprint_mismatch(self) -> None:
        verdict = validate_panel_row(_row(query_sequence_sha256="0" * 64))
        self.assertEqual("sequence_fingerprint_mismatch", verdict["review_reason"])
        verdict = validate_panel_row(_row(query_sequence_sha256=""))
        self.assertEqual("sequence_fingerprint_mismatch", verdict["review_reason"])

    def test_missing_provenance(self) -> None:
        verdict = validate_panel_row(_row(source=""))
        self.assertEqual("missing_provenance", verdict["review_reason"])

    def test_missing_or_mismatched_verification_evidence_cannot_promote_a_query(self) -> None:
        cases = (
            ({"verification_status": ""}, "missing_verification_evidence"),
            ({"verification_status": "pending"}, "verification_not_verified"),
            ({"verification_accession": "P93003"}, "verification_accession_mismatch"),
            ({"verification_sequence_sha256": "0" * 64}, "verification_sequence_fingerprint_mismatch"),
        )
        for overrides, reason in cases:
            verdict = validate_panel_row(_row(**overrides))
            self.assertEqual("review_required", verdict["review_status"])
            self.assertEqual(reason, verdict["review_reason"])
        # provider alone does not count as provenance; all four fields empty is.
        verdict = validate_panel_row(_row(source_url="", retrieved_at="", model="", prompt_fingerprint=""))
        self.assertEqual("missing_provenance", verdict["review_reason"])

    def test_explicit_rejected_input_stays_rejected(self) -> None:
        verdict = validate_panel_row(_row(review_status="rejected", review_reason="rejected_by_input"))
        self.assertEqual("rejected", verdict["review_status"])

    def test_input_review_required_is_sticky_even_when_checks_pass(self) -> None:
        # F1: an input review_required row keeps its status and reason even
        # when every mechanical check passes; it never upgrades to verified.
        verdict = validate_panel_row(_row(review_status="review_required", review_reason="manual_hold"))
        self.assertEqual("review_required", verdict["review_status"])
        self.assertEqual("manual_hold", verdict["review_reason"])
        self.assertEqual(len(S), verdict["query_length"])

    def test_input_not_evaluated_is_sticky_even_when_checks_pass(self) -> None:
        # F1: input not_evaluated stays not_evaluated regardless of checks.
        verdict = validate_panel_row(_row(review_status="not_evaluated", review_reason=""))
        self.assertEqual("not_evaluated", verdict["review_status"])

    def test_input_review_required_does_not_mask_mechanical_normalization(self) -> None:
        # F1: sticky status, but the sequence is still upper-cased and the
        # SHA-256 fingerprint is still recomputed.
        verdict = validate_panel_row(_row(review_status="review_required", query_sequence=S.lower(), query_sequence_sha256="wrong"))
        self.assertEqual("review_required", verdict["review_status"])
        self.assertEqual(S, verdict["query_sequence"])
        self.assertEqual(_row_sha256(S), verdict["query_sequence_sha256"])

    def test_input_review_status_outside_vocab_is_fail_safe_review_required(self) -> None:
        # F1: any string outside REVIEW_STATUS_VOCAB is a fail-safe review_required.
        verdict = validate_panel_row(_row(review_status="arbitrary_string"))
        self.assertEqual("review_required", verdict["review_status"])
        self.assertEqual("invalid_input_review_status", verdict["review_reason"])

    def test_empty_review_status_with_failing_check_is_review_required(self) -> None:
        # F1: empty input status lets the mechanical checks decide.
        verdict = validate_panel_row(_row(review_status="", query_sequence="", query_sequence_sha256=""))
        self.assertEqual("review_required", verdict["review_status"])
        self.assertEqual("empty_sequence", verdict["review_reason"])

    def test_verified_input_with_failing_check_is_review_required(self) -> None:
        # F1: a verified input that fails a mechanical check is demoted.
        verdict = validate_panel_row(_row(review_status="verified", query_sequence="", query_sequence_sha256=""))
        self.assertEqual("review_required", verdict["review_status"])
        self.assertEqual("empty_sequence", verdict["review_reason"])

    def test_fixture_bait4_review_required_is_not_promoted_to_verified(self) -> None:
        # F1: the committed fixture bait4 (review_required) must stay
        # review_required through the freeze, never verified.
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            result = freeze_blast_panel(str(FIXTURE / "panel_ok.tsv"), str(root / "out"), case_id="bait4_sticky")
            panel = pd.read_csv(result["panel"], sep="\t", dtype=str, keep_default_na=False).fillna("")
            bait4 = panel.loc[panel["panel_id"] == "bait4", ["review_status", "review_reason"]].iloc[0]
            self.assertEqual("review_required", bait4["review_status"])
            self.assertEqual("short_query", bait4["review_reason"])

    def test_trait_direction_normalized(self) -> None:
        self.assertEqual("unknown", validate_panel_row(_row(trait_direction="sideways"))["trait_direction"])
        self.assertEqual("increase", validate_panel_row(_row(trait_direction="INCREASE"))["trait_direction"])
        self.assertEqual("unknown", validate_panel_row(_row(trait_direction=""))["trait_direction"])
        self.assertEqual("loss", validate_panel_row(_row(trait_direction="loss"))["trait_direction"])


class FreezeInputRulesTests(unittest.TestCase):
    def test_missing_column_raises_value_error(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            rows = [_row()]
            del rows[0]["taxon"]
            path = root / "panel.tsv"
            header = [column for column in PANEL_COLUMNS if column != "taxon"]
            path.write_text(
                "\t".join(header) + "\n"
                + "\t".join(str(rows[0].get(column, "")) for column in header) + "\n",
                encoding="utf-8",
            )
            with self.assertRaises(ValueError) as ctx:
                freeze_blast_panel(str(path), str(root / "out"))
            self.assertIn("taxon", str(ctx.exception))

    def test_forbidden_column_fixture_raises_value_error(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            with self.assertRaises(ValueError) as ctx:
                freeze_blast_panel(str(FIXTURE / "panel_with_forbidden_column.tsv"), str(root / "out"))
            self.assertIn("answer_gene_id", str(ctx.exception))
            self.assertIn("forbidden", str(ctx.exception))

    def test_duplicate_panel_id_raises_value_error(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            with self.assertRaises(ValueError) as ctx:
                _freeze(root, [_row(panel_id="dup"), _row(panel_id="dup", query_sequence=PA)])
            self.assertIn("duplicate", str(ctx.exception))
            self.assertIn("dup", str(ctx.exception))

    def test_rank_contribution_non_none_raises_value_error(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            with self.assertRaises(ValueError) as ctx:
                _freeze(root, [_row(), _row(panel_id="row2", query_sequence=PA, rank_contribution="score")])
            self.assertIn("rank_contribution", str(ctx.exception))

    def test_forbidden_columns_bypass_forms_are_rejected(self) -> None:
        # F2: tokenizer bypass forms answer-gene_id / answerGeneId / gene_id
        # itself must all be refused by the normalized substring match.
        for extra in ("answer-gene_id", "answerGeneId", "gene_id", "candidate_gene_id", "score_interval_answer"):
            with self.subTest(extra=extra):
                with tempfile.TemporaryDirectory() as tempdir:
                    root = Path(tempdir)
                    rows = [_row()]
                    rows[0][extra] = "x"
                    path = root / "panel.tsv"
                    header = PANEL_COLUMNS + [extra]
                    path.write_text(
                        "\t".join(header) + "\n"
                        + "\t".join(str(rows[0].get(column, "")) for column in header) + "\n",
                        encoding="utf-8",
                    )
                    with self.assertRaises(ValueError) as ctx:
                        freeze_blast_panel(str(path), str(root / "out"))
                    self.assertIn("forbidden", str(ctx.exception))
                    self.assertIn(extra, str(ctx.exception))

    def test_reference_gene_name_and_rank_contribution_are_not_misdetected(self) -> None:
        # F2: reference_gene_name (-> referencegenename) and rank_contribution
        # must never trip the forbidden-column detector.
        frame = pd.DataFrame([_row()], columns=PANEL_COLUMNS)
        self.assertEqual([], _forbidden_columns(frame))
        # Harmless extra columns pass the check as well.
        frame = pd.DataFrame([[1, 2, 3]], columns=["a", "b", "c"])
        self.assertEqual([], _forbidden_columns(frame))

    def test_unsafe_panel_id_is_rejected(self) -> None:
        # F5: panel_ids must survive safe_id() unchanged; whitespace or other
        # characters fail the freeze.
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            for panel_id in ("bait 1", "bait:1", "bait/1"):
                with self.subTest(panel_id=panel_id):
                    with self.assertRaises(ValueError) as ctx:
                        _freeze(root, [_row(panel_id=panel_id), _row(panel_id="row2", query_sequence=PA)])
                    self.assertIn("unsafe panel_id", str(ctx.exception))

    def test_stable_id_namespace_panel_id_round_trips(self) -> None:
        # Issue #409: a canonical stable ID keeps its database prefix through
        # freeze and verify. Normalizing the namespace away would drop the
        # provenance the accession checks read, so the colon form must survive
        # verbatim in the frozen TSV.
        for panel_id in ("TAIR:AT1G69120", "UniProt:Q9TES1",
                         "RefSeq:NP_001185146.1", "Ensembl:AT1G69120.1"):
            with self.subTest(panel_id=panel_id), tempfile.TemporaryDirectory() as tempdir:
                root = Path(tempdir)
                result = _freeze(root, [_row(panel_id=panel_id),
                                        _row(panel_id="row2", query_sequence=PA)])
                panel = pd.read_csv(result["panel"], sep="\t", dtype=str).fillna("")
                self.assertIn(panel_id, list(panel["panel_id"]))
                verify_blast_panel_freeze(str(result["manifest"]))

    def test_unknown_namespace_panel_id_is_rejected(self) -> None:
        # Issue #409: the colon form is an allowlist, not a general escape.
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            for panel_id in ("gene:AT1G69120", "foo:bar"):
                with self.subTest(panel_id=panel_id):
                    with self.assertRaises(ValueError) as ctx:
                        _freeze(root, [_row(panel_id=panel_id),
                                       _row(panel_id="row2", query_sequence=PA)])
                    self.assertIn("unsafe panel_id", str(ctx.exception))

    def test_known_namespace_with_bad_local_part_is_rejected(self) -> None:
        # Issue #409: an allowed namespace does not excuse an empty or unsafe
        # local part, and a second colon still fails.
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            for panel_id in ("TAIR:", "TAIR:AT1G 69120", "UniProt:Q9TES1:extra", "TAIR:-x"):
                with self.subTest(panel_id=panel_id):
                    with self.assertRaises(ValueError) as ctx:
                        _freeze(root, [_row(panel_id=panel_id),
                                       _row(panel_id="row2", query_sequence=PA)])
                    self.assertIn("unsafe panel_id", str(ctx.exception))

    def test_namespace_match_is_case_insensitive(self) -> None:
        # Issue #409: sources spell the prefix inconsistently (TAIR / tair,
        # UniProt / uniprot); the namespace check must not depend on it.
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            result = _freeze(root, [_row(panel_id="tair:AT1G69120"),
                                    _row(panel_id="UNIPROT:Q9TES1", query_sequence=PA)])
            panel = pd.read_csv(result["panel"], sep="\t", dtype=str).fillna("")
            self.assertEqual(["tair:AT1G69120", "UNIPROT:Q9TES1"],
                             sorted(panel["panel_id"], key=str.casefold))

    def test_panel_ids_colliding_after_normalization_are_rejected(self) -> None:
        # Issue #409: admitting the colon form creates a collision the bare form
        # could not have. "UniProt:Q9TES1" and "UniProt_Q9TES1" are distinct
        # panel_ids that become one query FASTA id downstream, so the freeze
        # must refuse them instead of letting two baits share one record.
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            with self.assertRaises(ValueError) as ctx:
                _freeze(root, [_row(panel_id="UniProt:Q9TES1"),
                               _row(panel_id="UniProt_Q9TES1", query_sequence=PA)])
            self.assertIn("collide once normalized", str(ctx.exception))

    def test_empty_panel_id_is_rejected(self) -> None:
        # F5: a blank panel_id cannot become a query FASTA id or manifest key.
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            with self.assertRaises(ValueError) as ctx:
                _freeze(root, [_row(panel_id=""), _row(panel_id="row2", query_sequence=PA)])
            self.assertIn("empty panel_id", str(ctx.exception))

    def test_duplicate_panel_id_after_strip_is_rejected(self) -> None:
        # F5: duplicates are detected on the stripped value, so "bait1" and
        # " bait1 " cannot both survive into the frozen panel.
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            with self.assertRaises(ValueError) as ctx:
                _freeze(root, [_row(panel_id="bait1"), _row(panel_id=" bait1 ", query_sequence=PA)])
            self.assertIn("duplicate", str(ctx.exception))
            self.assertIn("bait1", str(ctx.exception))


class FreezeAndVerifyTests(unittest.TestCase):
    def test_fixture_panel_freezes_and_verifies(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            result = freeze_blast_panel(str(FIXTURE / "panel_ok.tsv"), str(root / "out"), case_id="route_a")
            self.assertEqual("verified", result["panel_status"])
            self.assertEqual(4, result["row_count"])
            self.assertEqual(3, result["verified_count"])
            self.assertEqual(1, result["review_required_count"])
            self.assertIs(result["interval_accessed"], False)

            panel = pd.read_csv(result["panel"], sep="\t", dtype=str, keep_default_na=False).fillna("")
            self.assertEqual(["bait1", "bait2", "bait3", "bait4"], list(panel["panel_id"]))
            self.assertEqual(["none"] * 4, list(panel["rank_contribution"]))
            self.assertEqual("verified", panel.loc[panel["panel_id"] == "bait1", "review_status"].iloc[0])
            self.assertEqual("short_query", panel.loc[panel["panel_id"] == "bait4", "review_reason"].iloc[0])
            self.assertEqual(_row_sha256(S), panel.loc[panel["panel_id"] == "bait1", "query_sequence_sha256"].iloc[0])

            verified = verify_blast_panel_freeze(str(result["manifest"]))
            self.assertTrue(verified["verified"])
            self.assertIs(verified["interval_accessed"], False)
            self.assertEqual(set(verified["checks"].values()), {"pass"})

    def test_manifest_content(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            result = freeze_blast_panel(str(FIXTURE / "panel_ok.tsv"), str(root / "out"), case_id="route_a")
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual("1.1", manifest["schema_version"])
            self.assertEqual("blast-panel-freeze-1.2.0", manifest["freeze_version"])
            self.assertEqual("blast-panel-freeze", manifest["command"])
            self.assertEqual("route_a", manifest["case_id"])
            self.assertIs(manifest["interval_accessed"], False)
            self.assertIs(manifest["truth_inputs_present"], False)
            self.assertEqual(TRUTH_ABSENCE_ASSERTION, manifest["truth_absence_assertion"])
            self.assertEqual(["panel"], manifest["input_allowlist"])
            self.assertIn("panel", manifest["inputs"])
            self.assertIn("frozen_bait_panel.tsv", manifest["inputs"])
            self.assertEqual(6, len(manifest["policy"]))
            self.assertEqual(dict(DEFAULT_THRESHOLDS), manifest["thresholds"])
            self.assertEqual("verified", manifest["rows"]["bait1"]["review_status"])
            self.assertEqual("verified", manifest["rows"]["bait1"]["verification_status"])
            self.assertEqual("review_required", manifest["rows"]["bait4"]["review_status"])
            self.assertEqual("short_query", manifest["rows"]["bait4"]["review_reason"])
            self.assertEqual("none", manifest["rows"]["bait4"]["rank_contribution"])
            self.assertEqual(64, len(manifest["panel_sha256"]))

    def test_freeze_is_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            first = freeze_blast_panel(str(FIXTURE / "panel_ok.tsv"), str(root / "out1"), case_id="same")
            second = freeze_blast_panel(str(FIXTURE / "panel_ok.tsv"), str(root / "out2"), case_id="same")
            self.assertEqual(Path(first["panel"]).read_bytes(), Path(second["panel"]).read_bytes())
            manifest1 = json.loads(Path(first["manifest"]).read_text(encoding="utf-8"))
            manifest2 = json.loads(Path(second["manifest"]).read_text(encoding="utf-8"))
            manifest1.pop("frozen_at_utc")
            manifest2.pop("frozen_at_utc")
            manifest1.pop("outputs")
            manifest2.pop("outputs")
            # The frozen-panel fingerprint embeds its output path and mtime,
            # which differ per out dir; the input panel fingerprint is byte-
            # identical and stays in the comparison.
            manifest1["inputs"].pop("frozen_bait_panel.tsv")
            manifest2["inputs"].pop("frozen_bait_panel.tsv")
            self.assertEqual(manifest1, manifest2)

    def test_verify_detects_panel_drift_nc5(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            result = freeze_blast_panel(str(FIXTURE / "panel_ok.tsv"), str(root / "out"), case_id="nc5")
            panel_path = Path(result["panel"])
            original = panel_path.read_text(encoding="utf-8")
            panel_path.write_text(original.replace("CDLYQ", "CDLYR", 1), encoding="utf-8")
            with self.assertRaises(ValueError) as ctx:
                verify_blast_panel_freeze(str(result["manifest"]))
            self.assertIn("inputs_fingerprints", str(ctx.exception))

    def test_verify_detects_input_panel_drift(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            panel = _panel_tsv(root, [_row(), _row(panel_id="row2", query_sequence=PA)])
            result = freeze_blast_panel(str(panel), str(root / "out"), case_id="input_drift")
            panel.write_text(panel.read_text(encoding="utf-8") + "\n", encoding="utf-8")
            with self.assertRaises(ValueError) as ctx:
                verify_blast_panel_freeze(str(result["manifest"]))
            self.assertIn("inputs_fingerprints", str(ctx.exception))

    def test_verify_rejects_bad_json_and_wrong_command(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            bad = root / "bad.json"
            bad.write_text("{not json", encoding="utf-8")
            with self.assertRaises(ValueError) as ctx:
                verify_blast_panel_freeze(str(bad))
            self.assertIn("manifest_json", str(ctx.exception))

            freeze = json.loads(Path(_freeze(root, [_row(), _row(panel_id="row2", query_sequence=PA)])["manifest"]).read_text(encoding="utf-8"))
            freeze["command"] = "other"
            wrong = root / "wrong.json"
            wrong.write_text(json.dumps(freeze, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            with self.assertRaises(ValueError) as ctx:
                verify_blast_panel_freeze(str(wrong))
            self.assertIn("command_and_freeze_version", str(ctx.exception))

    def test_verify_rejects_a_manifest_from_the_old_threshold_contract(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            freeze = json.loads(Path(_freeze(
                root, [_row(), _row(panel_id="row2", query_sequence=PA)]
            )["manifest"]).read_text(encoding="utf-8"))
            freeze["freeze_version"] = "blast-panel-freeze-1.1.0"
            freeze["thresholds"]["min_identity"] = 40.0
            old = root / "old-contract.json"
            old.write_text(json.dumps(freeze, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            with self.assertRaises(ValueError) as ctx:
                verify_blast_panel_freeze(str(old))
            self.assertIn("command_and_freeze_version", str(ctx.exception))

    def test_verify_rejects_non_bool_interval_accessed(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            freeze = json.loads(Path(_freeze(root, [_row(), _row(panel_id="row2", query_sequence=PA)])["manifest"]).read_text(encoding="utf-8"))
            freeze["interval_accessed"] = "false"
            wrong = root / "interval.json"
            wrong.write_text(json.dumps(freeze, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            with self.assertRaises(ValueError) as ctx:
                verify_blast_panel_freeze(str(wrong))
            self.assertIn("interval_accessed", str(ctx.exception))

    def test_min_panel_size_small_panel(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            rows = [_row(), _row(panel_id="row2", query_sequence=SHORT20, query_sequence_sha256=_row_sha256(SHORT20))]
            result = _freeze(root, rows)
            self.assertEqual("review_required_small_panel", result["panel_status"])
            self.assertEqual(1, result["verified_count"])
            manifest = json.loads(Path(result["manifest"]).read_text(encoding="utf-8"))
            self.assertEqual("review_required_small_panel", manifest["panel_status"])
            self.assertEqual(2, manifest["min_panel_size"])
            # freeze completes and verifies even when the panel is small.
            verified = verify_blast_panel_freeze(str(result["manifest"]))
            self.assertTrue(verified["verified"])

    def test_panel_arg_override_verifies_against_moved_panel(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            result = freeze_blast_panel(str(FIXTURE / "panel_ok.tsv"), str(root / "out"), case_id="moved")
            moved = root / "moved_panel.tsv"
            moved.write_bytes(Path(result["panel"]).read_bytes())
            verified = verify_blast_panel_freeze(str(result["manifest"]), panel_arg=str(moved))
            self.assertTrue(verified["verified"])
            self.assertEqual(str(moved), verified["panel"])

    def test_verify_rejects_missing_allowed_input_fingerprint(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            freeze = json.loads(Path(_freeze(root, [_row(), _row(panel_id="row2", query_sequence=PA)])["manifest"]).read_text(encoding="utf-8"))
            freeze["inputs"].pop("panel")
            wrong = root / "allowlist.json"
            wrong.write_text(json.dumps(freeze, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            with self.assertRaises(ValueError) as ctx:
                verify_blast_panel_freeze(str(wrong))
            self.assertIn("input_allowlist", str(ctx.exception))

    def test_verify_rejects_out_of_range_threshold(self) -> None:
        # F4: a value-range violation (min_query_coverage=60) is a ValueError.
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            freeze = json.loads(Path(_freeze(root, [_row(), _row(panel_id="row2", query_sequence=PA)])["manifest"]).read_text(encoding="utf-8"))
            freeze["thresholds"]["min_query_coverage"] = 60
            wrong = root / "threshold_range.json"
            wrong.write_text(json.dumps(freeze, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            with self.assertRaises(ValueError) as ctx:
                verify_blast_panel_freeze(str(wrong))
            self.assertIn("thresholds", str(ctx.exception))
            self.assertIn("exactly match", str(ctx.exception))

    def test_verify_rejects_in_range_threshold_mutation(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            freeze = json.loads(Path(_freeze(root, [_row(), _row(panel_id="row2", query_sequence=PA)])["manifest"]).read_text(encoding="utf-8"))
            freeze["thresholds"]["min_identity"] = 50.0
            wrong = root / "thresholds_mutated.json"
            wrong.write_text(json.dumps(freeze, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            with self.assertRaises(ValueError) as ctx:
                verify_blast_panel_freeze(str(wrong))
            self.assertIn("exactly match", str(ctx.exception))

    def test_verify_rejects_missing_or_wrong_threshold_keys(self) -> None:
        # F4: thresholds must be present and carry exactly the default keys.
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            freeze = json.loads(Path(_freeze(root, [_row(), _row(panel_id="row2", query_sequence=PA)])["manifest"]).read_text(encoding="utf-8"))
            freeze.pop("thresholds")
            missing = root / "thresholds_missing.json"
            missing.write_text(json.dumps(freeze, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            with self.assertRaises(ValueError) as ctx:
                verify_blast_panel_freeze(str(missing))
            self.assertIn("thresholds", str(ctx.exception))

            freeze = json.loads(Path(_freeze(root, [_row(), _row(panel_id="row2", query_sequence=PA)])["manifest"]).read_text(encoding="utf-8"))
            freeze["thresholds"] = {"min_query_length": 30}
            wrong = root / "thresholds_wrong_keys.json"
            wrong.write_text(json.dumps(freeze, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            with self.assertRaises(ValueError) as ctx:
                verify_blast_panel_freeze(str(wrong))
            self.assertIn("thresholds", str(ctx.exception))

    def test_verify_rejects_fingerprint_without_path(self) -> None:
        # F4: a fingerprint without a recorded path is a ValueError, never a
        # raw Path("") exception from hashing a directory.
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            freeze = json.loads(Path(_freeze(root, [_row(), _row(panel_id="row2", query_sequence=PA)])["manifest"]).read_text(encoding="utf-8"))
            freeze["inputs"]["panel"] = {"sha256": freeze["inputs"]["panel"]["sha256"]}
            wrong = root / "no_path.json"
            wrong.write_text(json.dumps(freeze, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            with self.assertRaises(ValueError) as ctx:
                verify_blast_panel_freeze(str(wrong))
            self.assertIn("no recorded path", str(ctx.exception))


class PreIntervalFirewallTests(unittest.TestCase):
    def test_candidate_argument_is_rejected_without_opening_the_file(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            with self.assertRaises(ValueError) as ctx:
                _freeze(root, [_row(), _row(panel_id="row2", query_sequence=PA)], leak_check_candidates_arg=str(root / "not-opened.tsv"))
            self.assertIn("pre-interval firewall", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
