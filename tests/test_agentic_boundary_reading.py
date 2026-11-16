"""Neutral real-CLI transport tests, without model calls or scientific decisions."""
import copy
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import pandas as pd

from scripts.agentic_boundary_material import restore_reading_shard
from scripts.agentic_boundary import merge
from scripts.agentic_boundary_reading import MAX_PAGE_BYTES, MAX_PAGE_LINES, MAX_LINE_CHARS
from scripts.agentic_route_b_material import _encoded, restore_display
from tests import test_agentic_boundary as boundary_fixture
from tests import test_agentic_boundary_material as material_fixture

sha = material_fixture.sha


ROOT = Path(__file__).resolve().parents[1]


class BoundaryReadingTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        boundary, self.request, self.ids = material_fixture.BoundaryMaterialTests().fixture(self.root, 2)
        long_text = ('日本語\\"\n\u2028' * 800) + 'NEUTRAL_FINAL_FACT'
        rows = [{'candidate_gene_id': cid, 'relevance_tier': 'T2_reported_similar_phenomenon',
                 'tier_aggregation': 'strongest_paper_v1',
                 'tier_evidence_json': json.dumps({'reason': long_text}, ensure_ascii=False),
                 'paper_evidence_summary_json': json.dumps([{'reason': long_text}], ensure_ascii=False),
                 'paper_tier_counts_json': '{"T2":1}', 'relevance_evidence_json': '[]'} for cid in self.ids]
        pd.DataFrame(rows).to_csv(boundary, sep='\t', index=False)
        request = json.loads(self.request.read_text())
        request['source_wet_sha256'] = request['candidate_input']['sha256'] = sha(boundary)
        panel = boundary_fixture.AgenticBoundaryTests.paper_panel('neutral', self.ids, 1)
        panel['papers'][0]['reports'][0]['paper_report']['source_excerpt'] = long_text
        request['primary_paper_evidence'] = [panel]
        self.request.write_text(json.dumps(request), encoding='utf-8')
        self.context = self.root / 'context.json'
        self.context.write_text(json.dumps({'case_id': 'case-m', 'phenotype': 'Neutral fixture',
            'candidate_material': {cid: {'annotation': long_text,
                'literal_marker': {'__cqtl_display_string_chunks__': ['original', 'data']}}
                for cid in self.ids + ['H0']}}), encoding='utf-8')
        self.originals = {p: p.read_bytes() for p in [boundary, self.request, self.context]}
        self.prepared = self.root / 'prepared'
        self.nested = self.root / 'nested/worker'
        self.nested.mkdir(parents=True)
        self.cli('prepare', '--request', str(self.request), '--out', str(self.prepared),
                 '--comparison-context', str(self.context), '--max-candidates', '1')
        self.manifest = json.loads((self.prepared / 'manifest.json').read_text())

    def cli(self, *args):
        env = {**os.environ, 'PYTHONPATH': os.pathsep.join([str(ROOT / 'src'), str(ROOT)]),
               'PYTHONIOENCODING': 'utf-8'}
        return subprocess.run([sys.executable, str(ROOT / 'scripts/agentic_boundary.py'), *args],
                              cwd=self.nested, env=env, capture_output=True, text=True, encoding='utf-8', check=True)

    def test_cli_pages_restore_all_fields_and_reach_real_consumer_from_nested_worker(self):
        delivery = self.manifest['reading_pages']
        source_names = [s['source_file'] for s in delivery['sources']]
        self.assertEqual(len(source_names), len(set(source_names)))
        self.assertEqual(4, len(source_names))  # Two views plus one shared catalog and index.
        for source in delivery['sources']:
            pages = [(self.prepared / p['file']).read_bytes() for p in source['pages']]
            original = json.loads((self.prepared / source['source_file']).read_bytes())
            self.assertEqual(original, restore_display(b''.join(pages), _encoded(original)))
            for page in pages:
                lines = page.decode('utf-8').splitlines()
                self.assertLessEqual(len(page), MAX_PAGE_BYTES)
                self.assertLessEqual(len(lines), MAX_PAGE_LINES)
                self.assertLessEqual(max(map(len, lines)), MAX_LINE_CHARS)
                # Simulated line-number overhead stays well below the observed ~56 KiB return ceiling.
                numbered = '\n'.join(f'{n}: {line}' for n, line in enumerate(lines, 1)).encode('utf-8')
                self.assertLess(len(numbered), 16 * 1024)
        replies = []
        for shard, task in zip(self.manifest['shards'], delivery['tasks']):
            task_path = self.prepared / task['file']
            reading_root = (task_path.parent / '..').resolve()
            self.assertEqual(self.prepared.resolve(), reading_root)
            task_text = task_path.read_text()
            self.assertIn('candidate_input AND primary_paper_evidence', task_text)
            for name in [self.manifest['material_catalog']['file'], self.manifest['comparison_index_artifact']['file'],
                         next(v['file'] for v in self.manifest['reading_views'] if v['shard_id'] == shard['shard_id'])]:
                source = next(s for s in delivery['sources'] if s['source_file'] == name)
                self.assertIn(source['pages'][0]['file'], task_text)
                self.assertIn(source['pages'][-1]['file'], task_text)
                self.assertTrue(all((reading_root / p['file']).is_file() for p in source['pages']))
            canonical = restore_reading_shard(self.prepared, self.manifest, shard['shard_id'])
            self.assertIn('NEUTRAL_FINAL_FACT', canonical['candidate_input'][0]['tier_evidence_json'])
            self.assertTrue(canonical['primary_paper_evidence'])
            path = self.root / (shard['shard_id'] + '-response.json')
            path.write_text(json.dumps({'case_id': 'case-m', 'source_wet_sha256': self.manifest['source_wet_sha256'],
                'comparison_context_sha256': self.manifest['comparison_context']['sha256'],
                'shard_id': shard['shard_id'], 'executor_id': 'fixture', 'executor_model': 'fixture/model',
                'status': 'completed', 'responses': [{'candidate_gene_id': cid, 'decision': 'not_evaluated',
                'reason': 'Neutral transport fixture; no scientific judgment.'} for cid in shard['candidate_gene_ids']]}), encoding='utf-8')
            replies += ['--response', str(path)]
        output = self.root / 'merged.tsv'
        self.cli('merge', '--prepared', str(self.prepared), *replies, '--role', 'proposal',
                 '--executor-id', 'fixture', '--executor-model', 'fixture/model', '--out', str(output))
        self.assertEqual(self.ids, pd.read_csv(output, sep='\t')['candidate_gene_id'].tolist())
        self.assertEqual(self.originals, {p: p.read_bytes() for p in self.originals})

    def test_missing_tampered_and_reordered_pages_are_rejected(self):
        source = self.manifest['reading_pages']['sources'][0]
        page = self.prepared / source['pages'][0]['file']
        original = page.read_bytes()
        for bad in (None, original + b' '):
            with self.subTest(missing=bad is None):
                if bad is None:
                    page.unlink()
                else:
                    page.write_bytes(bad)
                with self.assertRaisesRegex(ValueError, 'boundary reading page/task'):
                    restore_reading_shard(self.prepared, self.manifest, self.manifest['shards'][0]['shard_id'])
                page.write_bytes(original)
        changed = copy.deepcopy(self.manifest)
        changed['reading_pages']['sources'][0]['pages'].reverse()
        with self.assertRaisesRegex(ValueError, 'order/roster'):
            restore_reading_shard(self.prepared, changed, changed['shards'][0]['shard_id'])

    def test_self_consistent_display_and_task_tamper_reject_and_legacy_restores(self):
        changed = copy.deepcopy(self.manifest)
        source = changed['reading_pages']['sources'][0]
        page_info = next(p for p in source['pages'] if b'NEUTRAL_FINAL_FACT' in (self.prepared / p['file']).read_bytes())
        path = self.prepared / page_info['file']
        original = path.read_bytes()
        path.write_bytes(original.replace(b'NEUTRAL_FINAL_FACT', b'CHANGED_FINAL_FACT'))
        page_info['sha256'] = sha(path)
        source['display_sha256'] = hashlib.sha256(b''.join((self.prepared / p['file']).read_bytes() for p in source['pages'])).hexdigest()
        with self.assertRaisesRegex(ValueError, 'display differs from canonical'):
            restore_reading_shard(self.prepared, changed, changed['shards'][0]['shard_id'])
        path.write_bytes(original)
        task = changed['reading_pages']['tasks'][0]
        task_path = self.prepared / task['file']
        task_original = task_path.read_bytes()
        task_path.write_bytes(task_original.replace(b'Read every numbered page', b'Skip every numbered page'))
        task['sha256'] = sha(task_path)
        changed['reading_pages']['sources'] = self.manifest['reading_pages']['sources']
        with self.assertRaisesRegex(ValueError, 'task differs'):
            restore_reading_shard(self.prepared, changed, changed['shards'][0]['shard_id'])
        legacy = copy.deepcopy(self.manifest)
        pages = (self.prepared / 'reading_pages').resolve()
        self.assertTrue(pages.is_relative_to(self.root.resolve()))
        shutil.rmtree(pages)
        for task in legacy['reading_pages']['tasks']:
            (self.prepared / task['file']).unlink()
        del legacy['reading_pages']
        del legacy['reading_delivery_contract']
        self.assertEqual(self.ids[:1], restore_reading_shard(self.prepared, legacy, legacy['shards'][0]['shard_id'])['candidate_gene_ids'])

    def test_missing_or_null_metadata_cannot_downgrade_page_generation(self):
        page = self.prepared / self.manifest['reading_pages']['sources'][0]['pages'][0]['file']
        page.write_bytes(page.read_bytes() + b' ')
        for null in (False, True):
            for remove_marker in (False, True):
                for remove_catalog in (False, True):
                    with self.subTest(null=null, remove_marker=remove_marker, remove_catalog=remove_catalog):
                        changed = copy.deepcopy(self.manifest)
                        if null:
                            changed['reading_pages'] = None
                        else:
                            del changed['reading_pages']
                        if remove_marker:
                            del changed['reading_delivery_contract']
                        if remove_catalog:
                            del changed['material_catalog']
                        (self.prepared / 'manifest.json').write_text(json.dumps(changed), encoding='utf-8')
                        with self.assertRaisesRegex(ValueError, 'metadata missing'):
                            restore_reading_shard(self.prepared, changed, changed['shards'][0]['shard_id'])
                        with self.assertRaisesRegex(ValueError, 'metadata missing'):
                            merge(str(self.prepared), [], str(self.root / 'rejected.tsv'), role='proposal',
                                  executor_id='fixture', executor_model='fixture/model')
        changed = copy.deepcopy(self.manifest)
        del changed['reading_pages']
        del changed['reading_delivery_contract']
        pages = (self.prepared / 'reading_pages').resolve()
        self.assertTrue(pages.is_relative_to(self.root.resolve()))
        shutil.rmtree(pages)
        # A surviving task alone still identifies a page-bearing generation.
        with self.assertRaisesRegex(ValueError, 'metadata missing'):
            restore_reading_shard(self.prepared, changed, changed['shards'][0]['shard_id'])
        for task in self.manifest['reading_pages']['tasks']:
            (self.prepared / task['file']).unlink()
        changed['reading_delivery_contract'] = self.manifest['reading_delivery_contract']
        # The explicit marker still rejects missing metadata if all pages/tasks were lost.
        with self.assertRaisesRegex(ValueError, 'metadata missing'):
            restore_reading_shard(self.prepared, changed, changed['shards'][0]['shard_id'])


if __name__ == '__main__':
    unittest.main()
