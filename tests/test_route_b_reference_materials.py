"""Multi-reference investigation material through the real semantic consumer."""
from __future__ import annotations
import csv
import hashlib
import itertools
import json
import tempfile
import unittest
from pathlib import Path
from candidate_qtl.route_b_reference_materials import load_reference_materials
from candidate_qtl.route_b_semantic_bundle import build_route_b_semantic_bundle
from candidate_qtl.route_b_blast_evidence import run_local_blastp
from candidate_qtl.route_b_prepare import route_b_prepare


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def table(p, rows):
    with p.open('w', encoding='utf-8', newline='') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        w.writeheader(); w.writerows(rows)
    return p


class ReferenceMaterialTests(unittest.TestCase):
    def fixture(self, root, n=2):
        atlas=table(root/'interval_gene_atlas.tsv', [
            {'gene_id':f'g{i}', 'orthology_reference_gene':'refA;refZ' if i==0 else '',
             'orthology_relationship':'best_hit_homology', 'family_id':'', 'orthogroup':''}
            for i in range(n)])
        corpus=table(root/'reference.tsv', [
            {'reference_gene_id':'refA','reference_species':'species A','function_description':'uncharacterized',
             'source_id':'DB:A','source_url':'','mutant_phenotype':'','pathway':''},
            {'reference_gene_id':'refZ','reference_species':'species Z','function_description':'measured enzyme activity',
             'source_id':'PMID:123','source_url':'','mutant_phenotype':'short stem','pathway':'growth'},
        ])
        m={'command':'route-b-prepare','frozen_output_fingerprints':{'interval_gene_atlas':{'sha256':sha(atlas)}},
           'providers':{'reference_function':{'available':True,'path':str(corpus),'sha256':sha(corpus)}},
           'outputs':{'reference_function':str(corpus)}}
        (root/'route_b_preparation_manifest.json').write_text(json.dumps(m),encoding='utf-8')
        return atlas,corpus

    def test_n1000_all_refs_reach_real_prompt_without_membership_inflation(self):
        with tempfile.TemporaryDirectory() as t:
            r=Path(t);atlas,_=self.fixture(r,1000)
            out=build_route_b_semantic_bundle(str(atlas),str(r/'out'),phenotype='short stem')
            bundle=json.loads(Path(out['bundle']).read_text())
            self.assertEqual(1000,bundle['interval_gene_count'])
            self.assertEqual(1000,sum(len(u['member_interval_genes']) for u in bundle['units']))
            target=next(u for u in bundle['units'] if 'g0' in u['member_interval_genes'])
            cards=target['reference_material_cards']
            self.assertEqual({'refA','refZ'},{c['reference_gene'] for c in cards})
            self.assertIn('PMID:123',target['stable_source_ids'])
            prompt=Path(out['prompts']).read_text()
            self.assertIn('measured enzyme activity',prompt)
            self.assertIn('short stem',prompt)
            self.assertTrue(all(c['rank_contribution']=='none' for c in cards))
            self.assertFalse(any('verified_support' in c for c in cards))
            for c in cards:
                self.assertEqual('PMID:123' if c['reference_gene']=='refZ' else 'DB:A',c['source_id'])

    def test_standalone_atlas_does_not_guess_resource(self):
        with tempfile.TemporaryDirectory() as t:
            r=Path(t);p=table(r/'atlas.tsv',[{'gene_id':'g','orthology_reference_gene':'refZ'}])
            cards,meta=load_reference_materials(p,[{'gene_id':'g'}])
            self.assertEqual({},cards);self.assertEqual('not_supplied',meta['status'])

    def test_modified_corpus_is_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            r=Path(t);a,c=self.fixture(r);c.write_text(c.read_text()+'\n')
            with self.assertRaisesRegex(ValueError,'fingerprint'):
                build_route_b_semantic_bundle(str(a),str(r/'out'),phenotype='short stem')

    def test_modified_atlas_is_rejected(self):
        with tempfile.TemporaryDirectory() as t:
            r=Path(t);a,_=self.fixture(r);a.write_text(a.read_text()+'\n')
            with self.assertRaisesRegex(ValueError,'fingerprint'):
                load_reference_materials(a,[])

    def test_exact_identity_join_no_symbol_guess(self):
        with tempfile.TemporaryDirectory() as t:
            r=Path(t);a,_=self.fixture(r)
            cards,_=load_reference_materials(a,[{'gene_id':'g','orthology_reference_gene':'REFZ;aliasZ',
                'family_index_best_subject':'sp|REFZ|KNOWN_SYMBOL',
                'family_index_subjects':'tr|refz|CASE_MISMATCH'}])
            self.assertEqual({},cards)

    def test_route_b_prepare_family_endpoints_reach_source_cards_and_semantic_prompt(self):
        with tempfile.TemporaryDirectory() as t:
            r=Path(t)
            candidates=table(r/'candidates.tsv',[
                {'gene_id':'copyA','chrom':'chr1','start':'100','end':'200'},
                {'gene_id':'copyB','chrom':'chr1','start':'300','end':'400'},
            ])
            family=table(r/'family.tsv',[
                {'gene_id':'copyA','best_subject':'sp|P11111|KINASE_A','best_subject_symbol':'KINASE_A',
                 'best_subject_stem':'KINASE','family_subjects':'sp|P11111|KINASE_A;tr|P22222|KINASE_B',
                 'assignment_status':'assigned','best_identity':'88','best_query_coverage':'0.9',
                 'best_evalue':'1e-40','best_bitscore':'240','second_bitscore':'220','top2_bitscore_ratio':'0.92',
                 'symbol_source_rank':'2','symbol_source_bitscore_ratio':'0.97'},
                {'gene_id':'copyB','best_subject':'tr|P22222|KINASE_B','best_subject_symbol':'KINASE_B',
                 'best_subject_stem':'KINASE','family_subjects':'tr|P22222|KINASE_B',
                 'assignment_status':'assigned','best_identity':'84','best_query_coverage':'0.8',
                 'best_evalue':'1e-30','best_bitscore':'200','second_bitscore':'0','top2_bitscore_ratio':'0',
                 'symbol_source_rank':'1','symbol_source_bitscore_ratio':'1'},
            ])
            orthology=table(r/'orthology.tsv',[
                {'gene_id':'copyA','reference_gene_id':'ORTH_REF','relationship':'ortholog','provider':'fixture'}])
            corpus=table(r/'corpus.tsv',[
                {'reference_gene_id':'P11111','reference_species':'Species X','function_description':'P11111 activity card one',
                 'source_id':'SRC:P111-A','source_url':'','mutant_phenotype':'','pathway':''},
                {'reference_gene_id':'P11111','reference_species':'Species X','function_description':'P11111 activity card two',
                 'source_id':'SRC:P111-B','source_url':'','mutant_phenotype':'','pathway':''},
                {'reference_gene_id':'P22222','reference_species':'Species Y','function_description':'P22222 pathway observation',
                 'source_id':'SRC:P222','source_url':'','mutant_phenotype':'','pathway':''},
                {'reference_gene_id':'ORTH_REF','reference_species':'Species Z','function_description':'orthology-only observation',
                 'source_id':'SRC:ORTH','source_url':'','mutant_phenotype':'','pathway':''},
                {'reference_gene_id':'p22222','reference_species':'wrong case','function_description':'must not match',
                 'source_id':'SRC:CASE','source_url':'','mutant_phenotype':'','pathway':''},
            ])
            prepared=route_b_prepare(str(candidates),str(r/'prepared'),family_index_arg=str(family),
                orthology_arg=str(orthology),reference_function_args=[str(corpus)],auto_local_db=False)
            atlas=Path(prepared['interval_gene_atlas'])
            atlas_rows=list(csv.DictReader(atlas.open(encoding='utf-8'),delimiter='\t'))
            self.assertEqual(2,len(atlas_rows))
            out=build_route_b_semantic_bundle(str(atlas),str(r/'semantic'),phenotype='stem elongation')
            bundle=json.loads(Path(out['bundle']).read_text(encoding='utf-8'))
            by_gene={gene_id:unit for unit in bundle['units'] for gene_id in unit['member_interval_genes']}
            self.assertEqual({'copyA','copyB'},set(by_gene))
            cards=by_gene['copyA']['reference_material_cards']
            self.assertEqual({'copyA','copyB'},{c['interval_gene_id'] for c in cards})
            cards_a=[c for c in cards if c['interval_gene_id']=='copyA']
            cards_b=[c for c in cards if c['interval_gene_id']=='copyB']
            self.assertEqual({'P11111','P22222','ORTH_REF'},{c['reference_gene'] for c in cards_a})
            self.assertEqual({'P22222'},{c['reference_gene'] for c in cards_b})
            self.assertEqual({'SRC:P111-A','SRC:P111-B'},
                             {c['source_id'] for c in cards_a if c['reference_gene']=='P11111'})
            self.assertFalse(any(c['source_id']=='SRC:CASE' for c in cards))
            self.assertTrue(all(c['rank_contribution']=='none' and
                                c['material_status']=='investigation_context_not_verified_support'
                                for c in cards))
            prompt=Path(out['prompts']).read_text(encoding='utf-8')
            for text in ('P11111 activity card one','P11111 activity card two',
                         'P22222 pathway observation','orthology-only observation','SRC:P222'):
                self.assertIn(text,prompt)
            self.assertNotIn('must not match',prompt)

    def test_all_retrieved_qualified_references_order_invariant(self):
        import subprocess
        with tempfile.TemporaryDirectory() as t:
            r=Path(t)
            rows=['\t'.join(map(str,['g',ref,100,4,4,4,1,4,1,4,1e-20,score,'AAAA','AAAA','title']))
                  for ref,score in [('uncharacterized',300),('functional',299),('tie',300)]]
            outputs=[]
            for order in itertools.permutations(rows):
                def run(cmd,**kw):
                    Path(cmd[cmd.index('-out')+1]).write_text('\n'.join(order)+'\n')
                    return subprocess.CompletedProcess(cmd,0,'','')
                data,meta=run_local_blastp(r/'q',r/'db',r,1,30,.5,1e-5,
                    run_tool=run,find_binary=lambda _: 'blastp',normalise_subject=lambda x:x)
                outputs.append(data)
                self.assertEqual(3,meta['retained_reference_pair_count'])
                self.assertEqual(1,meta['genes_with_a_usable_hit'])
                self.assertTrue(all(x['relationship']=='best_hit_homology' for x in data))
            self.assertTrue(all(x==outputs[0] for x in outputs))
            self.assertEqual({'uncharacterized','functional','tie'},{x['reference_gene_id'] for x in outputs[0]})


if __name__=='__main__': unittest.main()
