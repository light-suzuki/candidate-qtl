"""Offline functional tests. A fake executable tests the subprocess boundary;
its alignments are fixtures, not real BLAST or biological accuracy evidence.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

from candidate_qtl import second_stage_recheck as recheck

FAKE = r'''
import sys,time
from pathlib import Path
a=sys.argv[1:]
if '-version' in a:
 print('FAKE blastp for interface tests, NOT NCBI BLAST');sys.exit(0)
def arg(x): return a[a.index(x)+1]
def fasta(p):
 records=[];name='';parts=[]
 for line in Path(p).read_text().splitlines():
  if line.startswith('>'):
   if name:records.append((name,''.join(parts)))
   name=line[1:].split()[0];parts=[]
  else:parts.append(line)
 if name:records.append((name,''.join(parts)))
 return records
reference='-db' in a or Path(arg('-subject')).name=='reference_proteins.fa'
if MODE=='timeout' and reference:time.sleep(5)
if MODE=='fail' and reference:
 print('fixture failure',file=sys.stderr);sys.exit(2)
if MODE=='empty' and reference:Path(arg('-out')).write_text('');sys.exit(0)
if MODE=='malformed' and reference:Path(arg('-out')).write_text('not an alignment\n');sys.exit(0)
qs=fasta(arg('-query'));ss=fasta(arg('-subject')) if '-subject' in a else [('db|native|id','MABCDE')]
rows=[]
for qid,q in qs:
 for sid,s in ss:
  n=min(len(q),len(s))
  pi=100*sum(x==y for x,y in zip(q[:n],s[:n]))/n
  rows.append('\t'.join(map(str,[qid,sid,pi,n,len(q),len(s),1,n,1,n,1e-5,50,q[:n],s[:n],'fixture title'])))
Path(arg('-out')).write_text('\n'.join(rows)+'\n')
'''


class FreshRecheckTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);self.out=self.root/'new'
        self.request=self.root/'adaptive_synthesis_request.json'
        self.request.write_text(json.dumps({'stage':'adaptive_synthesis_request','candidate_count':3,
            'candidate_gene_ids':['gA','gB','missing'], 'case_id':'synthetic','wet_budget':25}))
        self.proteins=self.root/'proteins.fa';self.proteins.write_text('>gA\nMABCDE\n>gB\nMABCD\n>unused\nMFFFF\n')
        self.references=self.root/'references.fa';self.references.write_text('>ref1 full original title\nMABCDEFGHI\n>ref2 another species\nMABCDE\n')
        self.fake_paths=set()
        self.binary=self.make_binary('normal')
        real_run=recheck.subprocess.run
        real_which=recheck.shutil.which
        def spawn(command, **kwargs):
            if command[0] in self.fake_paths:
                command=[sys.executable, '-S', *command]
            return real_run(command, **kwargs)
        patches=[mock.patch.object(recheck.subprocess,'run',side_effect=spawn),
                 mock.patch.object(recheck.shutil,'which',side_effect=lambda name: name if name in self.fake_paths else real_which(name))]
        for patch in patches:
            patch.start();self.addCleanup(patch.stop)
    def make_binary(self,mode):
        path=self.root/f'blastp-{mode}'
        path.write_text('#!'+sys.executable+' -S\nMODE='+repr(mode)+'\n'+FAKE)
        path.chmod(0o755)
        self.fake_paths.add(str(path))
        return str(path)
    def run_recheck(self,**kw):
        args=dict(phenotype='observed widened stem',species='Example plant',proteins_arg=str(self.proteins),
            reference_proteins_arg=str(self.references),blastp=self.binary)
        args.update(kw)
        return recheck.fresh_recheck(str(self.request),str(self.out),**args)
    def prompts(self):
        return [json.loads(x) for x in (self.out/'second_stage_recheck_prompts.jsonl').read_text().splitlines()]
    def test_fresh_searches_execute_and_first_stage_is_unchanged(self):
        before=self.request.read_bytes(),self.proteins.read_bytes()
        result=self.run_recheck()
        self.assertEqual('completed',result['reference_search']['state'])
        self.assertEqual('completed',result['cohort_search']['state'])
        self.assertEqual(3,result['candidate_count']);self.assertEqual(3,result['prompt_count'])
        self.assertEqual(before,(self.request.read_bytes(),self.proteins.read_bytes()))
        self.assertEqual('none',result['candidate_membership_effect'])
        self.assertEqual(0,result['model_calls']);self.assertEqual('not_executed',result['review_state'])
        self.assertFalse((self.out/'wet_candidates.tsv').exists())
    def test_multiple_references_retained_with_both_coverages(self):
        self.run_recheck();p=json.loads(self.prompts()[0]['user'])
        self.assertEqual({'ref1','ref2'},{r['subject_id'] for r in p['reference_hits']})
        first=next(r for r in p['reference_hits'] if r['subject_id']=='ref1')
        self.assertEqual(1,first['query_coverage']);self.assertAlmostEqual(6/10,first['subject_coverage'])
        self.assertEqual('ref1 full original title',first['subject_title'])
    def test_missing_sequence_is_still_a_research_target(self):
        self.run_recheck();p=json.loads(self.prompts()[-1]['user'])
        self.assertEqual('missing',p['gene_id'])
        self.assertEqual('not_evaluated_missing_sequence',p['reference_search_state'])
        self.assertEqual([],p['reference_hits'])
    def test_copy_comparison_scope_is_not_orthology_claim(self):
        self.run_recheck();p=json.loads(self.prompts()[0]['user'])
        self.assertEqual(['gB'],[r['subject_gene_id'] for r in p['cohort_comparisons']])
        self.assertIn('not_genome_wide',p['cohort_scope'])
        self.assertTrue(all(r['evidence_role']=='homology_context_only' for r in p['cohort_comparisons']))
    def test_duplicate_requests_are_rejected(self):
        self.request.write_text(json.dumps({'stage':'adaptive_synthesis_request','candidate_count':2,'candidate_gene_ids':['gA','gA']}))
        with self.assertRaises(ValueError):self.run_recheck()
        self.assertFalse(self.out.exists())
    def test_wrong_request_type_rejected(self):
        self.request.write_text('{}')
        with self.assertRaises(ValueError):self.run_recheck()
    def test_bad_candidate_count_rejected(self):
        d=json.loads(self.request.read_text());d['candidate_count']=100;self.request.write_text(json.dumps(d))
        with self.assertRaises(ValueError):self.run_recheck()
    def test_existing_run_not_overwritten(self):
        self.run_recheck();before=(self.out/'second_stage_recheck_manifest.json').read_bytes()
        with self.assertRaises(FileExistsError):self.run_recheck()
        self.assertEqual(before,(self.out/'second_stage_recheck_manifest.json').read_bytes())
    def test_missing_blastp_preserves_all_genes_and_emits_prompts(self):
        result=self.run_recheck(blastp=str(self.root/'not-installed'))
        self.assertEqual('not_evaluated_missing_blastp',result['reference_search']['state'])
        self.assertEqual(3,len(self.prompts()));self.assertEqual('fresh_research_partial',result['state'])
    def test_no_proteins_does_not_fake_zero_biology(self):
        result=self.run_recheck(proteins_arg=None)
        self.assertEqual('not_evaluated_missing_proteins',result['sequence_state'])
        self.assertEqual(3,result['candidate_count'])
        self.assertTrue(all(json.loads(p['user'])['reference_search_state']=='not_evaluated_missing_sequence' for p in self.prompts()))
    def test_no_reference_still_compares_remaining_copies(self):
        result=self.run_recheck(reference_proteins_arg=None)
        self.assertEqual('not_evaluated_no_reference_input',result['reference_search']['state'])
        self.assertEqual('completed',result['cohort_search']['state'])
    def test_empty_reference_is_not_searched_no_hit(self):
        self.references.write_text('')
        result=self.run_recheck()
        self.assertEqual('not_evaluated_empty_reference_panel',result['reference_search']['state'])
    def test_missing_reference_is_not_no_hit(self):
        result=self.run_recheck(reference_proteins_arg=str(self.root/'missing.fa'))
        self.assertEqual('not_evaluated_missing_reference_file',result['reference_search']['state'])
    def test_real_empty_search_output_distinguished(self):
        result=self.run_recheck(blastp=self.make_binary('empty'))
        self.assertEqual('completed',result['reference_search']['state'])
        self.assertEqual('searched_no_reported_hit',json.loads(self.prompts()[0]['user'])['reference_search_state'])
        self.assertEqual(3,result['candidate_count'])
    def test_provider_failure_does_not_stop_cohort_or_drop_genes(self):
        result=self.run_recheck(blastp=self.make_binary('fail'))
        self.assertEqual('not_evaluated_blast_failed',result['reference_search']['state'])
        self.assertEqual('completed',result['cohort_search']['state']);self.assertEqual(3,len(self.prompts()))
    def test_timeout_is_distinct(self):
        result=self.run_recheck(blastp=self.make_binary('timeout'),timeout=2)
        self.assertEqual('not_evaluated_timeout',result['reference_search']['state'])
        self.assertEqual('completed',result['cohort_search']['state'])
    def test_malformed_output_is_not_empty(self):
        result=self.run_recheck(blastp=self.make_binary('malformed'))
        self.assertEqual('not_evaluated_invalid_blast_output',result['reference_search']['state'])
        self.assertEqual(3,result['candidate_count'])
    def test_db_settings_are_retrieval_not_candidate_filters(self):
        result=self.run_recheck(reference_proteins_arg=None,protein_db_arg=str(self.root/'db'),max_targets=73,evalue=1)
        command=result['reference_search']['command']
        self.assertEqual('73',command[command.index('-max_target_seqs')+1])
        self.assertEqual('locator_only_not_content_pinned',result['inputs']['protein_db']['fingerprint_state'])
        self.assertEqual(3,result['candidate_count']);self.assertEqual('db|native|id',json.loads(self.prompts()[0]['user'])['reference_hits'][0]['subject_id'])
    def test_no_gene_quota_in_fasta_or_cohort(self):
        result=self.run_recheck(max_targets=1)
        self.assertEqual(2,result['reference_search']['retrieval_max_target_seqs'])
        self.assertEqual(2,result['cohort_search']['retrieval_max_target_seqs'])
        self.assertEqual(2,len(json.loads(self.prompts()[0]['user'])['reference_hits']))
    def test_1247_candidates_not_capped_to_benchmark_or_wet_target(self):
        ids=[f'g{i}' for i in range(1247)]
        self.request.write_text(json.dumps({'stage':'adaptive_synthesis_request','candidate_count':len(ids),'candidate_gene_ids':ids}))
        result=self.run_recheck(proteins_arg=None)
        self.assertEqual(1247,result['candidate_count']);self.assertEqual(1247,result['prompt_count'])
    def test_all_explicit_isoforms_retained(self):
        self.proteins.write_text('>gA.t1\nMABCDE\n>gA.t2\nMABCD\n>gB\nMABCDE\n')
        mapping=self.root/'map.tsv';mapping.write_text('protein_id\tgene_id\ngA.t1\tgA\ngA.t2\tgA\n')
        result=self.run_recheck(protein_gene_map_arg=str(mapping));p=json.loads(self.prompts()[0]['user'])
        self.assertEqual(['gA.t1','gA.t2'],p['protein_ids']);self.assertEqual(3,result['candidate_count'])
        self.assertTrue(any(r['subject_gene_id']=='gA' for r in p['cohort_comparisons']))
    def test_no_suffix_guessing(self):
        self.proteins.write_text('>gA.t1\nMABCDE\n')
        result=self.run_recheck()
        self.assertEqual('not_evaluated_no_mapped_proteins',result['sequence_state'])
    def test_conflicting_map_rejected(self):
        mapping=self.root/'map.tsv';mapping.write_text('protein_id\tgene_id\ngA\tgA\ngA\tgB\n')
        with self.assertRaisesRegex(ValueError,'conflicting'):self.run_recheck(protein_gene_map_arg=str(mapping))
    def test_duplicate_fasta_ids_rejected(self):
        self.proteins.write_text('>gA\nMAAA\n>gA\nMAAA\n')
        with self.assertRaisesRegex(ValueError,'duplicate protein'):self.run_recheck()
    def test_reference_or_db_exclusive(self):
        with self.assertRaisesRegex(ValueError,'OR'):self.run_recheck(protein_db_arg='db')
        self.assertFalse(self.out.exists())
    def test_no_observed_phenotype_fails_before_search(self):
        with self.assertRaisesRegex(ValueError,'observed phenotype'):self.run_recheck(phenotype='')
        self.assertFalse(self.out.exists())
    def test_no_sequences_or_scores_in_prompts(self):
        self.run_recheck()
        for p in self.prompts():
            packet=json.loads(p['user']);self.assertNotIn('sequence',packet)
            self.assertNotIn('rank',packet);self.assertNotIn('score',packet)
            self.assertEqual({'prompt_id','gene_id','system','user'},set(p))
    def test_every_produced_file_fingerprinted(self):
        result=self.run_recheck()
        for name,digest in result['files'].items():self.assertEqual(digest,recheck._sha(self.out/name))
    def test_cli_runs_without_model(self):
        with mock.patch('sys.stdout'):
            code=recheck.main(['--request',str(self.request),'--out-dir',str(self.out),'--phenotype','observed',
                               '--proteins',str(self.proteins),'--reference-proteins',str(self.references),'--blastp',self.binary])
        self.assertEqual(0,code);self.assertTrue((self.out/'reference.hits.tsv').exists())
    def test_guided_hook_obeys_existing_gate(self):
        m=self.root/'handoff.json';m.write_text(json.dumps({'second_stage_state':'skipped_at_or_below_25','adaptive_synthesis_request':str(self.request)}))
        args=argparse.Namespace(second_stage_recheck=True,profile='triple-route',execute=True)
        result={'canonical_route_b_wet':{'wet_handoff_manifest':str(m)}}
        with mock.patch.object(recheck,'fresh_recheck') as run:
            recheck.attach_guided_recheck(result,args);run.assert_not_called()
        self.assertNotIn('second_stage_recheck',result)
    def test_guided_hook_nontriple_and_dry_run(self):
        result={};recheck.attach_guided_recheck(result,argparse.Namespace(profile='screening',execute=True))
        self.assertEqual({},result)
        recheck.attach_guided_recheck(result,argparse.Namespace(profile='triple-route',execute=False))
        self.assertEqual({},result)
    def test_guided_requested_reaches_real_recheck(self):
        m=self.root/'handoff.json';m.write_text(json.dumps({'second_stage_state':'requested','adaptive_synthesis_request':str(self.request)}))
        args=argparse.Namespace(profile='triple-route',execute=True,out_dir=str(self.root/'guided'),
            phenotype='observed',focal_species='plant',proteins=str(self.proteins),protein_db=str(self.root/'db'))
        result={'canonical_route_b_wet':{'wet_handoff_manifest':str(m)}}
        with mock.patch.object(recheck.shutil,'which',return_value=self.binary):
            recheck.attach_guided_recheck(result,args)
        self.assertEqual('completed',result['second_stage_recheck']['reference_search']['state'])
        self.assertTrue(Path(result['second_stage_recheck']['prompts']).is_file())
    def test_same_inputs_reuse_completed_report_without_rerunning(self):
        first=self.run_recheck()
        with mock.patch.object(recheck.subprocess,'run') as run:
            second=self.run_recheck(reuse_existing=True)
        run.assert_not_called();self.assertTrue(second['reused_previous_report'])
        self.assertEqual(first['files'],second['files'])
    def test_reuse_rejects_changed_input(self):
        self.run_recheck();self.proteins.write_text('>gA\nMAAAA\n')
        with self.assertRaisesRegex(ValueError,'changed'):self.run_recheck(reuse_existing=True)
    def test_reuse_rejects_changed_artifact(self):
        self.run_recheck();(self.out/'candidate_status.tsv').write_text('corrupted')
        with self.assertRaisesRegex(ValueError,'artifact'):self.run_recheck(reuse_existing=True)



class AlignmentSummaryTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);self.raw=self.root/'raw.tsv'
    def summarize(self,lines,q='AAAACC',subject_length=10):
        self.raw.write_text('\n'.join(lines)+'\n')
        return recheck._summarise(self.raw,{'q':{'gene_id':'g','protein_id':'g.t1','title':'g','length':len(q),'sequence':q}},
                                  {'s':{'protein_id':'ref','title':'ref title','length':subject_length}})
    def test_duplicate_overlap_and_both_coverages(self):
        first='q\ts\t100\t4\t6\t10\t1\t4\t1\t4\t1e-9\t10\tAAAA\tAAAA\tref'
        second='q\ts\t50\t4\t6\t10\t3\t6\t3\t6\t1e-9\t10\tAACC\tAADD\tref'
        a=self.summarize([first,second,first]);b=self.summarize([second,first])
        self.assertEqual(a,b);self.assertAlmostEqual(4/6*100,a[0]['query_residue_identity'])
        self.assertEqual(1,a[0]['query_coverage']);self.assertAlmostEqual(.6,a[0]['subject_coverage'])
    def test_query_gap_subject_gap_have_distinct_denominators(self):
        line='q\ts\t50\t4\t3\t4\t1\t3\t1\t3\t1e-8\t10\tAA-C\tAAT-\tref'
        r=self.summarize([line],q='AAC',subject_length=4)[0]
        self.assertEqual(1,r['query_coverage']);self.assertEqual(.75,r['subject_coverage'])
        self.assertAlmostEqual(2/3*100,r['query_residue_identity'])
    def test_unknown_query_or_sequence_mismatch_refused(self):
        line='q\ts\t100\t6\t6\t10\t1\t6\t1\t6\t1e-9\t10\tAAAAAA\tAAAAAA\tref'
        with self.assertRaisesRegex(ValueError,'does not match'):self.summarize([line])
    def test_malformed_alignment_refused_not_no_hit(self):
        with self.assertRaisesRegex(ValueError,'invalid BLAST'):self.summarize(['garbage'])


if __name__=='__main__':unittest.main()
