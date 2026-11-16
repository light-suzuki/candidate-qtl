"""BENCH_ONLY: the saved Sol shared-material reader, portable and retry-scoped.

One lossless chunk per invocation; no model calls or scientific default decisions.
An ACK attests transport reading, not semantic correctness. Historical cursors
and scientific responses must stay untouched; initialize a new cursor instead.
"""
import argparse,hashlib,json,sys,uuid
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'scripts'),str(ROOT/'src')]
from agentic_route_b_material import restore_reading_shard
LIMIT=8000

def local_path(value):
    """Interpret recorded paths natively without modifying frozen metadata."""
    text=str(value)
    if sys.platform=='win32' and text.startswith('/mnt/') and len(text)>7 and text[6]=='/':
        text=text[5].upper()+':/'+text[7:]
    return Path(text)

def raw(p):return local_path(p).read_bytes()
def digest(b):return hashlib.sha256(b).hexdigest()
def read(p):return json.loads(raw(p).decode('utf-8-sig'))
def encoded(v):return (json.dumps(v,ensure_ascii=False,sort_keys=True,separators=(',',':'))+'\n').encode('utf-8')
def save(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    temp=p.with_name(p.name+'.tmp');temp.write_bytes(encoded(v));temp.replace(p)
def safe(root,name):
    p=(root/name).resolve()
    if not p.is_relative_to(root):raise ValueError('material outside frozen reading view')
    return p

def material_ids(v,known):
    found=set()
    if isinstance(v,dict):
        for k,x in v.items():
            if k.endswith('material_id') and isinstance(x,str) and x not in known:raise ValueError('unknown material reference')
            found.update(material_ids(x,known))
    elif isinstance(v,list):
        for x in v:found.update(material_ids(x,known))
    elif isinstance(v,str) and v in known:found.add(v)
    return found

def output(v):
    line=json.dumps(v,ensure_ascii=False,separators=(',',':'))
    if len(line)+1>LIMIT:raise ValueError('output exceeds hard 8000-character bound')
    sys.stdout.write(line+'\n')

def scope(s):
    if s.get('reading_shard'):return {'reading_shard':s['reading_shard']['file']}
    return {'range':[s['start'],s['stop']]}

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('command',choices=['init','next','status','commit','export'])
p.add_argument('--state',type=Path,required=True)
p.add_argument('--view-root',type=Path);p.add_argument('--start',type=int);p.add_argument('--stop',type=int)
p.add_argument('--reading-shard',help='init with exactly this manifest reading assignment, preserving its order; exclusive with start/stop')
p.add_argument('--expected-source-sha256');p.add_argument('--executor');p.add_argument('--ack')
p.add_argument('--checkpoint',type=Path,help='validated canonical partial merge; deliver only unfinished IDs')
p.add_argument('--allow-partial',action='store_true',help='export completed decisions and explicit unread states')
p.add_argument('--model',help='actual requested model for this new execution; never inferred')
p.add_argument('--decision-file',type=Path);p.add_argument('--attest-genuine-authorship',action='store_true')
p.add_argument('--authorship-note');p.add_argument('--out',type=Path);p.add_argument('--validation-only',action='store_true')
a=p.parse_args();folder=a.state.resolve();statefile=folder/'state.json'
if a.command=='init':
    if folder.exists():raise ValueError('state directory already exists; resume existing cursor or use a new path')
    if a.view_root is None or not a.executor or not a.expected_source_sha256:raise ValueError('init requires view-root/executor/expected-source-sha256')
    viewroot=a.view_root.resolve();manifest=read(viewroot/'prepare_manifest.json')
    if not manifest['reading_view']['exact_roundtrip_verified']:raise ValueError('view lacks exact roundtrip proof')
    if manifest['source_sha256']!=a.expected_source_sha256:raise ValueError('frozen source hash differs')
    source=local_path(manifest['source_prompts'])
    if digest(raw(source))!=manifest['source_sha256']:raise ValueError('original prompt bytes changed')
    ordered=[pid for sh in manifest['shards'] for pid in sh['prompt_ids']]
    if len(ordered)!=len(set(ordered)):raise ValueError('duplicate canonical prompt ID')
    reading_shards=manifest['reading_view']['shards']
    selected=None
    if a.reading_shard is not None:
        if a.start is not None or a.stop is not None:raise ValueError('reading-shard and start/stop are mutually exclusive')
        reading_shards=[sh for sh in reading_shards if sh['file']==a.reading_shard]
        if len(reading_shards)!=1:raise ValueError('unknown or duplicate reading shard')
        selected=reading_shards[0]
        positions={pid:i for i,pid in enumerate(ordered)}
        roster=selected['prompt_ids']
        if not roster or len(roster)!=len(set(roster)) or any(pid not in positions for pid in roster):raise ValueError('invalid reading shard prompt roster')
        wanted={pid:positions[pid] for pid in roster}
    else:
        if a.start is None or a.stop is None or not (0<=a.start<a.stop<=len(ordered)):raise ValueError('invalid range')
        wanted={pid:i for i,pid in enumerate(ordered) if a.start<=i<a.stop}
    units={}
    checkpoint_binding=None
    if a.checkpoint:
        from practical_packet_run import _route_b_execution_context
        from candidate_qtl.route_b_response_validation import route_b_response_completed
        checkpoint=a.checkpoint.resolve()
        execution_path=checkpoint.with_name(checkpoint.stem+'_manifest.json')
        folder.mkdir(parents=True)
        _route_b_execution_context(execution_path,{'case_id':manifest['case_id'],'bundle':{'prompts':str(source)}},checkpoint,folder)
        healthy={row['prompt_id'] for row in read(checkpoint)['candidates'] if route_b_response_completed(row)}
        wanted={pid:index for pid,index in wanted.items() if pid not in healthy}
        checkpoint_binding={'response':str(checkpoint),'response_sha256':digest(raw(checkpoint)),
                            'execution_manifest':str(execution_path),'execution_sha256':digest(raw(execution_path))}
    for sh in reading_shards:
        # An all-completed checkpoint still validates the selected assignment.
        if selected is None and not wanted.keys() & set(sh['prompt_ids']):continue
        assignmentpath=safe(viewroot,sh['file']);assignment=read(assignmentpath)
        if digest(raw(assignmentpath))!=sh['sha256']:raise ValueError('assignment bytes changed')
        # Existing canonical restoration checks every table, linked raw value and binding.
        restore_reading_shard(viewroot,sh['file'])
        materials={m['material_id']:m for m in assignment['materials']}
        for view in assignment['prompts']:
            pid=view['prompt']['prompt_id']
            if pid not in wanted:continue
            links=[materials[mid] for mid in sorted(material_ids(view,materials))]
            for link in links:
                if digest(raw(safe(viewroot,link['file'])))!=link['sha256']:raise ValueError('raw material hash changed')
            index=wanted[pid];unitfile=folder/'unit_inputs'/f'{index:04d}.json'
            if index in units:raise ValueError('duplicate assigned prompt ID')
            wrapper={'contract':assignment['contract'],'source_sha256':manifest['source_sha256'],'instructions':assignment['instructions'],'prompts':[view],'materials':links}
            save(unitfile,wrapper)
            priority=lambda m:0 if m['material_id'].startswith('system-') else 1 if m['material_id'].startswith('user-prefix-') else 2
            first=[m for m in links if priority(m)<2];rest=[m for m in links if priority(m)==2]
            artifact=lambda m:{'path':str(safe(viewroot,m['file'])),'sha256':m['sha256'],'kind':'complete_shared_material','material_id':m['material_id'],'dedupe':True}
            pieces=[artifact(m) for m in sorted(first,key=priority)]
            pieces.append({'path':str(unitfile),'sha256':digest(raw(unitfile)),'kind':'complete_unit_tables_bindings','dedupe':False})
            pieces.extend(artifact(m) for m in rest)
            units[index]={'source_index':index,'prompt_id':pid,'route_b_unit_id':view['prompt']['route_b_unit_id'],'canonical_prompt_sha256':view['prompt_sha256'],'reading_assignment':sh['file'],'reading_assignment_sha256':sh['sha256'],'artifacts':pieces}
    if set(units)!=set(wanted.values()):raise ValueError('assigned prompt coverage mismatch')
    state={'schema_version':1,'case_id':manifest['case_id'],'executor':a.executor,'validation_only':a.validation_only,'view_root':str(viewroot),'manifest_sha256':digest(raw(viewroot/'prepare_manifest.json')),'source_sha256':manifest['source_sha256'],'start':a.start,'stop':a.stop,'units':[units[i] for i in wanted.values()],'unit_cursor':0,'artifact_cursor':0,'char_cursor':0,'pending':None,'acknowledged_chunks':[],'shared_fully_acknowledged':{},'unit_ready':False,'decisions':[],'checkpoint':checkpoint_binding,'model':a.model}
    if selected is not None:
        state['reading_shard']={'file':selected['file'],'sha256':selected['sha256'],'prompt_ids':selected['prompt_ids']}
    save(statefile,state);output({'state':'initialized','case_id':state['case_id'],**scope(state),'units':len(state['units']),'max_output_characters':LIMIT,'validation_only':a.validation_only});sys.exit()
s=read(statefile)
if a.executor and a.executor!=s['executor']:raise ValueError('executor differs; never reuse another worker shared-read state')

def status():return {'state':'all_decisions_checkpointed' if s['unit_cursor']==len(s['units']) else 'unit_ready_requires_genuine_decision' if s['unit_ready'] else 'awaiting_chunk_ack' if s['pending'] else 'reading','case_id':s['case_id'],**scope(s),'decision_count':len(s['decisions']),'source_index':s['units'][s['unit_cursor']]['source_index'] if s['unit_cursor']<len(s['units']) else None,'pending_token':s['pending']['token'] if s['pending'] else None,'acknowledged_chunk_count':len(s['acknowledged_chunks'])}
if a.command=='status':output(status());sys.exit()
if a.command=='next':
    if a.ack:
        pending=s['pending']
        if not pending or pending['token']!=a.ack:raise ValueError('ack token does not match last delivered chunk')
        delivered=read(folder/pending['delivery_file'])
        if digest(delivered['content'].encode('utf-8'))!=pending['content_sha256']:raise ValueError('delivery audit changed')
        s['acknowledged_chunks'].append(pending);s['char_cursor']=pending['end_character'];s['pending']=None
        save(statefile,s)
    elif s['pending']:
        # Repeat the same single bounded chunk; never silently skip unacknowledged text.
        output(read(folder/s['pending']['delivery_file']));sys.exit()
    if s['unit_cursor']==len(s['units']) or s['unit_ready']:output(status());sys.exit()
    unit=s['units'][s['unit_cursor']]
    while s['artifact_cursor']<len(unit['artifacts']):
        artifact=unit['artifacts'][s['artifact_cursor']]
        if artifact['dedupe'] and artifact['sha256'] in s['shared_fully_acknowledged']:
            s['artifact_cursor']+=1;s['char_cursor']=0;continue
        data=raw(artifact['path'])
        if digest(data)!=artifact['sha256']:raise ValueError('frozen input changed during reading')
        text=data.decode('utf-8');begin=s['char_cursor']
        if begin==len(text):
            if artifact['dedupe']:s['shared_fully_acknowledged'][artifact['sha256']]={'material_id':artifact['material_id'],'first_completed_source_index':unit['source_index']}
            s['artifact_cursor']+=1;s['char_cursor']=0;continue
        number=len(s['acknowledged_chunks']);token=uuid.uuid4().hex
        common={'state':'READ_ONE_CHUNK','source_index':unit['source_index'],'prompt_id':unit['prompt_id'],'route_b_unit_id':unit['route_b_unit_id'],'kind':artifact['kind'],'raw_artifact_sha256':artifact['sha256'],'ack_token':token,'begin_character':begin,'total_characters':len(text),'instruction':'Read this entire content; next invocation acknowledges ONLY this chunk. Never batch multiple next invocations into one tool output.'}
        lo,hi=1,min(len(text)-begin,LIMIT)
        while lo<hi:
            mid=(lo+hi+1)//2;trial=dict(common,end_character=begin+mid,content=text[begin:begin+mid])
            if len(json.dumps(trial,ensure_ascii=False,separators=(',',':')))<=LIMIT-1:lo=mid
            else:hi=mid-1
        end=begin+lo;message=dict(common,end_character=end,content=text[begin:end])
        file=f'deliveries/{number:06d}.json';save(folder/file,message)
        s['pending']={'token':token,'source_index':unit['source_index'],'artifact_sha256':artifact['sha256'],'begin_character':begin,'end_character':end,'content_sha256':digest(text[begin:end].encode('utf-8')),'delivery_file':file,'output_characters':len(json.dumps(message,ensure_ascii=False,separators=(',',':')))+1}
        save(statefile,s);output(message);sys.exit()
    s['unit_ready']=True;save(statefile,s)
    output(dict(status(),prompt_id=unit['prompt_id'],route_b_unit_id=unit['route_b_unit_id'],instruction='ALL pieces delivered and acknowledged. Independently interpret THIS unit, write its genuine decision, then commit. No default-empty generator; no advance without authored decision.'));sys.exit()
if a.command=='commit':
    if not s['unit_ready'] or s['pending']:raise ValueError('unit not ready; unread pieces remain')
    if not a.attest_genuine_authorship or not a.authorship_note or not a.decision_file or not a.executor:raise ValueError('commit requires decision-file/executor/genuine-authorship attestation/note')
    if s['validation_only']:raise ValueError('mechanical validation cursor cannot create scientific decisions')
    unit=s['units'][s['unit_cursor']];decision=read(a.decision_file)
    if decision.get('prompt_id')!=unit['prompt_id'] or decision.get('route_b_unit_id')!=unit['route_b_unit_id']:raise ValueError('decision IDs differ')
    if decision.get('status')!='completed' or not isinstance(decision.get('raw_text'),str):raise ValueError('require actual authored completed raw JSON decision; no unread defaults')
    claims=json.loads(decision['raw_text'])
    if not isinstance(claims,dict) or not isinstance(claims.get('claims'),list) or any(not isinstance(c,dict) for c in claims['claims']):raise ValueError('claims must be literal array of objects')
    if 'model_json' in decision and decision['model_json']!=claims:raise ValueError('raw/model JSON mismatch')
    dest=folder/'decisions'/f'{unit["source_index"]:04d}.json'
    if dest.exists():raise ValueError('immutable decision already exists')
    dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(raw(a.decision_file))
    s['decisions'].append({'source_index':unit['source_index'],'prompt_id':unit['prompt_id'],'route_b_unit_id':unit['route_b_unit_id'],'file':str(dest),'sha256':digest(raw(dest)),'executor':a.executor,'genuine_authorship_attested':True,'authorship_note':a.authorship_note,'acknowledged_chunk_count_at_commit':len(s['acknowledged_chunks'])})
    s['unit_cursor']+=1;s['artifact_cursor']=0;s['char_cursor']=0;s['unit_ready']=False;save(statefile,s)
    output(dict(status(),committed_unit=unit['route_b_unit_id']));sys.exit()
if a.command=='export':
    if s['validation_only']:raise ValueError('mechanical validation cursor cannot export scientific response')
    if not a.out or (s['unit_cursor']!=len(s['units']) and not a.allow_partial):raise ValueError('all genuine decisions required before export unless --allow-partial')
    if a.out.exists():raise ValueError('immutable export already exists')
    model=s.get('model') or a.model
    if not model:raise ValueError('export requires the recorded requested model')
    responses=[]
    for d in s['decisions']:
        if digest(raw(d['file']))!=d['sha256']:raise ValueError('decision checkpoint changed')
        responses.append(read(d['file']))
    for unit in s['units'][s['unit_cursor']:]:
        reason='awaiting_independent_interpretation' if unit is s['units'][s['unit_cursor']] and s['unit_ready'] else 'reading_incomplete'
        responses.append({'prompt_id':unit['prompt_id'],'route_b_unit_id':unit['route_b_unit_id'],
                          'status':'not_evaluated','reason':reason})
    complete=s['unit_cursor']==len(s['units'])
    assigned={'assigned_reading_shard':s['reading_shard']} if s.get('reading_shard') else {'assigned_range':[s['start'],s['stop']]}
    save(a.out,{'case_id':s['case_id'],'executor':s['executor'],'model':model,'responses':responses,'reading_provenance':{'state_file':str(statefile),'state_sha256':digest(raw(statefile)),'source_sha256':s['source_sha256'],**assigned,'genuine_decision_count':len(s['decisions']),'all_chunks_acknowledged':complete,'completion_status':'complete' if complete else 'partial','checkpoint':s.get('checkpoint'),'no_automatic_scientific_decisions':True}})
    output({'state':'exported','responses':len(responses),'out':str(a.out),'sha256':digest(raw(a.out))})
