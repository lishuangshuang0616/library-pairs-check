import gzip,json,re,subprocess,time,os,tempfile,argparse
from pathlib import Path
parser=argparse.ArgumentParser(description="Generate synthetic R1 data in sample directories and measure the release binary on macOS. Uses several hundred MB of temporary disk space.")
parser.add_argument('--output-dir',type=Path)
options=parser.parse_args()
if options.output_dir:
    root=options.output_dir.resolve();root.mkdir(parents=True,exist_ok=False)
else:
    root=Path(tempfile.mkdtemp(prefix='library-pair-benchmark-'))
print('Benchmark directory:',root,flush=True)
project=Path(__file__).resolve().parents[1]
white=(project/'assets/barcode_whitelist.txt').read_text().splitlines()
binary=str(project/'target/release/library-pair-check')
rows=[]
def write_file(path,barcodes,copies,kind,compressed):
    opener=gzip.open if compressed else open
    with opener(path,'wb') as out:
        chunk=[]
        for i in range(barcodes):
            left,right=white[i//len(white)],white[i%len(white)]
            seq=(f'TCTGCG{left}CCTTCC{right}CGATGACGTACGTACTTTTTTTTTTT' if kind=='cdna' else f'{left}{right}')
            rec=f'@r{i}\n{seq}\n+\n'+('I'*len(seq))+'\n'
            chunk.append(rec.encode()*copies)
            if len(chunk)==2000:
                out.write(b''.join(chunk));chunk=[]
        out.write(b''.join(chunk))
def measure(name,cdna,oligo,records):
    base=root/name
    for role,paths in [('cDNA',cdna),('oligo',oligo)]:
        target=base/'SampleA'/role
        target.mkdir(parents=True,exist_ok=True)
        for idx,path in enumerate(paths):
            suffix='.fq.gz' if path.name.endswith('.gz') else '.fq'
            os.link(path,target/f'L{idx+1:03d}_R1_001{suffix}')
    output=root/(name+'.json')
    args=['/usr/bin/time','-l',binary,'--input',str(base),'--threads','1','--output',str(output)]
    t=time.perf_counter();r=subprocess.run(args,capture_output=True,text=True);wall=time.perf_counter()-t
    (root/(name+'.stderr')).write_text(r.stderr)
    if r.returncode:raise RuntimeError(r.stderr)
    report=json.loads(output.read_text())
    rss=int(re.search(r'(\d+)\s+maximum resident set size',r.stderr).group(1))
    item={'scenario':name,'records':records,'wall_seconds':round(wall,3),'million_records_per_second':round(records/wall/1e6,3),'peak_rss_mib':round(rss/1024**2,2),'supported_barcodes_per_library':report['cdna'][0]['supported_barcodes'],'observed_barcodes_per_library':report['cdna'][0]['observed_barcodes'],'pairs':len(report['pairs'])}
    rows.append(item);(root/'results.json').write_text(json.dumps(rows,indent=2));print(json.dumps(item),flush=True)
for name,barcodes,copies,compressed in [('small_plain',20000,5,False),('medium_plain',200000,5,False),('medium_gzip',200000,5,True),('large_gzip',600000,3,True),('singletons_gzip',1000000,1,True)]:
    paths=[]
    for kind in ['cdna','oligo']:
        p=root/(name+'_'+kind+('.fq.gz' if compressed else '.fq'))
        write_file(p,barcodes,copies,kind,compressed);paths.append(p)
    measure(name,[paths[0]],[paths[1]],2*barcodes*copies)
# Multiple lanes: several R1 files inside each library directory.
paths=[]
for kind in ['cdna','oligo']:
    group=[]
    source=root/('medium_gzip_'+kind+'.fq.gz')
    for i in range(4):
        p=root/f'multi_{kind}_{i}.fq.gz';os.link(source,p);group.append(p)
    paths.append(group)
measure('four_by_four_gzip',paths[0],paths[1],8000000)
