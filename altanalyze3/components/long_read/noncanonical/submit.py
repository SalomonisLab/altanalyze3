"""Submit the same complete pipeline as source-local LSF jobs and checked barriers."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from .manifest import load_config,atomic_json,digest


def submit(config_path,queue='normal',memory=32000,walltime='8:00',dry_run=False):
    config=load_config(config_path)
    root=Path(config['output']);logs=root/'logs';logs.mkdir(parents=True,exist_ok=True)
    namespace='aa3nc_'+config['_sha256'][:12]
    records=[]
    def job(stage,name,dependencies=None,source=None):
        command=['bsub','-q',queue,'-n','1','-R',f'rusage[mem={memory}] span[hosts=1]',
                 '-M',str(memory),'-W',walltime,'-J',name,'-o',str(logs/(name+'.out')),'-e',str(logs/(name+'.err'))]
        if dependencies:command+=['-w',' && '.join(f'done({n})' for n in dependencies)]
        command += [sys.executable,'-m','altanalyze3.components.long_read.noncanonical.pipeline',
                    '--config',str(Path(config_path).resolve()),'--stage',stage]
        if source:command+=['--source',source]
        record=dict(stage=stage,source_id=source,name=name,command=command)
        if not dry_run:
            env=dict(os.environ,PYTHONHASHSEED='0',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',
                     ISOFORM_MAX_PROCESSES='1',ALTANALYZE3_PSI_JOBS='1')
            env['PYTHONPATH']=str(Path(__file__).resolve().parents[4])
            response=subprocess.check_output(command,env=env,text=True)
            match=re.search(r'Job <(\d+)>',response)
            if not match:raise RuntimeError('Could not record submitted job: '+response)
            record['job_id']=match[1]
        records.append(record)
        return record.get('job_id',name)
    # Source identities, not biological UIDs/library labels, define independent jobs.
    jobs=[job('extract',namespace+'_'+s['source_id'],source=s['source_id']) for s in config['sources']]
    cohort=job('aggregate',namespace+'_aggregate',jobs)
    final=job('standard',namespace+'_standard',[cohort])
    atomic_json(root/('submission_plan.json' if dry_run else 'SUBMITTED.json'),
                dict(config_sha256=digest(config_path),jobs=records,final_job=final,status='PLANNED' if dry_run else 'SUBMITTED_NOT_COMPLETE'))
    return records


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',required=True)
    p.add_argument('--queue',default='normal');p.add_argument('--memory',type=int,default=32000)
    p.add_argument('--walltime',default='8:00');p.add_argument('--dry-run',action='store_true')
    a=p.parse_args();print(json.dumps(submit(a.config,a.queue,a.memory,a.walltime,a.dry_run),indent=2))

if __name__=='__main__':main()
