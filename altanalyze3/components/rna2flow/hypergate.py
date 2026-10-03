"""Optional bridge to the official HyperGate R package; no substitute implementation."""
import csv
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import numpy as np


def hypergate_conditions(X, target, features, beta=1., timeout=120):
    executable=os.environ.get('RNA2FLOW_RSCRIPT') or shutil.which('Rscript')
    if not executable:
        candidate=Path('/Library/Frameworks/R.framework/Resources/bin/Rscript')
        if candidate.exists():executable=str(candidate)
    if not executable:raise ValueError('HyperGate requires Rscript and the official hypergate R package')
    with tempfile.TemporaryDirectory(prefix='rna2flow_hypergate_') as tmp:
        root=Path(tmp)
        # Safe numeric column names avoid collisions and quoting of arbitrary antibody names.
        with (root/'matrix.tsv').open('w') as handle:
            handle.write('\t'.join('m'+str(i) for i in range(len(features)))+'\n')
            np.savetxt(handle,X,delimiter='\t',fmt='%.17g')
        np.savetxt(root/'target.tsv',np.asarray(target,int),fmt='%d')
        script='''args <- commandArgs(trailingOnly=TRUE)
local_lib <- path.expand('~/.local/lib/R')
if(dir.exists(local_lib)) .libPaths(c(local_lib,.libPaths()))
suppressPackageStartupMessages(library(hypergate))
x <- as.matrix(read.delim(args[1],check.names=FALSE))
t <- scan(args[2],quiet=TRUE)
g <- hypergate(xp=x,gate_vector=t,level=1,beta=as.numeric(args[4]),verbose=FALSE)
write.table(hgate_info(g),args[3],sep='\\t',row.names=FALSE,quote=FALSE)
'''
        (root/'run.R').write_text(script)
        try:
            result=subprocess.run([executable,str(root/'run.R'),str(root/'matrix.tsv'),str(root/'target.tsv'),str(root/'rules.tsv'),str(beta)],capture_output=True,text=True,timeout=timeout)
        except subprocess.TimeoutExpired as e:raise ValueError('Official HyperGate exceeded the search time limit') from e
        if result.returncode:raise ValueError('Official HyperGate failed: '+result.stderr[-800:])
        conditions=[]
        with (root/'rules.tsv').open() as handle:
            for row in csv.DictReader(handle,delimiter='\t'):
                j=int(row['channels'][1:]);threshold=float(row['threshold']);op='>' if row['sign']=='+' else '<='
                # Official HyperGate uses >= for positive gates. Represent it exactly
                # in this workflow's > convention for the original floating dtype.
                if op=='>':threshold=float(np.nextafter(np.asarray(threshold,dtype=X.dtype),np.asarray(-np.inf,dtype=X.dtype)))
                conditions.append(dict(channel=features[j],operator=op,threshold=threshold))
        return conditions
