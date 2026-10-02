"""Download lipid mzMLs with a MassIVE browser session and verified byte counts."""
from __future__ import annotations

import argparse
import hashlib
from http.cookiejar import CookieJar
import json
from pathlib import Path
import re
import time
from urllib.parse import urlencode
from urllib.request import build_opener, HTTPCookieProcessor

ACCESSION='MSV000081973'
TASK='1c494cc549464750b8b13358241dae21'
BASE='https://massive.ucsd.edu/ProteoSAFe/'


def download(out):
    out=Path(out);directory=out/'mzML';directory.mkdir(parents=True,exist_ok=True)
    opener=build_opener(HTTPCookieProcessor(CookieJar()))
    opener.addheaders=[('User-Agent','Mozilla/5.0'),('Referer',BASE+'dataset_files.jsp?task='+TASK)]
    with opener.open(BASE+'dataset_files.jsp?task='+TASK,timeout=45) as response:
        html=response.read().decode()
    found=re.search(r'var dataset_files = (.*?);',html)
    if not found:raise ValueError('MassIVE file inventory unavailable')
    inventory=json.loads(found[1]);(out/'dataset_inventory.json').write_text(json.dumps(inventory,indent=2)+'\n')
    rows=[r for r in inventory['row_data'] if r['collection']=='peak' and '_L_' in r['name'] and r['name'].endswith('.mzML')]
    if len(rows)!=30:raise ValueError(f'Expected 30 original lipid runs, got {len(rows)}')
    records=[]
    for r in rows:
        path=directory/r['name'];url=BASE+'DownloadResultFile?'+urlencode({'file':r['file_descriptor'],'forceDownload':'true'})
        if not path.exists() or path.stat().st_size!=r['size']:
            temporary=path.with_suffix('.download')
            for attempt in range(3):
                try:
                    with opener.open(url,timeout=180) as response,temporary.open('wb') as handle:
                        while chunk:=response.read(1024*1024):handle.write(chunk)
                    if temporary.stat().st_size!=r['size']:raise ValueError('Byte count mismatch for '+r['name'])
                    temporary.rename(path);break
                except Exception:
                    if attempt==2:raise
                    time.sleep(3)
        with path.open('rb') as handle:
            prefix=handle.read(1000)
        if b'mzML' not in prefix or b'<html' in prefix.lower():raise ValueError('Not an mzML file: '+str(path))
        records.append({'name':r['name'],'file_descriptor':r['file_descriptor'],'bytes':path.stat().st_size,
                        'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'url':url})
        (out/'download_manifest.json').write_text(json.dumps({'accession':ACCESSION,'task':TASK,'files':records},indent=2)+'\n')
        print(len(records),r['name'],flush=True)
    return records


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',required=True)
    download(parser.parse_args().out)
