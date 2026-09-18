"""Check schemas, wrapper definitions and required data inside the built wheel."""
import argparse
import json
import re
import zipfile
from pathlib import Path
import xml.etree.ElementTree as ET
from jsonschema import Draft7Validator

root=Path(__file__).resolve().parents[2]
schema=json.loads((root/'components/snaf/nextflow/nextflow_schema.json').read_text())
Draft7Validator.check_schema(schema)
validator=Draft7Validator(schema)
assert not list(validator.iter_errors({'input':'samples.csv','mode':'bam','min_reads':20}))
assert list(validator.iter_errors({'input':'samples.csv','min_reads':-1}))
assert list(validator.iter_errors({'input':'samples.csv','with_pyneoquant':'yes'}))
config=(root/'components/snaf/nextflow/nextflow.config').read_text().split('params {',1)[1].split('\n}',1)[0]
params=set(re.findall(r'^\s*(\w+)\s*=',config,re.M))
assert params==set(schema['properties']), (params-set(schema['properties']),set(schema['properties'])-params)
for xml in (root/'deployment/galaxy').glob('*.xml'):
    tool=ET.parse(xml).getroot()
    assert tool.find('tests/test') is not None
    for test in tool.findall('tests/test'):
        for param in test.findall('param'):
            declared=tool.find(f"inputs/param[@name='{param.attrib['name']}']")
            if declared is not None and declared.get('type')=='data':
                assert (xml.parent/'test-data'/param.attrib['value']).is_file()
for workflow_path in (root/'deployment/galaxy/workflows').glob('*.ga'):
    workflow=json.loads(workflow_path.read_text())
    steps=workflow['steps']
    for step in steps.values():
        if step.get('type')!='tool':
            continue
        tool=next(ET.parse(xml).getroot() for xml in (root/'deployment/galaxy').glob('*.xml')
                  if ET.parse(xml).getroot().get('id')==step['tool_id'])
        assert tool.get('version')==step['tool_version']
        assert {o['output_name'] for o in step['workflow_outputs']}=={o.get('name') for o in tool.find('outputs')}
        for name,connection in step['input_connections'].items():
            assert tool.find(f"inputs/param[@name='{name}']") is not None
            assert str(connection['id']) in steps
p=argparse.ArgumentParser();p.add_argument('--wheel');args=p.parse_args()
if args.wheel:
    with zipfile.ZipFile(args.wheel) as z:
        names=set(z.namelist())
        for suffix in ['altanalyze3/components/neoantigen/cli.py',
                       'altanalyze3/components/neoantigen/data/ipepgen_contract.json',
                       'altanalyze3/components/bam/bam2hla/build/signatures_hg38.json.gz',
                       'altanalyze3/components/bam/bam2hla/build/signatures_hg19.json.gz',
                       'altanalyze3/components/snaf/BayesTS/GTEx_BayesTS.tsv.gz',
                       'altanalyze3/components/snaf/deepimmuno/models/cnn_model_331_3_7/variables.index',
                       'altanalyze3/components/snaf/deepimmuno/data/after_pca.txt']:
            assert suffix in names,suffix
        assert not any(n.startswith('pyneoquant/') for n in names)
print('Schemas, Galaxy fixtures and packaged resources verified')
