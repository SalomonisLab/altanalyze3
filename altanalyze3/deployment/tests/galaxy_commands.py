"""Execute commands rendered from Galaxy XML test cases, without a Galaxy server.

This checks the real Cheetah templates and output assertions. It is not a Galaxy
job-runner or Tool Shed integration test. Use planemo test for server acceptance.
"""
import argparse
import os
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path
from Cheetah.Template import Template

p=argparse.ArgumentParser();p.add_argument('tools',nargs='+');a=p.parse_args()
for name in a.tools:
    path=Path(name).resolve();tool=ET.parse(path).getroot()
    for i,test in enumerate(tool.findall('tests/test')):
        values={}
        for param in tool.findall('inputs/param'):
            key=param.get('name')
            if param.get('type')=='select':values[key]=param.find('option').get('value')
            else:values[key]=param.get('value','')
        for param in test.findall('param'):
            key=param.get('name');value=param.get('value')
            spec=tool.find(f"inputs/param[@name='{key}']")
            if spec.get('type')=='data':value=str(path.parent/'test-data'/value)
            values[key]=value
        command=str(Template(tool.find('command').text,searchList=[values]))
        with tempfile.TemporaryDirectory(prefix='snaf-galaxy-') as td:
            result=subprocess.run(command,shell=True,executable='/bin/bash',cwd=td,text=True,capture_output=True)
            failure=test.get('expect_failure')=='true'
            assert bool(result.returncode)==failure,(path.name,result.stdout,result.stderr)
            for assertion in test.findall('assert_stderr/has_text'):
                assert assertion.get('text') in result.stderr,result.stderr
            for expected in test.findall('output'):
                out=tool.find(f"outputs/data[@name='{expected.get('name')}']")
                actual=Path(td,out.get('from_work_dir')).read_text()
                if expected.get('file'):
                    assert actual==(path.parent/'test-data'/expected.get('file')).read_text()
                for assertion in expected.findall('assert_contents/has_text'):
                    assert assertion.get('text') in actual,(path,actual)
        print(f'{path.name}: test {i+1} passed')
