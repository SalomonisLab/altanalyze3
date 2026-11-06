"""Evaluate a small, capability-based plotting subset of Python, never exec/eval.

Generated programs receive tables, plotting functions and artist styling methods.
They cannot access files, network, environment, arbitrary imports or Python object
internals. The interpreter also runs in a separate process with bounded resources.
"""
from __future__ import annotations

import ast
from contextvars import ContextVar
import io
import inspect
import json
import operator
import os
import re
import shutil
import signal
from pathlib import Path
import subprocess
import sys
import tempfile
import textwrap
import warnings

_EVIDENCE = ContextVar('plot_evidence', default=None)


def bind_literals(value):
    """Reconcile verbatim numeric vectors only against a complete unique source.

    A model may copy a small table into its plotting program. Accept that only
    after exact ordered full-column comparison, never a preview/subset match.
    """
    context=_EVIDENCE.get()
    if isinstance(value,Data) or context is None:return value
    if isinstance(value,set):
        matches=[key for key,values in context['groups'].items()
                 if value==set(values) or value==set(context['aliases'].get(key,[]))]
        if matches and all(set(context['groups'][key])==set(context['groups'][matches[0]]) for key in matches):
            key=matches[0]
            context['bindings'].append(dict(group=key,rows=len(context['groups'][key]),verification='exact complete verified group identity/alias match'))
            return Data(set(context['groups'][key]))
    if isinstance(value,dict) and value and all(isinstance(v,(int,float)) or v is None for v in value.values()):
        matches=[(key,column) for key,table in context['tables'].items() for column,values in table.items()
                 if len(values)==len(value) and list(value.values())==values]
        if len(matches)==1:
            context['bindings'].append(dict(table=matches[0][0],column=matches[0][1],rows=len(value),
                keys=list(value),verification='exact ordered complete-column match for mapping values'))
            return Data(value)
    if isinstance(value,(list,tuple)) and value and all(isinstance(v,(int,float)) or v is None for v in value):
        matches=[(key,column) for key,table in context['tables'].items() for column,values in table.items()
                 if len(values)==len(value) and list(value)==values]
        if len(matches)==1:
            context['bindings'].append(dict(table=matches[0][0],column=matches[0][1],rows=len(value),verification='exact ordered complete-column match'))
            return Data(value)
    if isinstance(value,(list,tuple)):return type(value)(bind_literals(v) for v in value)
    if isinstance(value,dict):return {k:bind_literals(v) for k,v in value.items()}
    return value


class Data:
    def __init__(self, value): self.value = value


def raw(value):
    if isinstance(value, Data): return value.value
    if isinstance(value, API): return value.value
    if isinstance(value, (list, tuple)): return type(value)(raw(v) for v in value)
    if isinstance(value, dict): return {k: raw(v) for k, v in value.items()}
    return value


def grounded(value):
    if isinstance(value, Data): return True
    if isinstance(value, (list, tuple)): return any(grounded(v) for v in value)
    if isinstance(value, dict): return any(grounded(v) for v in value.values())
    return False


class Call:
    def __init__(self, function, data=False): self.function, self.data = function, data
    def __call__(self, args, kwargs):
        if self.data:
            args,kwargs=bind_literals(args),bind_literals(kwargs)
        if self.data and not grounded((args, kwargs)):
            context=_EVIDENCE.get() or {}
            columns=', '.join(key+'_'+re.sub(r'\W+','_',str(column)) for key,table in context.get('tables',{}).items() for column in table)
            raise ValueError('Plot values must come from tables, not invented literal measurements. '
                'Use the complete provided column lists ('+columns+') or groups[label] directly. '
                + self.usage())
        if self.data and _EVIDENCE.get() is not None:_EVIDENCE.get()['plots']+=1
        if any(k in kwargs for k in ('url', 'urls', 'usetex', 'data')):
            raise ValueError('External resources and alternate data sources are unavailable.')
        if 'figsize' in kwargs:
            size = raw(kwargs['figsize'])
            if len(size) != 2 or not all(1 <= float(v) <= 16 for v in size):
                raise ValueError('Choose a figure size between 1 and 16 inches.')
        if 'dpi' in kwargs and not 50 <= float(raw(kwargs['dpi'])) <= 150:
            raise ValueError('Use 50–150 dpi.')
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('error',UserWarning)
                value = self.function(*raw(args), **raw(kwargs))
        except (ValueError, TypeError, KeyError, AttributeError, UserWarning) as exc:
            raise ValueError(str(exc)+' Native API: '+self.usage()) from None
        return wrap(value, grounded((args, kwargs)))


    def usage(self):
        # Ground repairs in the installed API, without a question-specific recipe.
        try: signature=str(inspect.signature(self.function))
        except (TypeError,ValueError):signature=''
        return self.function.__name__+signature+'\n'+(inspect.getdoc(self.function) or '')[:2400]


class API:
    def __init__(self, value, names): self.value, self.names = value, names
    def attribute(self, name):
        if name.startswith('_') or name not in self.names:
            raise ValueError(f'Plotting attribute {name} is unavailable.')
        return Call(getattr(self.value, name), name in PLOTS)


PLOTS = {'plot', 'scatter', 'bar', 'barh', 'hist', 'boxplot', 'violinplot', 'imshow', 'pie', 'stackplot',
         'fill_between', 'errorbar', 'histplot', 'boxenplot', 'heatmap', 'stripplot', 'swarmplot',
         'venn2', 'venn3', 'venn2_circles', 'venn3_circles'}
STYLE = {'set_title', 'set_xlabel', 'set_ylabel', 'set_xlim', 'set_ylim', 'set_xticks', 'set_yticks',
         'set_xticklabels', 'set_yticklabels', 'set_xscale', 'set_yscale', 'legend', 'text', 'annotate', 'grid',
         'axis', 'axhline', 'axvline', 'suptitle', 'tight_layout', 'subplots_adjust', 'set_color', 'set_alpha',
         'set_fontsize', 'get_label_by_id', 'get_patch_by_id'}


def wrap(value, evidence=False):
    if value is None: return None
    if value.__class__.__module__.startswith(('matplotlib.', 'matplotlib_venn.')):
        return API(value, PLOTS | STYLE)
    if evidence: return Data(value)
    if isinstance(value, (list, tuple)): return type(value)(wrap(v) for v in value)
    return value


class Program:
    def __init__(self, tables, groups=None, aliases=None):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib_venn
        import seaborn
        self.plt = plt
        self.literal_bindings=[]
        groups,aliases=groups or {},aliases or {}
        _EVIDENCE.set({'tables':tables,'bindings':self.literal_bindings,'groups':groups,'aliases':aliases,'plots':0})
        plt.rcParams['text.usetex'] = False
        plt.rcParams['svg.fonttype'] = 'none'
        self.modules = {'matplotlib.pyplot': API(plt, PLOTS | {'figure', 'subplots', 'title', 'xlabel', 'ylabel', 'legend', 'text', 'annotate', 'axis', 'xlim', 'ylim', 'grid', 'tight_layout', 'show'}),
                        'matplotlib_venn': API(matplotlib_venn, {'venn2', 'venn3', 'venn2_circles', 'venn3_circles'}),
                        'seaborn': API(seaborn, PLOTS)}
        self.values = {'tables': Data(tables),'groups':Data({k:set(v) for k,v in groups.items()}),
            'plt':self.modules['matplotlib.pyplot'],'sns':self.modules['seaborn'],
            'matplotlib_venn':self.modules['matplotlib_venn']}
        for module in self.modules.values():
            for name in module.names:
                if name not in self.values and hasattr(module.value,name):
                    self.values[name] = module.attribute(name)
        for key, table in tables.items():
            self.values[key] = Data(table)
            for column, values in table.items():
                name = key + '_' + re.sub(r'\W+', '_', str(column))
                self.values[name] = Data(values)
        if tables:self.values['data'] = Data(next(reversed(tables.values())))
        for name, function in {'len': len, 'tuple': tuple, 'list': list, 'dict':dict, 'set':set, 'str': str, 'float': float, 'int': int, 'zip':zip,
                               'sum': sum, 'min': min, 'max': max, 'sorted': sorted}.items():
            self.values[name] = Call(function)
        self.steps = 0

    def value(self, node):
        self.steps += 1
        if self.steps > 5000: raise ValueError('Plot program is too large.')
        if isinstance(node, ast.Constant): return node.value
        if isinstance(node, ast.Name):
            if node.id.startswith('_') or node.id not in self.values: raise ValueError('Unknown plotting name: ' + node.id)
            return self.values[node.id]
        if isinstance(node, (ast.Tuple, ast.List)):
            values = [self.value(v) for v in node.elts]
            return tuple(values) if isinstance(node, ast.Tuple) else values
        if isinstance(node,ast.Set):return set(raw(self.value(v)) for v in node.elts)
        if isinstance(node, ast.Dict): return {raw(self.value(k)): self.value(v) for k, v in zip(node.keys, node.values)}
        if isinstance(node, ast.JoinedStr):
            return ''.join(str(raw(self.value(v.value))) if isinstance(v,ast.FormattedValue) else str(v.value) for v in node.values)
        if isinstance(node, ast.Subscript):
            value, key = self.value(node.value), raw(self.value(node.slice))
            if isinstance(value, Data): return Data(value.value[key])
            if not isinstance(value, (list, tuple, dict)): raise ValueError('Only table/list indexing is available.')
            return value[key]
        if isinstance(node, ast.Slice):
            return slice(*(raw(self.value(v)) if v is not None else None for v in (node.lower, node.upper, node.step)))
        if isinstance(node, ast.Attribute):
            value = self.value(node.value)
            if not isinstance(value, API): raise ValueError('Only declared plotting API attributes are available.')
            return value.attribute(node.attr)
        if isinstance(node, ast.Call):
            function = self.value(node.func)
            if not isinstance(function, Call): raise ValueError('Only declared plotting functions can be called.')
            if any(k.arg is None for k in node.keywords): raise ValueError('Expanded keyword arguments are unavailable.')
            return function([self.value(a) for a in node.args], {k.arg: self.value(k.value) for k in node.keywords})
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            value=self.value(node.operand); out=(1 if isinstance(node.op,ast.UAdd) else -1)*raw(value)
            return Data(out) if grounded(value) else out
        if isinstance(node, ast.BinOp) and type(node.op) in {ast.Add, ast.Sub, ast.Mult, ast.Div,ast.BitAnd,ast.BitOr,ast.BitXor}:
            left, right = self.value(node.left), self.value(node.right)
            if isinstance(node.op, ast.Mult):
                for a,b in ((raw(left),raw(right)),(raw(right),raw(left))):
                    if isinstance(a,(str,list,tuple)) and isinstance(b,int) and len(a)*abs(b)>2_000_000:
                        raise ValueError('Plot allocation is too large.')
            out={ast.Add:operator.add,ast.Sub:operator.sub,ast.Mult:operator.mul,ast.Div:operator.truediv,
                 ast.BitAnd:operator.and_,ast.BitOr:operator.or_,ast.BitXor:operator.xor}[type(node.op)](raw(left),raw(right))
            return Data(out) if grounded((left,right)) else out
        raise ValueError('Unsupported plotting expression: ' + type(node).__name__)

    def assign(self, target, value):
        if isinstance(target, ast.Name) and not target.id.startswith('_') and target.id not in self.modules:
            self.values[target.id]=bind_literals(value)
        elif isinstance(target,(ast.Tuple,ast.List)):
            values=raw(value) if isinstance(value,Data) else value
            if len(values)!=len(target.elts):raise ValueError('Assignment dimensions do not match.')
            for key,item in zip(target.elts,values):self.assign(key,Data(item) if isinstance(value,Data) else item)
        else:raise ValueError('Only local names can be assigned.')

    def run(self, code, caption):
        if len(code)>12000:raise ValueError('Plot program exceeds 12 KB.')
        tree=ast.parse(code)
        if sum(1 for _ in ast.walk(tree))>1000:raise ValueError('Plot program is too large.')
        for statement in tree.body:
            if isinstance(statement,ast.Import):
                for alias in statement.names:
                    if alias.name not in self.modules:raise ValueError('Import is unavailable: '+alias.name)
                    self.values[alias.asname or alias.name]=self.modules[alias.name]
            elif isinstance(statement,ast.ImportFrom):
                if statement.module not in self.modules or statement.level:raise ValueError('Import is unavailable.')
                for alias in statement.names:self.values[alias.asname or alias.name]=self.modules[statement.module].attribute(alias.name)
            elif isinstance(statement,ast.Assign):
                value=self.value(statement.value)
                for target in statement.targets:self.assign(target,value)
            elif isinstance(statement,ast.Expr):self.value(statement.value)
            else:raise ValueError('Only imports, assignments and plotting calls are available.')
        if not self.plt.get_fignums():raise ValueError('The program did not create a figure.')
        if not _EVIDENCE.get()['plots']:raise ValueError('The program did not plot any retrieved data.')
        figure=self.plt.gcf()
        width,height=figure.get_size_inches()
        if not 1<=width<=16 or not 1<=height<=16:raise ValueError('Figure exceeds size bounds.')
        figure.set_dpi(120)
        figure.subplots_adjust(bottom=.24)
        figure.text(.025,.025,'\n'.join(textwrap.wrap(caption,100)),fontsize=9,va='bottom',wrap=True)
        buffer=io.StringIO()
        # A bounded full canvas avoids the tight-bbox allocation failure.
        figure.savefig(buffer,format='svg',dpi=120)
        self.plt.close('all')
        if len(buffer.getvalue())>4*1024**2:raise ValueError('Vector plot exceeds 4 MB; use the established interactive view or an explicit aggregate. No rows were sampled.')
        return buffer.getvalue()


def render(code, tables, caption, with_provenance=False,groups=None,aliases=None):
    payload=json.dumps(dict(code=code,tables=tables,caption=caption,groups=groups,aliases=aliases),default=str)
    if len(payload)>32*1024**2:raise ValueError('Selected plotting data exceed 32 MB; use an established view or explicit aggregate.')
    with tempfile.TemporaryDirectory(prefix='scalable-plot-') as directory:
        # An isolated HOME must not trigger an OS-wide font inventory on every
        # request (system_profiler on macOS can exceed the render timeout).
        # Seed only native font JSON, never the host's matplotlibrc/configuration.
        import matplotlib
        for font_cache in Path(matplotlib.get_cachedir()).glob('fontlist*.json'):
            try: shutil.copyfile(font_cache,Path(directory)/font_cache.name)
            except OSError: pass  # A concurrently replaced cache can be rebuilt.
        env={'PATH':os.environ.get('PATH',''),'HOME':directory,'MPLCONFIGDIR':directory,
             'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
        process=subprocess.Popen([sys.executable,str(Path(__file__).resolve()),'--worker'],text=True,
            stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd=directory,env=env,
            start_new_session=True)
        try:
            stdout,stderr=process.communicate(payload,timeout=30)
        except subprocess.TimeoutExpired:
            # Font discovery/native plotting may spawn children; kill this
            # private process group rather than leaving descendants running.
            if os.name=='posix':
                try:os.killpg(process.pid,signal.SIGKILL)
                except ProcessLookupError:pass
            else:process.kill()
            process.communicate()
            raise ValueError('Plotting exceeded its 30-second time limit.') from None
    if process.returncode:raise ValueError('Plotting failed: '+stderr[-2000:])
    result=json.loads(stdout)
    if result.get('error'):raise ValueError(result['error'])
    return result if with_provenance else result['svg']


if __name__=='__main__':
    try:
        import resource
        resource.setrlimit(resource.RLIMIT_CPU,(20,20))
        if sys.platform.startswith('linux'):resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
        request=json.loads(sys.stdin.read(32*1024**2+1))
        program=Program(request['tables'],request.get('groups'),request.get('aliases'))
        from importlib.metadata import version
        print(json.dumps({'svg':program.run(request['code'],request['caption']),'bindings':program.literal_bindings,
            'versions':{name:version(name)for name in ('matplotlib','matplotlib-venn','seaborn')}}))
    except Exception as exc:
        print(json.dumps({'error':f'{type(exc).__name__}: {exc}'}))
