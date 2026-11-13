"""CLI hooks; all modalities are requested unless explicitly disabled in config."""
from pathlib import Path
from .manifest import load_config
from .pipeline import main as run_pipeline

def run(args):
    argv=['--config',args.config,'--stage',args.stage]
    if getattr(args,'source',None):argv+=['--source',args.source]
    run_pipeline(argv)

def add_arguments(parser):
    parser.add_argument('--config',required=True,help='Pinned source/reference/baseline-method JSON manifest')
    parser.add_argument('--stage',choices=['extract','aggregate','standard','all'],default='all')
    parser.add_argument('--source',default=None,help='Exact source ID for per-source extract jobs')

def phase1(args):
    """Called by the existing sclr entry point before declaring phase 1 complete."""
    if getattr(args,'legacy_only',False):return
    candidate=getattr(args,'noncanonical_config',None) or Path(args.metadata).with_name('noncanonical.json')
    if not Path(candidate).is_file():
        if getattr(args,'noncanonical_config',None):
            raise FileNotFoundError(candidate)
        # sclr is one legacy phase, not the complete workflow. Its established
        # inputs remain sufficient. lr-integrate / the full launcher default to all.
        return
    config=load_config(candidate)
    from .pipeline import extract_source
    sources=[s for s in config['sources'] if not getattr(args,'sample',None) or s['uid']==args.sample]
    if not sources:raise ValueError('No source in modality manifest matches --sample')
    for source in sources:extract_source(config,source,standard=False)
