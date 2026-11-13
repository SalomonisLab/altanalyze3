"""Shim. The rule module lives in the fork:
forked_altanalyze3/altanalyze3/components/bam/cross_chromosome_rules.py (loaded by file path so the
package __init__ files are not imported)."""
import importlib.util
from pathlib import Path

_path = Path(__file__).resolve().parents[1] / 'forked_altanalyze3/altanalyze3/components/bam/cross_chromosome_rules.py'
_spec = importlib.util.spec_from_file_location('cross_chromosome_rules', _path)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
globals().update({k: v for k, v in vars(_mod).items() if not k.startswith('__')})
