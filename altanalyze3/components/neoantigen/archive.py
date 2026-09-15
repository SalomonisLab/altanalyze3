"""Extract data-only reference archives without links or traversal."""
import shutil
import tarfile
from pathlib import Path


def unpack_reference(archive, outdir):
    root = Path(outdir).resolve(); root.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive) as tf:
        members = tf.getmembers()
        for member in members:
            dest = (root / member.name).resolve()
            if not dest.is_relative_to(root) or not (member.isfile() or member.isdir()):
                raise ValueError(f'Unsafe reference archive member: {member.name}')
        for member in members:
            dest = root / member.name
            if member.isdir():
                dest.mkdir(parents=True, exist_ok=True)
            else:
                dest.parent.mkdir(parents=True, exist_ok=True)
                with tf.extractfile(member) as inp, dest.open('wb') as out:
                    shutil.copyfileobj(inp, out)
    if not (root/'Alt91_db').is_dir() or not (root/'controls').is_dir():
        raise ValueError('Reference archive must contain Alt91_db and controls at its root')
