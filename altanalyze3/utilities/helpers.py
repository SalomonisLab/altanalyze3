import time
import random
import string
import logging
import hashlib
import pkg_resources


def get_version():
    """
    Returns current version of the package if it's installed.

    Must NOT hard-fail on dependency mismatches. ``pkg_resources.require()`` resolves and enforces
    EVERY pinned dependency in the package metadata, so a single version skew (e.g. an installed
    tqdm newer than the pinned ``tqdm==4.62.3``) raises ``VersionConflict`` and crashes the whole CLI
    at startup -- even though the version string itself is trivially available. Look up just the
    distribution version, and fall back gracefully if anything goes wrong.
    """
    try:
        return pkg_resources.get_distribution("altanalyze3").version
    except Exception:
        try:
            # Last resort: importlib.metadata (py3.8+), independent of pkg_resources.
            from importlib.metadata import version as _il_version
            return _il_version("altanalyze3")
        except Exception:
            return "unknown version"


def get_tmp_suffix(length=None):
    """
    Returns random filename extension as a string of specified length
    """
    length = 10 if length is None else length
    return "." + "".join(random.choices(string.ascii_uppercase + string.digits, k=length))


def get_md5_sum(location, block_size=2**20):
    md5_sum = hashlib.md5()
    with open(location , "rb") as input_stream:
        while True:
            buf = input_stream.read(block_size)
            if not buf:
                break
            md5_sum.update(buf)
    return md5_sum.hexdigest()

class TimeIt():
    """
    Prints elapsed time to console log (info level)
    """

    def __enter__(self):
        self.start = time.time()

    def __exit__(self, type, value, traceback):
        logging.info (f"""Elapsed time: {round(time.time() - self.start)} sec""")
