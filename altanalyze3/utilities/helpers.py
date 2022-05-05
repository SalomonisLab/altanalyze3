import time
import random
import string
import logging
import pkg_resources


def get_version():
    """
    Returns current version of the package if it's installed
    """
    pkg = pkg_resources.require("altanalyze3")
    return pkg[0].version if pkg else "unknown version"


def get_tmp_suffix(length=None):
    """
    Returns random filename extension as a string of specified length
    """
    length = 10 if length is None else length
    return "." + "".join(random.choices(string.ascii_uppercase + string.digits, k=length))


class TimeIt():
    """
    Prints elapsed time to console log (info level)
    """

    def __enter__(self):
        self.start = time.time()

    def __exit__(self, type, value, traceback):
        logging.info (f"""Elapsed time: {round(time.time() - self.start)} sec""")
