import os,sys
from altanalyze3.utilities.helpers import get_version

from altanalyze3.utilities.parser import ArgsParser


def run(command: str, **kwargs):
    """
    Programmatically run AltAnalyze3 commands.

    Parameters:
        command (str): One of "juncount", "intcount", "aggregate", "index".
        kwargs (dict): All CLI arguments as keyword args.
    """
    args_list = [command]
    for k, v in kwargs.items():
        if isinstance(v, bool):
            if v:
                args_list.append(f"--{k}")
        elif isinstance(v, list):
            args_list.append(f"--{k}")
            args_list.extend(map(str, v))
        else:
            args_list.append(f"--{k}")
            args_list.append(str(v))
    
    parser = ArgsParser(args_list)
    parser.args.func(parser.args)
