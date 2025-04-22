"""
    author: Suhas Vittal
"""

import re

retype1 = r'--([0-9A-Za-z_-]+)'
retype2 = r'-([0-9A-Za-z_-]+)'

def parse(args: list[str]) -> dict[str, str]:
    arg_list = {}
    ptr = 0
    while ptr < len(args):
        s = args[ptr]
        m = re.match(retype1, s)
        if m:
            identifier = m.group(1)
            data = args[ptr+1]
            arg_list[identifier] = data
            ptr += 2
            continue
        m = re.match(retype2, s)
        if m:
            identifier = m.group(1)
            arg_list[identifier] = None
        ptr += 1
    return arg_list

def try_get_string(arg_list: dict[str, str], arg: str) -> str:
    if arg not in arg_list:
        exit(1)
    return arg_list[arg]

def try_get_int(arg_list: dict[str, str], arg: str) -> int:
    return int(try_get_string(arg_list, arg))

def try_get_float(arg_list: dict[str, str], arg: str) -> float:
    return float(try_get_string(arg_list, arg))
