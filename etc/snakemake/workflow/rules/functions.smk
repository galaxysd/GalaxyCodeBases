import sys
if sys.version_info < (3, 10):
    raise ValueError("[x]This workflow requires at least Python 3.10.")

import os
import pprint
pp = pprint.PrettyPrinter(indent=4)

def eprint(*args, **kwargs) -> None:
    print(*args, **kwargs, file=sys.stderr, flush=True)

def get_all_sha2(wildcards):
    targets = list()
    for cz in FileList:
        fn = 'results/'+os.path.splitext(cz)[0]+'.sha2'
        targets.append(fn)
    return targets

from prettytable import PrettyTable
rule help:
    '''
    Print help menu.
    <https://www.biostars.org/p/220268/#238766>
    '''
    localrule: True
    run:
        tbl = PrettyTable(['Rules', 'Description'])
        tbl.align = 'l'
        for rule in workflow.rules:
            if not rule.name.startswith('_'):
                if rule.docstring:
                    tbl.add_row([rule.name, rule.docstring.splitlines()[0]])
                    #tbl.add_row([rule.name, rule.docstring])
        print(tbl)
