#!/usr/bin/env python3

try:
    from .core import *
except:
    from core import *

import sys
import click

import time

click_File_settings = dict(
    exists=True,
    file_okay=True, 
    dir_okay=False,
)
click_Path_settings = dict(
    exists=True,
    file_okay=False, 
    dir_okay=True,
    writable=True,
)
@click.command(
    context_settings=dict(
        ignore_unknown_options=True,
        help_option_names=['-h', '--help'],
    ),
)
@click.option('-n', '--nfo', prompt='Sample Info. file', help='Sample Info. xlsx file', required=True, type=click.Path(click_File_settings))
@click.option('-m', '--ms', prompt='MS file', help='MS result xlsx file', required=True, type=click.Path(click_File_settings))
@click.option('-l', '--hla', prompt='HLA file', help='HLA result xlsx file', required=True, type=click.Path(click_File_settings))
@click.argument('meds', metavar='<MED.xlsx files>', nargs=-1, required=True, type=click.Path(click_File_settings))
@click.option('-o', '--output', prompt='Output path', help='Output path', required=True, type=click.Path(click_Path_settings))
@click.option('-g', '--gui', is_flag=True, help='Force GUI window', default=False)
@click.option('-v', '--verbose', is_flag=True, help='Enables verbose mode', default=False)
def tcli(nfo,ms,hla,meds,output,gui,verbose):
    print(verbose)
    mainfunc(nfofile=nfo,msfile=ms,hlafile=hla,medfiles=meds,outpath=output)


def main():
    if len(sys.argv) > 1:
        click.echo('Hello World !')
        click.echo('[E]Hello World!', err=True)
        click.secho('Hello World!', fg='green')
        click.secho('Some more text', bg='blue', fg='white')
        click.secho('ATTENTION', blink=True, bold=True)
        tcli()

if __name__ == '__main__':
    main()
