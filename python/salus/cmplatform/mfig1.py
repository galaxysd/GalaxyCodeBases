#!/usr/bin/env python3
import sys
import os
from typing import NamedTuple

PlatformTuple = ('Visium', 'Salus')
"""
赛陆鼠脑数据：/share/result/spatial/CYL10/YL1025/E1new/results/ 配对HE图：（）  bin400
10X鼠脑数据：/share/result/spatial/data/BoAo_sp/sub2/srp_Illumina_mbrain/outs 
赛陆鼠肾数据：/share/result/spatial/CYL10/YL1025/D1new/results/  bin400
10X鼠肾数据：/share/result/spatial/data/BoAo_sp/sub2/srp_Illumina_mkidney/outs
赛陆睾丸数据：/share/result/spatial/HY0606_mouse/HY0606B2_supplementary/results/  配对HE图   bin40,not 400.
"""
SamplesDict = {
    'mbrainS': {
        'type': 'Salus',
        'prefix': '/share/result/spatial/CYL10/YL1025/E1new/results'
        'bin': 400,
    },
    'mbrainV': {
        'type': 'Visium',
        'prefix': '/share/result/spatial/data/BoAo_sp/sub2/srp_Illumina_mbrain/outs'
        'bin': 400,
    },
    'mkidneyS': {
        'type': 'Salus',
        'prefix': '/share/result/spatial/CYL10/YL1025/D1new/results'
        'bin': 400,
    },
    'mkidneyV': {
        'type': 'Visium',
        'prefix': '/share/result/spatial/data/BoAo_sp/sub2/srp_Illumina_mkidney/outs'
        'bin': 400,
    },
    'mtesticleS': {
        'type': 'Salus',
        'prefix': '/share/result/spatial/HY0606_mouse/HY0606B2_supplementary/results'
        'bin': 40,
    },
}
