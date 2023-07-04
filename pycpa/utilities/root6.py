import sys
import pathlib
import subprocess
import numpy as np
import pandas as pd

import string
import random
from typing import Literal

try:
    import ROOT
except:
    path_root = subprocess.check_output('which root', shell=True)
    path_root = str(path_root, encoding='utf-8').strip()
    path_root = pathlib.Path(str(path_root)).parent.parent
    path_root = pathlib.Path(path_root, 'lib')
    sys.path.append(str(path_root))
    import ROOT

class RandomNameGenerator:
    def __init__(self):
        self.used_name = set()

    def generate(self, 
        len:int=10,
        pattern:Literal['char','digit','mixed'] = 'mixed'
    ) -> str:
        """generate a random string of length 10.
        Parameter
        ---------
        pattern:
            - `mixed` : string starts with character followed by 9 char/digit. some
            - `char` : character only
            - `digit` : digits only
        len : int
            length of the sequence
        """
        while True:
            key = self._generate(pattern=pattern, len=len)
            if not key in self.used_name:
                self.used_name.add(key)
                break
        return key

    def _generate(self, pattern, len):
        if pattern == 'char':
            entries = [random.choice(string.ascii_letters) for _ in range(len)]
        elif pattern == 'digit':
            entries = [random.choice(string.digits) for _ in range(len)]
        elif pattern == 'mixed':
            first_entry = random.choice(string.ascii_letters)
            entries = [first_entry] + [random.choice(string.ascii_letters + string.digits) for _ in range(1, len)]
        return ''.join(entries)

rdnkey_generator = RandomNameGenerator()

class TFile:
    def __init__(self, path, mode : Literal['READ', 'RECREATE']='READ'):
        self.path = pathlib.Path(path)
        self.mode = mode.upper()

    def __enter__(self):
        self.file = ROOT.TFile(str(self.path), self.mode)
        return self

    def __getitem__(self, key):
        rdnkey = rdnkey_generator.generate()
        obj = self.file.Get(key)
        obj.SetName(rdnkey)
        return obj

    def __exit__(self, *args):
        self.file.Close()

    def keys(self):
        return [key.GetName() for key in self.file.GetListOfKeys()]

class HistogramReader:
    def __init__(self):
        pass
    
    def __call__(self, hist, *args, **kwargs):
        return self.histogram_to_df(hist, *args, **kwargs)

    def histogram_to_df(self, hist, *args, **kwargs):
        acceptable_types = [
            ROOT.TH1C, ROOT.TH1S, ROOT.TH1I, ROOT.TH1F, ROOT.TH1D,
            ROOT.TH2C, ROOT.TH2S, ROOT.TH2I, ROOT.TH2F, ROOT.TH2D,
        ]

        if not isinstance(hist, tuple(acceptable_types)):
            raise TypeError(f'histogram must be one of {acceptable_types}')
        
        if isinstance(hist, ROOT.TH2):
            return self._hist2d_to_df(hist, *args, **kwargs)
        
        elif isinstance(hist, ROOT.TH1):
            return self._hist1d_to_df(hist, *args, **kwargs)
        
        else:
            raise TypeError(f'histogram must be one of {acceptable_types}')
        
    def _hist1d_to_df(self, hist, xname='x', yname='y'):
        df = dict()
        columns = [xname, yname, f'{yname}_err']
        getters = [
            hist.GetXaxis().GetBinCenter, 
            hist.GetBinContent,
            hist.GetBinError
        ]
        for col, get in zip(columns, getters):
            df[col] = [get(b) for b in range(1, hist.GetNbinsX() + 1)]

        df = pd.DataFrame(df)
        condition = (df[yname] != 0.0)
        df[f'{yname}_ferr'] = np.where(
            condition, 
            np.abs(df[f'{yname}_err']/df[yname]), 
            0.0
        )
        return df

    def _hist2d_to_df(self, histo, xname='x', yname='y', zname='z', keep_zeros=True):
        x = np.array([
            histo.GetXaxis().GetBinCenter(b) for b in range(1, histo.GetNbinsX() + 1)
        ])
        y = np.array([
            histo.GetYaxis().GetBinCenter(b) for b in range(1, histo.GetNbinsY() + 1)
        ])

        content = np.array(histo)
        content = content.reshape(len(x) + 2, len(y) + 2, order='F')
        content = content[1:-1, 1:-1]

        error = np.array([histo.GetBinError(b)
                         for b in range((len(x) + 2) * (len(y) + 2))])
        error = error.reshape(len(x) + 2, len(y) + 2, order='F')
        error = error[1:-1, 1:-1]

        xx, yy = np.meshgrid(x, y, indexing='ij')
        df = pd.DataFrame({
            xname: xx.flatten(),
            yname: yy.flatten(),
            zname: content.flatten(),
            f'{zname}_err': error.flatten(),
        })
        mask = (df[zname] != 0.0)
        df[f'{zname}_ferr'] = np.where(
            mask,
            np.abs(df[f'{zname}_err'] / df[zname]),
            0.0
        )
        return df if keep_zeros else df.query(f'{zname} != 0.0').reset_index(drop=True)

def histogram_to_df(hist, *args, **kwargs):
    return HistogramReader()(hist, *args, **kwargs)