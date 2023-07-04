import pathlib
import re
import pandas as pd
import numpy as np
from astropy import units
from collections import defaultdict
from typing import Literal

from pycpa.utilities import ame

class Particle:
    ALIAS = {
        'n': 'n1',
        'p': 'H1',
        'd': 'H2',
        't': 'H3',
        '3He': 'He3',
        '4He': 'He4',
        'neutron' : 'n1',
        'proton': 'H1',
        'deuteron': 'H2',
        'triton': 'H3',
        'Helium3': 'He3',
        'Helium4': 'He4',
        'alpha' : 'He4',
        'Alpha' : 'He4',
    }
    def __init__(self, name):
        if name in self.ALIAS:
            name = self.ALIAS[name]
        self.name = name
        self.mass = ame.get_mass(name)
        self.N, self.Z = ame.get_NZ(symbol=self.name)
    
class RunLog:
    COLUMNS = {
        'reaction' : str,
        'start-id' : int, 
        'end-id' : int,
        'badmap' : str,
        'shadowbar' : int,
        'trigger' : str,
    }
    def __init__(self, path : str = None):
        pth = pathlib.Path(__file__).parent.parent / 'database' / 'e15190' / 'RunInfo.dat' if path is None else path

        with open(str(pth), 'r') as f:
            # read non-empty lines that do not start with '#'
            content = [line.split() for line in f.readlines() if line.strip() and not line.startswith('#')]
            df = pd.DataFrame(content, columns=list(self.COLUMNS.keys())).astype(self.COLUMNS)
            df['nruns'] = df['end-id'] - df['start-id'] + 1
        self.df = df

    def get_number_of_runs(self, reaction : str) -> int:
        return self.df[self.df['reaction'] == reaction]['nruns'].sum()

class CollisionReaction:
    @staticmethod
    def dissemble_reaction(reaction, dtype:Literal['dict', 'list']='dict'):
        beamA, targetA, beam_energy = [
            int(m) for m in re.compile('[0-9]+').findall(reaction)]
        beam, target = [
            m for m in re.compile('[A-Za-z]{2}').findall(reaction)]
        
        if dtype == 'dict':
            return {
                'beam' : beam,
                'target' : target,
                'beamA' : beamA,
                'targetA' : targetA,
                'beam_energy' : beam_energy,
            } 
        elif dtype == 'list':
            return beam, target, beamA, targetA, beam_energy
        else:
            raise Exception('`dtype` must be `dict` or `list`.')

    @staticmethod
    def get_betacms(reaction):
        d = CollisionReaction.dissemble_reaction(reaction)
        beam_mass = ame.get_mass(d['beam'] + str(d['beamA'])) 
        target_mass = ame.get_mass(d['target'] + str(d['targetA'])) 
        beam_ke = d['beam_energy'] * d['beamA'] * units.MeV

        beam_energy_tot = beam_ke + beam_mass
        mom_beam = np.sqrt(beam_ke ** 2 + 2 * beam_ke * beam_mass)
        gamma = beam_energy_tot / beam_mass
        return (mom_beam / (gamma * beam_mass + target_mass)).value

    @staticmethod
    def get_rapidity_beam(reaction):
        d = CollisionReaction.dissemble_reaction(reaction)
        beam_mass = ame.get_mass(d['beam'] + str(d['beamA'])) 
        beam_ke = d['beam_energy'] * d['beamA'] * units.MeV

        mom_beam = np.sqrt(beam_ke ** 2 + 2 * beam_ke * beam_mass)
        return 0.5 * np.log((beam_ke + beam_mass + mom_beam) / (beam_ke + beam_mass - mom_beam)).value

    @staticmethod
    def get_rapidity_beam_and_target(reaction):
        d = CollisionReaction.dissemble_reaction(reaction)
        beam_mass = ame.get_mass(d['beam'] + str(d['beamA'])) 
        beam_ke = d['beam_energy'] * d['beamA'] * units.MeV

        mom_beam = np.sqrt(beam_ke ** 2 + 2 * beam_ke * beam_mass)

        target_mass = ame.get_mass(d['target'] + str(d['targetA']))
        
        return (0.5 * np.log((beam_ke + beam_mass + target_mass + mom_beam) / (beam_ke + beam_mass + target_mass - mom_beam))).value

        


        