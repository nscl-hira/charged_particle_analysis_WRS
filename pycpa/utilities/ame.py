import re
import pathlib
import requests
import pandas as pd
from astropy import units, constants
from pycpa import PROJECT_DIR

class AME:
    DEFAULT_PATH = PROJECT_DIR / 'database' / 'ame' / 'mass20.txt'
    def __init__(self, url=None, dl_path=None):
        if url is None:
            url = 'https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt'

        r = requests.get(url)
        dl_path = self.DEFAULT_PATH if dl_path is None else pathlib.Path(dl_path)
        dl_path.parent.mkdir(exist_ok=True, parents=True)

        with open(str(dl_path), 'w') as file:
            file.write(r.text)
        
        self.dl_path = dl_path
        self._parse_mass_table()

    def _parse_mass_table(self):
        with open(str(self.dl_path), 'r') as file:
            content = file.readlines()

        # record fortran format
        for line in content:
            if 'format' in line:
                fortran_format = line.split(':')[1].strip().split(',')
                break

        # remove headers
        for i, line in enumerate(content):
            if line.count('keV') == 3 and line.count('micro-u') == 1:
                break
        
        content = content[i+1:]

        # '#' means estimated value (non-experimental)
        # '*' means unable to calculate
        content = [re.sub('[*#]', ' ', line) for line in content]

        columns_names = {
            '1' : str,
            'N-Z' : int,
            'N' : int,
            'Z' : int,
            'A' : int,
            'EL' : str, 
            'O' : str,
            'mass_excess' : float,
            'mass_excess_err' : float,
            'binding_energy_per_nucleon' : float,
            'binding_energy_per_nucleon_err' : float,
            'B-' : str,
            'beta_decay_energy' : float,
            'beta_decay_energy_err' : float,
            'au' : int,
            'atomic_mass_au' : float,
            'atomic_mass_err' : float,
        }

        spacing = [int(re.findall(r'\d+', format_string)[0]) for format_string in fortran_format]
        
        parsed_content = []
        for line in content:
            start = 0
            parsed_line = []
            for space in spacing:
                parsed_line.append(line[start : start + space])
                start += space
            parsed_content.append(parsed_line)

        df = pd.DataFrame(parsed_content)
        df = df.loc[:, ~df.apply(lambda x: x.str.isspace()).all()]
        df.columns = list(columns_names.keys())

        for column, dtype in columns_names.items():
            df[column] = df[column].str.strip()
            try: 
                df[column] = df[column].astype(dtype)
            except ValueError:
                df[column] = pd.to_numeric(df[column], errors='coerce')

        # parse by physics
        df['symbol'] = list(map(lambda s: s.lower().strip(), df['EL']))
        df['symbol'] = df['symbol'] + df['A'].astype(str)
        
        df['mass_excess'] = df['mass_excess'] * units.keV.to('MeV', 1.)
        df['mass_excess_err'] = df['mass_excess_err'] * units.keV.to('MeV', 1.)

        df['au'] = df['A'] + df['atomic_mass_au'] * units.micron.to('m', 1.)
        df['au_err'] = df['atomic_mass_err'] * units.micron.to('m', 1.)

        df['mass'] = df['au'] * (units.u * constants.c ** 2).to('MeV')
        df['mass_err'] = df['au_err'] * (units.u * constants.c ** 2).to('MeV')

        df['binding_energy_per_nucleon'] = df['binding_energy_per_nucleon'] * \
            units.keV.to('MeV', 1.)
        df['binding_energy_per_nucleon_err'] = df['binding_energy_per_nucleon_err'] * \
            units.keV.to('MeV', 1.)

        df.drop(columns=['1', 'N-Z', 'O', 'B-', 'EL', 'atomic_mass_au', 'atomic_mass_err'], inplace=True)
        self.df = df

    def _get_row(self, symbol=None, N=None, Z=None):
        if N is not None and Z is not None:
            return self._get_row_from_NZ(N, Z)
        elif symbol is not None:
            return self._get_row_from_symbol(symbol)
        else:
            raise ValueError('Either symbol or N and Z must be specified.')

    def _get_row_from_symbol(self, symbol):
        return self.df.loc[self.df['symbol'] == symbol.lower()]

    def _get_row_from_NZ(self, N, Z):
        return self.df.loc[((self.df['N'] == N) & (self.df['Z'] == Z))]

    def get_symbol(self, N, Z):
        return self._get_row(N=N, Z=Z)['symbol'].values[0]

    def get_NZ(self, symbol):
        row = self._get_row(symbol)
        return (row['N'].values[0], row['Z'].values[0])

    def get_mass(self, symbol='', N=None, Z=None, unit='MeV'):
        row = self._get_row(symbol, N, Z)
        return row['mass'].values[0] * units.MeV.to(unit)

    def get_binding_energy(self, symbol=None, N=None, Z=None, unit='MeV'):
        return self._get_row(symbol, N, Z)['binding_energy_per_nucleon'].values[0] * units.MeV.to(unit)

    def get_reduced_dataframe(
            self, 
            fname='ame_mass.txt', 
            usecols=['symbol', 'Z', 'A', 'mass', 'binding_energy_per_nucleon']
        ) -> pd.DataFrame:
        """ Output the mass table to a format convenient for reading in C++.
        Parameters
        ----------
        fname : str
            Output file name.
        usecols : list
            Columns to be output.
        """
        if not set(usecols).issubset(self.df.columns):
            raise ValueError(
                'usecols must be a subset of the columns of the AME table.')

        df = self.df[usecols].copy()
        df.to_csv(fname, sep=' ', columns=usecols, index=False)

        return df