import numpy as np
import pandas as pd


class histogram_handler:
    @staticmethod
    def projection(df, cut=None, range=None, bins=None, axis='x', normalize=True):
        """ Project a 2D histogram (in terms of pd.DataFrame) along an axis
        Parameters
        ----------
        df : pd.DataFrame
            Dataframe containing the histogram to be projected. 
        cut : tuple, optional
            Range of the axis to be cut. If None, no cut is applied.
        range : tuple, optional
            Range of the axis to be projected. If None, the range of the dataframe is used.
        bins : int, optional
            Number of bins to be used. If None, the number of bins of the dataframe is used.
        axis : str, optional
            Axis along which the projection is performed. Must be 'x' or 'y'.
        normalize : bool, optional
            If True, the projection is normalized to the bin width.
        Returns
        -------
        pd.DataFrame
            Dataframe containing the projected histogram.
        """

        df = df.copy()
        if df.shape[1] < 4:
            raise ValueError('Dataframe must have at least 3 columns')
        elif df.shape[1] == 4:
            df.columns = ['x', 'y', 'z', 'z_err']
        elif df.shape[1] == 5:
            df.columns = ['x', 'y', 'z', 'z_err', 'z_ferr']
        else:
            raise ValueError('Dataframe must have at most 5 columns')

        cut_axis = 'y' if axis == 'x' else 'x'
        if cut is not None:
            df = histogram_handler.query(df, cut, cut_axis)

        range = (df[axis].min(), df[axis].max()) if range is None else range
        bins = df.shape[0] if bins is None else bins

        hist, edges = np.histogram(
            df[axis],
            weights=df['z'],
            range=range,
            bins=bins
        )

        histerr, _ = np.histogram(
            df[axis],
            weights=df['z_err']**2,
            range=range,
            bins=bins
        )
        histerr = np.sqrt(histerr)
        norm = (np.diff(range) / bins)[0] if normalize else 1.

        return pd.DataFrame({
            'x': (edges[1:] + edges[:-1])/2,
            'y': hist / norm,
            'y_err': histerr / norm,
            'y_ferr': np.divide(histerr, hist, out=np.zeros_like(hist), where=hist != 0)
        })

    @staticmethod
    def findbin(df, x, axis='x'):
        """ Find the bin for a given array in the dataframe-histogram in the same way as FindBin in ROOT.
        Parameters
        ----------
        df : pd.DataFrame
            Dataframe containing the histogram, with the bin centers in the specified axis.
        x : array
            Array of values to be binned.
        axis : str, optional    
            Axis along which the projection is performed. Must be 'x' or 'y'.
        Returns
        -------
        array
            Array of bin indices.
        """
        bincenters = np.unique(df[axis])
        binwidth = bincenters[1] - bincenters[0]
        bins = np.array([bincenters[0] - binwidth / 2] +
                        list(bincenters + binwidth / 2))
        return np.digitize(x=x, bins=bins, right=False)

    @staticmethod
    def query(df, range, axis='x'):
        df = df.copy()
        # find the bin (start from 1), subtract 1 to get the index of array
        bin1, bin2 = histogram_handler.findbin(df, range) - 1
        x = np.unique(df[axis])
        range = (x[bin1], x[bin2])
        return df.query(f'{axis} >= {range[0]} & {axis} <= {range[1]}')

    @staticmethod
    def rebin1D(df, range, bins=None, normalize=True):
        df = df.copy()
        if df.shape[1] > 4:
            raise ValueError('Dataframe must have at most 4 columns')
        elif df.shape[1] == 4:
            df.columns = ['x', 'y', 'y_err', 'y_ferr']
        elif df.shape[1] == 3:
            df.columns = ['x', 'y', 'y_err']
        elif df.shape[1] == 2:
            df.columns = ['x', 'y']
            df['y_err'] = 0
        else:
            raise ValueError('Dataframe must have at least 2 columns')

        bins = df.shape[0] if bins is None else bins
        hist, edges = np.histogram(
            df['x'],
            weights=df['y'],
            range=range,
            bins=bins
        )

        histerr, _ = np.histogram(
            df['x'],
            weights=df['y_err']**2,
            range=range,
            bins=bins
        )
        histerr = np.sqrt(histerr)
        norm = (np.diff(range) / bins)[0] if normalize else 1.

        return pd.DataFrame({
            'x': (edges[1:] + edges[:-1])/2,
            'y': hist / norm,
            'y_err': histerr / norm,
            'y_ferr': np.divide(histerr, hist, out=np.zeros_like(hist), where=hist != 0)
        })

    @staticmethod
    def mul(a, b, axis='y'):
        """ Multiply two 1D histograms (in terms of pd.DataFrame). Currently only implemented for 1D histograms as we rarely multiply multi-dimensional histograms.
        Parameters
        ----------
        a : pd.DataFrame
            Dataframe containing the first histogram.
        b : pd.DataFrame
            Dataframe containing the second histogram.
        axis : str, optional    
            Axis along which the multiplication is performed, depending on the name of the columns. Default is 'y'.
        Returns
        -------
        pd.DataFrame
            Dataframe containing the product of the two histograms.
        """
        if not set(['x', 'y', 'y_err']).issubset(set(a.columns)) and set(['x', 'y', 'y_err']).issubset(set(b.columns)):
            raise ValueError('Dataframes must have at least 3 columns')

        return pd.DataFrame({
            'x': a['x'],
            axis: a[axis] * b[axis],
            f'{axis}_err': np.sqrt((a[f'{axis}_err'] * b[axis])**2 + (a[axis] * b[f'{axis}_err'])**2),
            f'{axis}_ferr': np.sqrt((a[f'{axis}_ferr'] * b[axis])**2 + (a[axis] * b[f'{axis}_ferr'])**2)
        })

    @staticmethod
    def div(a, b, axis='y'):
        """ Divide two 1D histograms (in terms of pd.DataFrame). Currently only implemented for 1D histograms as we rarely divide multi-dimensional histograms.
        Parameters
        ----------
        a : pd.DataFrame
            Dataframe containing the first histogram.
        b : pd.DataFrame
            Dataframe containing the second histogram.
        axis : str, optional
            Axis along which the division is performed, depending on the name of the columns. Default is 'y'.
        Returns
        -------
        pd.DataFrame
            Dataframe containing the quotient of the two histograms.
        """
        if not set(['x', 'y', 'y_err']).issubset(set(a.columns)) and set(['x', 'y', 'y_err']).issubset(set(b.columns)):
            raise ValueError('Dataframes must have at least 3 columns')

        if not 'y_ferr' in a.columns:
            a['y_ferr'] = np.divide(
                a['y_err'], a['y'], out=np.zeros_like(a['y']), where=a['y'] != 0)
        if not 'y_ferr' in b.columns:
            b['y_ferr'] = np.divide(
                b['y_err'], b['y'], out=np.zeros_like(b['y']), where=b['y'] != 0)

        y = np.divide(a[axis], b[axis], out=np.zeros_like(
            a[axis]), where=b[axis] != 0)
        y_err = y * np.sqrt(a[f'{axis}_ferr'] ** 2. + b[f'{axis}_ferr'] ** 2.)
        return pd.DataFrame({
            'x': a['x'],
            axis: y,
            f'{axis}_err': y_err,
            f'{axis}_ferr': np.divide(y_err, y, out=np.zeros_like(y), where=y != 0)
        })

    @staticmethod
    def sum(a, b, axis='y', sign=1.):
        """ Get the sum of two 1D histograms (in terms of pd.DataFrame). Currently only implemented for 1D histograms as we rarely add multi-dimensional histograms.
        Parameters
        ----------
        a : pd.DataFrame
            Dataframe containing the first histogram.
        b : pd.DataFrame
            Dataframe containing the second histogram.
        axis : str, optional
            Axis along which the addition is performed, depending on the name of the columns. Default is 'y'.
        sign : float, optional
            Sign of the second histogram. Default is 1. 
        Returns
        -------
        pd.DataFrame
            Dataframe containing the sum of the two histograms.
        """

        if not set(['x', 'y', 'y_err']).issubset(set(a.columns)) and set(['x', 'y', 'y_err']).issubset(set(b.columns)):
            raise ValueError('Dataframes must have at least 3 columns')

        y = a[axis] + b[axis] * sign
        y_err = np.sqrt(a[f'{axis}_err']**2 + b[f'{axis}_err']**2)
        return pd.DataFrame({
            'x': a['x'],
            axis: y,
            f'{axis}_err': y_err,
            f'{axis}_ferr': np.divide(y_err, y, out=np.zeros_like(y), where=y != 0)
        })

    @staticmethod
    def add(a, b, axis='y'):
        return histogram_handler.sum(a, b, axis, sign=1.)

    @staticmethod
    def sub(a, b, axis='y'):
        return histogram_handler.sum(a, b, axis, sign=-1.)
