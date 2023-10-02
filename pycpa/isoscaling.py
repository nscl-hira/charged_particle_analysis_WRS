import sys
import numpy as np
import warnings
from iminuit import Minuit
from iminuit.cost import LeastSquares
from typing import Literal

class isoscaling:
    def __init__(self, method:Literal['iterate','free']='iterate', **kwargs):
        self.method = method
        kw = {
            # minuit fitter settings
            'minuit_ncall' : None, 
            'minuit_iterate' : 5, 
            'minuit_verbose' : 1,
            'minuit_throw_nan' : True,

            # works only for method='iterate'
            'norm_step_size' : 0.01,
            'norm_range' : (0.70, 1.2),

            # for debugging purpose
            'catch_runtime_warning' : False,
        }
        kw.update(kwargs)
        for key, value in kw.items():
            setattr(self, key, value)

    @staticmethod
    def model(X:np.ndarray, norm:float, alpha:float, beta:float):
        """
        Parameters
        -----------
        X : np.ndarray of shape (2, nparticles)
            [n, z] of each particle used in the fit
        """
        return norm * np.exp(X[0] * alpha + X[1] * beta)
    
    @staticmethod
    def model_error(X, norm, alpha, beta, norm_err, alpha_err, beta_err):
        return isoscaling.model(X, norm, alpha, beta) * np.sqrt(X[0] ** 2 * alpha_err ** 2 + X[1] ** 2 * beta_err**2 + norm_err ** 2 / norm ** 2)

    def _fit_single(
        self, 
        X:np.ndarray, 
        y:np.ndarray, 
        yerr:np.ndarray, 
        initvals:list[float]=[1.0,0.5,-0.5], 
        fixed:list[bool]=[True,False,False],
        limits:list[list[float]]=[[0.,None],[None,None],[None,None]],
    ) -> Minuit:
        
        y, yerr = y.flatten(), yerr.flatten()
        assert y.shape == yerr.shape == (X.shape[1],), 'y and yerr must have the same shape'

        if self.catch_runtime_warning:
            # custom cost function to handle RuntimeWarning which usually occurs 
            # when the minuit.limits is too borad
            def cost_fcn(norm, alpha, beta):
                try : 
                    warnings.filterwarnings("error")
                    d = y - isoscaling.model(X, norm, alpha, beta)
                    e = yerr
                    return np.sum(np.divide(d, e, where=(e != 0.), out=np.zeros_like(y)) ** 2)
                except RuntimeWarning:
                    raise RuntimeError(f'RuntimeWarning, {alpha}, {beta}, {norm}')
                finally:
                    warnings.filterwarnings("default")
        else:
            cost_fcn = LeastSquares(X, y, yerr, isoscaling.model)

        minuit = Minuit(cost_fcn, norm=initvals[0], alpha=initvals[1], beta=initvals[2])
        minuit.limits = limits
        for idx, para in enumerate(minuit.params):
            minuit.fixed[para.name] = fixed[idx]

        minuit.throw_nan = self.minuit_throw_nan
        minuit.print_level = self.minuit_verbose

        minuit.migrad(self.minuit_ncall, self.minuit_iterate)
        minuit.hesse()
        return minuit

    def _find_optimum_normalization(
        self, 
        X:np.ndarray, 
        y:np.ndarray, 
        yerr:np.ndarray,
        **kwargs
    ) -> float:
        optimum_chi2 = sys.float_info.max
        optimum_norm = self.norm_range[0]

        kw = {
            'initvals' : [1.0, 0.5, -0.5],
            'fixed' : [True,False,False],
            'limits' : [[None,None],[None,None],[None,None]],
        }
        kw.update(kwargs)

        for n in np.arange(*self.norm_range, self.norm_step_size):
            kw['initvals'][0] = n
            kw['fixed'][0] = True
            chi2 = np.sum([
                self._fit_single(
                    X, 
                    y[idx],
                    yerr[idx],
                    **kw
                ).fval for idx in range(y.shape[0])
            ])
            if chi2 < optimum_chi2:
                optimum_chi2 = chi2
                optimum_norm = n

        return optimum_norm

    def _fit_iterate(self, X, y, yerr, **kwargs):
        norm = self._find_optimum_normalization(X, y, yerr, **kwargs)
        kwargs['initvals'][0] = norm
        kwargs['fixed'][0] = True
        kwargs['limits'][0] = [norm, norm]
        return self._fit(X, y, yerr, **kwargs)

    def _fit(self, X:np.ndarray, y:np.ndarray, yerr:np.ndarray, **kwargs):
        minuits = [self._fit_single(X, y[n], yerr[n], **kwargs) for n in range(self.npoints)]
        self.chi2 = [m.fval for m in minuits]
        self.param_names = minuits[0].parameters
        self.param_values = {
            name : [m.values[self.param_names.index(name)] for m in minuits] for name in self.param_names
        }
        self.param_errors = {
            name : [m.errors[self.param_names.index(name)] for m in minuits] for name in self.param_names
        }
        self.minuits = minuits
        self.param_covariance = [np.array(m.covariance) for m in minuits]
        self.valid = [m.valid for m in minuits]
        self.accurate = [m.accurate for m in minuits]
        for n in range(self.npoints):
            if not self.valid[n]:
                warnings.warn(f'fit for point {n} is not valid')
            if not self.accurate[n]:
                warnings.warn(f'fit for point {n} is not accurate')
        return self

    def fit(
        self, 
        X:np.ndarray, 
        y:np.ndarray, 
        yerr:np.ndarray, 
        initvals:list[float]=[1.0,0.5,-0.5], 
        fixed:list[bool]=[False,False,False],
        limits:list[list[float]]=[[None,None],[None,None],[None,None]],
    ):
        """ perform a global isoscaling fit for all data points. This function calls `self._fit` and will create the following attributes:
        - `self.chi2` : list of chi2 values
        - `self.param_names` : names of the parameters
        - `self.param_values` : dict[str, list], values of the parameters
        - `self.param_errors` : dict[str, list], errors of the parameters
        - `self.minuits` : list of Minuit objects
        - `self.param_covariance` : list of covariance matrices
        - `self.valid` : list of bool, whether the fit is valid
        - `self.accurate` : list of bool, whether the fit is accurate

        Parameters
        -----------
        X : np.ndarray of shape (2, nparticles)
            [n, z] of each particle used in the fit
        y : np.ndarray of shape (npoints, nparticles)
            y values of each particle used in the fit
        yerr : np.ndarray of shape (npoints, nparticles)
            yerr values of each particle used in the fit
        initvals : list of float
            initial values of [norm, alpha, beta]
        fixed : list of bool
            whether to fix [norm, alpha, beta]. If `self.method == 'iterate'`, `fixed[0]` will be forced to be `True`.
        limits : list of list of float
            limits of [norm, alpha, beta], default to be [None, None, None]. If `self.method == 'iterate'`, `limits[0]` will be forced to be [norm, norm]. where `norm` is the optimum normalization found by `self._find_optimum_normalization`.
        """
        assert y.shape == yerr.shape and y.shape[1] == X.shape[1], 'y, yerr and X must have the same shape'
        self.nparticles = X.shape[1]
        self.npoints = y.shape[0]
        method = self._fit_iterate if self.method == 'iterate' and self.npoints > 1 else self._fit
        return method(X, y, yerr, initvals=initvals, fixed=fixed, limits=limits)
    
    def _sample_error(self, X:np.ndarray, nsamples:int=1000):
        """ sample the error of the model parameters using the covariance matrix
        Parameters
        ----------
        X : np.ndarray of shape (2, nparticles)
            [n, z] of each particle used in the fit
        """
        samples = []
        for n in range(self.npoints):
            samples_ = np.random.multivariate_normal(
                mean=[self.param_values[name][n] for name in self.param_names], 
                cov=self.param_covariance[n],
                size=nsamples
            )
            # get_range = lambda x : (np.min(x, axis=0), np.max(x, axis=0))
            get_range = lambda x : 0.5 * (np.max(x, axis=0) - np.min(x, axis=0))
            samples.append(get_range([isoscaling.model(X, *sample) for sample in samples_]))
        return np.array(samples)

    def predict(self, X:np.ndarray, return_std=False, method:Literal['bootstrap','propagate']='bootstrap', **kwargs) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
        
        vals = np.transpose([vals for vals in self.param_values.values()])
        y = np.array([isoscaling.model(X, *val) for val in vals])
        
        if return_std:
            assert method in ['bootstrap', 'propagate'], 'method must be either `bootstrap` or `propagate`'
            if method == 'bootstrap':
                std = self._sample_error(X, **kwargs)
            else:
                errs = np.transpose([errs for errs in self.param_errors.values()])
                std = np.array([isoscaling.model_error(X, *val, *err) for val,err in zip(vals, errs)])
            return y, std
        return y
