"""

results.py

Author: Jordan Mirocha
Affiliation: McGill University
Created on: Tue 17 May 2022 11:13:01 EDT

Description:

"""

import pickle
import numpy as np

try:
    import h5py
except ImportError:
    pass

metadata = ['model_name', 'variant_ids', 'variants', 'chains']

class ResultsFinal(object):
    def __init__(self, data):
        self.raw_data = data

class ResultsMCMC(object):
    def __init__(self, data):
        self.raw_data = data

    @property
    def is_chain(self):
        if not hasattr(self, '_is_chain'):
            self._is_chain = self.data['chains']
        return self._is_chain

    @property
    def data(self):
        """
        Re-format chains: convert to dictionary if they aren't already.
        """

        if not hasattr(self, '_data'):
            if type(self.raw_data) is dict:
                self._data = self.raw_data
            else:
                raise NotImplemented('Only know how to handle dictionaries so far.')

        return self._data

    @property
    def variants(self):
        """
        Names of variants of main model result.
        """
        if not hasattr(self, '_variants'):
            if 'variants' not in self.raw_data.keys():
                self._variants = [None]
            elif self.raw_data['variants'] is None:
                self._variants = [None]
            else:
                self._variants = self.raw_data['variants']

        return self._variants

    @property
    def params(self):
        if not hasattr(self, '_params'):
            params = []
            for redshift in self.redshifts:
                for key in self.data[redshift].keys():
                    if key in params:
                        continue

                    params.append(key)

            self._params = params

        return self._params

    @property
    def redshifts(self):
        if not hasattr(self, '_redshifts'):
            all_keys = list(self.data.keys())
            _redshifts = []
            for key in all_keys:
                if key in metadata:
                    continue
                _redshifts.append(key)

            self._redshifts = np.sort(_redshifts)

        return self._redshifts

    def get_samples(self, param, redshift):
        """
        Return 1-D array of samples for input `param`.
        """

        if param not in self.params:
            raise KeyError("No param {} in this dataset! Options={}".format(
                param, self.params
            ))

        rdata = self.data[redshift]

        if self.is_chain and type(rdata) is str:
            with open(rdata, 'rb') as f:
                data = pickle.load(f)

            raise NotImplemented('help')
            pars, z = data['pinfo']
            #k =
            return data['chain'][k]
        else:
            return rdata[param]

    def get_limits(self, param, redshift, two_tailed=True):
        """
        Return list of limits.
        """

        if self.data[redshift] is None:
            return None, None

        if self.is_chain:
            d = self.get_samples(param, redshift)

            conf = [2.5, 16] if two_tailed else [5, 32]
            lims = [np.percentile(d, _c_) for _c_ in conf]
        else:
            _lims = self.data[redshift][param]

            if len(_lims) == 4:
                lims = [_lims[0], _lims[2]]
            elif len(_lims) == 2:
                lims = _lims
            else:
                raise ValueError("Do not understand length of data ({})".format(
                    len(_lims)
                ))

        return lims
