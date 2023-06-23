import pandas as pd
import numpy as np

from SEAT.Serpent2UncertaintyPropagation.MomentumPropagation.conversions import *

__all__ = [
    "Covariance"
    ]

def get_cov(df):  ## sets numerical types to the column and index levels
    df_ = df.copy()
    MAT = df_.columns.get_level_values(0).astype(int)
    MT = df_.columns.get_level_values(1).astype(int)
    ENE = np.array([np.array(i.split(", ")).astype(float) for
                    i in df_.columns.get_level_values(2).str.replace('(', "").str.replace(']', "")])[:, 1]
    idx = pd.MultiIndex.from_arrays([MAT, MT, ENE]).set_names(["MAT", "MT", "E [eV]"]) ## E as float (right of interval)
    df_.columns = idx
    df_.index = idx
    return df_

def reindex_cov(df):
    df_ = df.copy()
    idx = pd.MultiIndex.from_arrays([df_.index.get_level_values(0).map(MATs),
                                     df_.index.get_level_values(1).map(MTs),
                                     df_.index.get_level_values(2).map(ECCO33)]).set_names(["N", "MT", "E [MeV]"])
    df_.columns = idx
    df_.index = idx
    return df_

class Covariance(pd.DataFrame):
    @classmethod
    def from_csv(cls, paths: list, *args, **kwargs):
        data = pd.concat([get_cov(pd.read_csv(path, *args, **kwargs)) for path in paths])
        return cls(reindex_cov(data).fillna(0))

    @property
    def idx(self):
        return self.index
