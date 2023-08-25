import pandas as pd
import numpy as np
import serpentTools as sts
import SEAT.nuclides

import uncertainties
from uncertainties import unumpy

from SEAT.Serpent2UncertaintyPropagation.MomentumPropagation.conversions import *

__all__ = [
    "Sensitivity"
    ] 

def get_zais(reader) -> list:
    """
    Parameters
    ----------
    reader: sts.SensitivityReader
        the reader where the sensitivity results are stored.

    Returns
    -------
    list of the perturbed zais as nuclide symbols

    """
    zais = [SEAT.nuclides.zam2nuclide(i) for i in reader.zais]
    return zais

def get_perts_mt(reader) -> list:
    """
    Parameters
    ----------
    reader: sts.SensitivityReader
        the reader where the sensitivity results are stored.

    Returns
    -------
    list of the perturbed reactions as reaction names

    """
    perts = [MT2NAME[MTLIST2MT[i]] for i in reader.perts]
    return perts

def get_perts_rea(reader):
    """
    Parameters
    ----------
    reader: sts.SensitivityReader
        the reader where the sensitivity results are stored.

    Returns
    -------
    list of the perturbed reactions as reaction names

    """
    perts = [MT2NAME[REALIST2MT[i]] for i in reader.perts]
    return perts

def get_ene(reader):
    """
    Parameters
    ----------
    reader: sts.SensitivityReader

    Returns
    -------
    list of the sensitivity energy groups in [MeV]

    """
    ene = reader.energies[1:]
    return ene

def sts2df(reader, drop_total: bool=True, observe=None) -> pd.DataFrame:
    """
    Transforms the serpentTools sensitivity reader to a
    pandas.DataFrame object.

    Parameters
    ----------
    reader: sts.SensitivityReader
        the reader where the sensitivity results are stored.
    derop_total: bool, optional
        removes total sensitivities from the data.
    observe: str or list, optional
        the observable of the output sensitivity.

    Returns
    -------
    pd.DataFrame of the organized sensitivity vector

    """
    zai = get_zais(reader)
    try:
        pert = get_perts_mt(reader)
    except KeyError:
        pert = get_perts_rea(reader)        
    ene = get_ene(reader)

    out = []
    for obs, sens in reader.sensitivities.items():
        if obs in observe:
            out.append(pd.DataFrame(sens.reshape(len(zai) * len(pert) * len(ene), 2),
                       index=pd.MultiIndex.from_product([zai, pert, ene], names=["N", "MT", "E [MeV]"]),
                       columns=["S", "stat. err."]).reset_index().assign(Observable=obs))
    # k, b = None, None
    # if 'keff' in reader.sensitivities.keys() and 'keff' == observable:
    #     k = pd.DataFrame(reader.sensitivities["keff"].reshape(len(zai) * len(pert) * len(ene), 2),
    #                      index=pd.MultiIndex.from_product([zai, pert, ene], names=["N", "MT", "E [MeV]"]),
    #                      columns=["S", "stat. err."]).reset_index().assign(Observable="keff")
    # if 'beff' in reader.sensitivities.keys() and 'keff' == observable:
    #     b = pd.DataFrame(reader.sensitivities["beff"].reshape(len(zai) * len(pert) * len(ene), 2),
    #                      index=pd.MultiIndex.from_product([zai, pert, ene], names=["N", "MT", "E [MeV]"]),
    #                      columns=["S", "stat. err."]).reset_index().assign(Observable="beff")
    
    out = pd.concat(out, ignore_index=True)
    if drop_total:
        out = out.query('~MT.str.contains("total").values')
    return out.set_index(["Observable", "N", "MT", "E [MeV]"])


class Sensitivity:
    """
    Container of sensitivity information.

    Names will be:
        Index: 'Observable', 'N', 'MT', 'E [MeV]'
        Columns: 'S', 'stat. err.' (in relative terms)
        Columns: 'nom_val' and 'abs_sdev' will be created when accessing 'data'
        Other columns can be added without posing problems
    """
    _csv = False
    _csv_index = None
    _df = None
    _idx_names = ['Observable', 'N', 'MT', 'E [MeV]']
    _col_names = ['S', 'stat. err.', 'nom_val', 'abs_sdev']  # colum names not necessarily limited to this

    def __init__(self, file: str, drop_total: bool=True, observable: str=None, concat=None):
        self.file = file
        self.drop_total = drop_total
        self.observable = observable
        self.concat = concat  # things to add to the sensitivity vector in S.data

    @classmethod
    def from_csv(cls, *args, index: list, **kwargs):
        instance = cls(*args, **kwargs)
        instance._csv = True
        instance._csv_index = index
        return instance

    @classmethod
    def from_df(cls, df: pd.DataFrame):
        instance = cls(file=None)
        instance._df = df.copy()
        return instance

    @property
    def idx(self):
        if self._df is not None:
            idx = self._df.index
        elif self._csv:
            idx = self.data.index
        else:
            obs = [self.observable]
            zai = get_zais(self.reader)
            try:
                pert = get_perts_mt(self.reader)
            except KeyError:
                pert = get_perts_rea(self.reader)        
            ene = get_ene(self.reader)
            idx = pd.MultiIndex.from_product([obs, zai, pert, ene],
                                          names=self._idx_names)
        if self.concat is not None:
            idx = idx.append(self.concat.index)
        if self.drop_total:
            idx = pd.MultiIndex.from_frame(
                idx.to_frame(index=False).query("MT != 'total'"))
        return idx

    @property
    def data(self):
        if self._df is not None:
            df = self._df.copy()
        elif self._csv:
            df = pd.read_csv(self.file).drop("Unnamed: 0", axis=1)
            if self.drop_total:
                df = df.query('~MT.str.contains("total").values')
            if self._csv_index is not None:
                df = df.set_index(self._csv_index)
        else:
            df = sts2df(sts.read(self.file), self.drop_total, [self.observable])
        if self.concat is not None:
            df = pd.concat([df, self.concat])
        # uncertain sensitivities with uncertainties.py if if not already
        if not isinstance(df[self._col_names[0]].iloc[0],
                          uncertainties.core.Variable):
            df[self._col_names[0]] = unumpy.uarray(df[self._col_names[0]],
                                                   df[self._col_names[1]] \
                                                       * df[self._col_names[0]])
        # create nom_val e abs_sdev columns if not already there
        if self._col_names[2] not in df.columns:
            df[self._col_names[2]] = df[self._col_names[0]].apply(lambda x: x.n)
        if self._col_names[3] not in df.columns:
            df[self._col_names[3]] = df[self._col_names[2]] * df[self._col_names[1]]
        return df

    @property
    def reader(self):
        return sts.read(self.file)

    @property
    def observables(self):
        return sts.read(self.file).sensitivities.keys()

    @property
    def upper(self):
        up = self.data.copy()
        up[self._col_names[0]] = (up[self._col_names[0]] +
                                  up[self._col_names[3]]).apply(
                                      AffineScalarFun_2_Variable)
        return self.from_df(up)

    @property
    def lower(self):
        low = self.data.copy()
        low[self._col_names[0]] = (low[self._col_names[0]] -
                                   low[self._col_names[3]]).apply(
                                       AffineScalarFun_2_Variable)
        return self.from_df(low)

    def observe(self, observable: str):
        df = self.data.reset_index().query(
                self._idx_names[0] + " == @observable").set_index(self._idx_names)
        return self.__class__.from_df(df)

    def query(self, query):
        return self.__class__.from_df(self.data.query(query))
