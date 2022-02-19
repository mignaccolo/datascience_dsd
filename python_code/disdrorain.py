# This program contains the implementation of the disdrorain and disdrorain_2dvd classes
# to process drop count disdrometer data and 2dvd drop by drop disdrometer data
#
# Author: Massimiliano Ignaccolo
#         massimiliano.ignaccolo [at] tutanota.com. massimiliano.ignaccolo [at] sas.com
# Last modified: Feb 01 2022


# ------------- Necessary Python packages -START
import pandas as pd
import copy as cp
import numpy as np
# ------------- Necessary Python packages -END

# class used to allow dynamic initialization of attributes
class LazyProperty(object):
    def __init__(self, func):
        self._func = func
        self.__name__ = func.__name__
        self.__doc__ = func.__doc__

    def __get__(self, obj, klass=None):
        if obj is None:
            return None
        result = obj.__dict__[self.__name__] = self._func(obj)
        return result


# class to process disdrometer data that are divided in classes (e.g. non 2DVD data)
class disdrorain(object):
    def __init__(self, classpath=None, datapath=None, dataframe=None, fieldsep=' ',
                 instrument_area=5000, time_interval=60,
                 Aplawspeed=3.776, Bplawspeed=0.67,
                 Aexpospeed=9.65, Bexpospeed=10.3, Cexpospeed=0.6,
                 ncells=25):
        """
        didrorain class
        # classpath = path to disdrometer class limits file
            if None by default the RD80 Valdvogel class limits are used
        # datapath = path to file
        # dataframe = name of data frame
        # fieldsep = string separating field data
        # instrument_area = catchement area of disrometer in mm^2.
            Default value is 5000 = 50cm^2
        # time_interval = time interval of disdrometer recording in seconds.
            Default value is 60 = 1 minute
        # v(D) = Aplawspeed*D^Bplawspeed -- terminal velocity of drop of diameter D (power law velocity)
        # v(D) = Aexpospeed-Bexpospeed*exp(-Cexpospeed*D) -- terminal velocity of drop of diameter D ("exponenetial" law velocity)
        # ncells = number of cells in which to divide a disdrometer class when calculating moments for a "cloud" distribution
            under the "exponential" law velocity. This division is necessary to ensure that inside each cell the speed can
            be considered  constant
        """

        # default diameter classes: RD80 Valdvogel
        default_d = {'C1': [0.313, 0.405], 'C2': [0.405, 0.505], 'C3': [0.505, 0.596], 'C4': [0.596, 0.715],
                     'C5': [0.715, 0.827], 'C6': [0.827, 0.999], 'C7': [0.999, 1.232], 'C8': [1.232, 1.429],
                     'C9': [1.429, 1.582], 'C10': [1.582, 1.748], 'C11': [1.748, 2.077], 'C12': [2.077, 2.441],
                     'C13': [2.441, 2.727], 'C14': [2.727, 3.011], 'C15': [3.011, 3.385], 'C16': [3.385, 3.704],
                     'C17': [3.704, 4.127], 'C18': [4.127, 4.573], 'C19': [4.573, 5.145], 'C20': [5.145, 5.601]}
        # if we are reading disdrometer data from csv file
        if datapath is not None:
            self.data = pd.read_csv(datapath, sep=fieldsep, header=None)
        # if disdrometer data are already in a data frame
        if dataframe is not None:
            self.data = dataframe.copy()
        self.data.rename(columns=lambda x: 'C' + str(x + 1), inplace=True)
        self.data.index.name = 'record number'
        # if path to a csv file containing the limits of each class is not specified then use
        # default (RD80) class limits
        if classpath is None:  # use default diameter classes
            self.classlimits = pd.DataFrame(data=default_d)
            self.classlimits.rename(index={0: 'left', 1: 'right'}, inplace=True)
        else:  # load class limits from file
            self.classlimits = pd.read_csv(classpath, sep=fieldsep, header=None)
            self.classlimits.rename(columns=lambda x: 'C' + str(x + 1), index={0: 'left', 1: 'right'}, inplace=True)
        self.classlimits.index.name = 'class borders'
        # initialize remaining class attributes
        self.instrument_area = instrument_area
        self.time_interval = time_interval
        self.Aplawspeed = Aplawspeed
        self.Bplawspeed = Bplawspeed
        self.Aexpospeed = Aexpospeed
        self.Bexpospeed = Bexpospeed
        self.Cexpospeed = Cexpospeed
        self.ncells = ncells

    # this attribute is not initialized at time of object class creation
    # but only if requested
    @LazyProperty
    # Bulk Variable attribute: obtain bulk variables values using the "plaw" drop velocity
    def bulkvar_vplaw(self):
        return self.bulk_variables()

    # this attribute is not initialized at time of object class creation
    # but only if requested
    @LazyProperty
    # Bulk Variable attribute: obtain bulk variables values using the "exponential" drop velocity
    def bulkvar_vexpo(self):
        return self.bulk_variables(_speed_='expo')

    # this attribute is not initialized at time of object class creation
    # but only if requested
    @LazyProperty
    # phase space parameters attribute for flux representation:
    # obtain the satistical moments (phase space parameteres) of the flux pdf
    def psp_flux(self):
        return self.flux_phase_space_parameters()

    # this attribute is not initialized at time of object class creation
    # but only if requested
    @LazyProperty
    # phase space parameters attribute for cloud representation under the "plaw" drop velocity:
    # obtain the satistical moments (phase space parameteres) of the cloud-plaw pdf
    def psp_cloud_vplaw(self):
        return self.cloud_phase_space_parameters()

    # this attribute is not initialized at time of object class creation
    # but only if requested
    @LazyProperty
    # phase space parameters attribute for cloud representation under the "exponential" drop velocity:
    # obtain the satistical moments (phase space parameteres) of the cloud-expo pdf
    def psp_cloud_vexpo(self):
        return self.cloud_phase_space_parameters(_speed_='expo')

    # this attribute is not initialized at time of object class creation
    # but only if requested
    @LazyProperty
    # spectrum object for flux representation: return the probability density (flux pdf)
    def spectrum_flux(self):
        return self.flux_drop_pdf()

    # Moment Calculator method for the flux pdf
    def flux_moment_calculator(self, _alphalist_):
        """
        Purpose: calculate for each record all the moments in the list of moments _alphalist_ (flux pdf)

        Return: data frame with moments as columns

        _alphalist_: list of moments order (needs to be a python list)
        """

        moments_dict = dict()
        drops_per_record = self.data.sum(axis=1)
        class_span = self.classlimits.loc['right', :] - self.classlimits.loc['left', :]
        for elem in _alphalist_:
            class_span_alphap1 = pow(self.classlimits.loc['right', :], elem + 1) - pow(self.classlimits.loc['left', :], elem + 1)
            moments_dict[f"M{elem}"] = self.data.dot(class_span_alphap1 / (class_span * (elem + 1))) / drops_per_record

        _df_ = pd.DataFrame(moments_dict)
        return _df_

    # Moment Calculator method for the cloud pdf:
    # plaw" or "exponential" drop velocity are impelmenterd
    def cloud_moment_calculator(self, _alphalist_, _speed_='plaw'):
        """
        Purpose: calculate for each record all the moments in the list of moments _alphalist_ (cloud pdf)

        Return: data frame with moments as columns

        _alphalist_: list of moments order (needs to be a python list)
        _speed_: specifiy the law to use for the drop velocity deafult values is "plaw"
        """

        # the zero-th moment is necessary for calculation
        # we need to add to the list of calculated moments and then drop it if the zero-th
        #   moment is not explicitely requested
        zero_drop = False
        _alpha_ = _alphalist_.copy()

        # if _speed_ is anything but 'expo' the plaw formula for drop velocity is considered
        if (_speed_ == 'expo'):
            if 0 not in _alpha_:
                _alpha_.append(0)
                zero_drop = True
            moments_dict = dict()
            drops_per_record = self.data.sum(axis=1)
            class_df = self.classlimits
            _array_ = np.zeros(class_df.shape[1])
            for elem in _alpha_:
                for diamclass in range(0, class_df.shape[1]):
                    class_span = round(class_df.loc['right', class_df.columns[diamclass]] -
                                       class_df.loc['left', class_df.columns[diamclass]], 3)
                    _delta_ = class_span/self.ncells
                    _sum_ = 0.0
                    for j in range(0, self.ncells):
                        d_effective = class_df.loc['left', class_df.columns[diamclass]] + (_delta_/2) + j*_delta_
                        _sum_ = _sum_ + \
                            pow(d_effective, elem) / (self.Aexpospeed - self.Bexpospeed*np.exp(-self.Cexpospeed*d_effective))*_delta_
                    _array_[diamclass] = _sum_ / class_span
                _series_ = pd.Series(_array_, index=class_df.columns)
                moments_dict[f"X{elem}"] = self.data.dot(_series_) / drops_per_record

            #
            _df_ = pd.DataFrame(moments_dict)
            columns_drop = list()
            for elem in _alpha_:
                columns_drop.append(f"X{elem}")

            #
            if zero_drop is True:
                columns_drop.append('M0')
                for elem in _alpha_:
                    _df_[f"M{elem}"] = _df_[f"X{elem}"]/_df_.X0
            else:
                for elem in _alpha_:
                    if elem != 0:
                        _df_[f"M{elem}"] = _df_[f"X{elem}"]/_df_.X0
                    else:
                        _df_[f"M{elem}"] = _df_[f"X{elem}"]

            _df_.drop(columns=columns_drop, inplace=True)

        # the drop velocity is power law
        else:
            if 0 not in _alpha_:
                _alpha_.append(0)
                zero_drop = True
            _alpha_new = _alpha_.copy()
            i = 0
            for elem in _alpha_:
                _alpha_new[i] = _alpha_[i]-self.Bplawspeed
                i = i+1
            _df_ = self.flux_moment_calculator(_alpha_new)/self.Aplawspeed

            #
            columns_drop = list()
            for elem in _alpha_new:
                columns_drop.append(f"M{elem}")

            #
            if zero_drop is True:
                columns_drop.append('X-0.67')
                for elem in _alpha_new:
                    _df_[f"X{elem}"] = _df_[f"M{elem}"]/_df_['M-0.67']
            else:
                for elem in _alpha_new:
                    if elem != -0.67:
                        _df_[f"X{elem}"] = _df_[f"M{elem}"]/_df_['M-0.67']
                    else:
                        _df_[f"X{elem}"] = _df_[f"M{elem}"]

            _df_.drop(columns=columns_drop, inplace=True)

            #
            columns_name = list()
            if zero_drop is True:
                for elem in _alpha_:
                    if elem != 0:
                        columns_name.append(f"M{elem}")
            else:
                for elem in _alpha_:
                    columns_name.append(f"M{elem}")

            _df_.columns = columns_name

        return _df_

    # Method for remove isolted drop counts (outliers)
    def outlier_deletion(self, change_original=False):
        """
        Purporse: eliminate outliers counts. We look for the class with
        maximun count. From this class, we move to the left (right) until a
        class of null count is found or the minimum (maximum) diameter class
        is reached: this class is the left most (right most) class. All counts
        left of the leftmost and right of the right most class are set to zero.
        E.g. : the following 20 class disdrometer count
        [60 157 121 124 68 67 74 44 14 10 18 11 0 2 0 0 1 0 0 0] is reduced to
        [60 157 121 124 68 67 74 44 14 10 18 11 0 0 0 0 0 0 0 0]

        Return: an element of the disdrorain class (A or B).
        A) Change input disdrorain class element by eliminating outliers
        B) Create new element of disdrorain class equal to input except for
        records were outliers are eliminatedself.

        change_original: flag to decide if disdrorain class element output is (A <-> True) or (B <-> False)
        """

        # the clean up function finds
        # 1) the class with max count
        # 2) the class span of consencutive non zero counts which include the class with max count
        # All other counts will be considered outliers
        def cleanup(x):
            imax = x.argmax()
            zeroin = np.where(x == 0)[0]
            if len(np.where(zeroin > imax)[0]) > 0:
                rbi = np.where(zeroin > imax)[0].min()
                rb = zeroin[rbi]
                for i in range(rb, len(x)):
                    x[i] = 0
            if len(np.where(zeroin < imax)[0]) > 0:
                lbi = np.where(zeroin < imax)[0].max()
                lb = zeroin[lbi]
                for i in range(0, lb):
                    x[i] = 0
            return x

        # check if we keep original disdrorain class element untouched or not
        if change_original is False:
            # self_no_outliers = cp.deepcopy(self)
            rainobj = cp.deepcopy(self)
        else:
            rainobj = self

        sum_prior = rainobj.data.sum(axis=1)  # total number of drops in each record prior to outlier removal
        _matrix_ = rainobj.data.values
        for i, row in enumerate(_matrix_):
            _matrix_[i] = cleanup(row)
        _tempdata_ = pd.DataFrame(_matrix_, columns=rainobj.data.columns, index=rainobj.data.index)
        sum_after = _tempdata_.sum(axis=1)  # total number of drops in each racord after outlier removal

        # create summary data frame -- START
        _df1_ = pd.DataFrame({'ndrops_prior': sum_prior.values})
        _df2_ = pd.DataFrame({'ndrops_after': sum_after.values})
        summary = _df1_.join(_df2_)
        summary_final = cp.deepcopy(summary.loc[summary.ndrops_prior != summary.ndrops_after, :])
        summary_final.loc[:, 'delta_drops'] = summary_final.ndrops_after - summary_final.ndrops_prior
        summary_final.loc[:, 'perc_delta_drops'] = (summary_final.ndrops_after -
                                                    summary_final.ndrops_prior) / summary_final.ndrops_prior * 100
        summary_final.reset_index(inplace=True)
        summary_final.rename(columns={'index': 'record number'}, inplace=True)
        # create summary data frame -- STOP

        # return rainobj
        if change_original is False:
            return rainobj, summary_final
        else:
            return summary

    # Method for removing records with a too small drop count
    def remove_counts_below_threshold(self, _countth_=60):
        """
        Purporse: Remove records with a too small drop count

        Return: object of class disdrorain without the records with too small count

        _countth_: the count threshold
        """

        disdroclass_good = cp.deepcopy(self)
        disdroclass_good.data['ndrops'] = disdroclass_good.data.sum(axis=1)
        disdroclass_good.data = disdroclass_good.data[disdroclass_good.data['ndrops'] >= _countth_]
        disdroclass_good.data.reset_index(inplace=True)
        disdroclass_good.data.drop(columns=['ndrops', 'record number'], inplace=True)
        disdroclass_good.data.index.names = ['record number']

        return (disdroclass_good)

    # Method for removing "too narrow" counts: we keep records where at least _nclmin_ are occupied
    def remove_narrow(self, _nclmin_=4):
        """
        Purporse: Remove records where too few classes are occupied

        Return: XXXX

        _nclmin_: minumum number of class to be occupied
        """

        disdroclass_notnarrow = cp.deepcopy(self)
        disdroclass_narrow = cp.deepcopy(self)
        disdroclass_narrow.data['nzero_ncl'] = disdroclass_narrow.data[disdroclass_narrow.data != 0].count(axis=1)
        disdroclass_notnarrow.data['nzero_ncl'] = disdroclass_notnarrow.data[disdroclass_notnarrow.data != 0].count(axis=1)
        disdroclass_notnarrow.data = \
            disdroclass_notnarrow.data.loc[disdroclass_notnarrow.data.nzero_ncl >= _nclmin_, :].drop(columns=['nzero_ncl'])
        disdroclass_narrow.data = disdroclass_narrow.data.loc[disdroclass_narrow.data.nzero_ncl < _nclmin_, :].drop(columns=['nzero_ncl'])

        _row_ = [[disdroclass_notnarrow.data.shape[0], disdroclass_notnarrow.data.sum().sum(),
                  round(disdroclass_notnarrow.bulkvar_vplaw.R.sum()), 'legit']]
        _row_.append([disdroclass_narrow.data.shape[0], disdroclass_narrow.data.sum().sum(),
                      round(disdroclass_narrow.bulkvar_vplaw.R.sum()), 'narrow'])
        _summary_ = pd.DataFrame(_row_, columns=['nrecords', 'ndrops', 'rainfall_total', 'type'])

        disdroclass_notnarrow.data.reset_index(inplace=True)
        disdroclass_notnarrow.data.drop(columns=['record number'], inplace=True)
        disdroclass_notnarrow.data.index.names = ['record number']

        return (disdroclass_notnarrow, disdroclass_narrow, _summary_)

    # Method for calculating the phase space parameters (statistical moments) of the pdf in the flux representation
    def flux_phase_space_parameters(self):
        """
        Purpose: calculate for each record the following "flux" quantity:
            drop count (N), mean (mu), standard deviation (sigma), skewness (gamma), kurtosis (kappa),
            fifth central moment (eta), sixth central moment (omega)

        Return: data frame with all the above quantities as columns
        """
        list_moments_order = list([1, 2, 3, 4, 5, 6])
        _df_ = self.flux_moment_calculator(list_moments_order)
        _df_['N'] = self.data.sum(axis=1)
        _df_['mu'] = _df_.M1
        _df_['sigma'] = pow((_df_.M2 - pow(_df_.M1, 2)), 0.5)
        _df_['gamma'] = (_df_.M3 + (2 * pow(_df_.M1, 3)) - (3 * _df_.M1 * _df_.M2)) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 1.5))
        _df_['kappa'] = (_df_.M4 - (3 * pow(_df_.M1, 4)) + (6 * _df_.M2 * pow(_df_.M1, 2))
                         - (4 * _df_.M1 * _df_.M3)) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 2))
        _df_['eta'] = (_df_.M5 + (4 * pow(_df_.M1, 5)) + (10 * _df_.M3 * pow(_df_.M1, 2)) - (10 * _df_.M2 * pow(_df_.M1, 3))
                       - (5 * _df_.M4 * _df_.M1)) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 2.5))
        _df_['omega'] = (_df_.M6 - (6 * _df_.M5 * _df_.M1) + (15 * _df_.M4 * pow(_df_.M1, 2)) - (20 * _df_.M3 * pow(_df_.M1, 3))
                         + (15 * _df_.M2 * pow(_df_.M1, 4)) - (5 * pow(_df_.M1, 6))) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 3))
        _df_.drop(columns=['M1', 'M2', 'M3', 'M4', 'M5', 'M6'], inplace=True)
        return _df_

    # Method for calculating the phase space parameters (statistical moments) of the pdf in the cloud representation
    def cloud_phase_space_parameters(self, _speed_='plaw'):
        """
        Purpose: calculate for each record the following "cloud" quantity:
            drop count (Nv), mean (mu), standard deviation (sigma), skewness (gamma), kurtosis (kappa),
            fifth central moment (eta), sixth central moment (omega)

        Return: data frame with all the above quantities as columns

        _speed_: specifiy the law to use for the drop velocity deafult values is "plaw"
        """
        list_moments_order = list([0, 1, 2, 3, 4, 5, 6])
        if (_speed_ == 'expo'):
            _df_ = self.cloud_moment_calculator(list_moments_order, _speed_='expo')
        else:
            _df_ = self.cloud_moment_calculator(list_moments_order, _speed_='plaw')
        _df_['N'] = self.data.sum(axis=1)
        _df_['Nv'] = (round(1 / ((self.instrument_area / 1000000) * self.time_interval) * _df_.N * _df_.M0)).astype(int)
        _df_['mu'] = _df_.M1
        _df_['sigma'] = pow((_df_.M2 - pow(_df_.M1, 2)), 0.5)
        _df_['gamma'] = (_df_.M3 + (2 * pow(_df_.M1, 3)) - (3 * _df_.M1 * _df_.M2)) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 1.5))
        _df_['kappa'] = (_df_.M4 - (3 * pow(_df_.M1, 4)) + (6 * _df_.M2 * pow(_df_.M1, 2))
                         - (4 * _df_.M1 * _df_.M3)) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 2))
        _df_['eta'] = (_df_.M5 + (4 * pow(_df_.M1, 5)) + (10 * _df_.M3 * pow(_df_.M1, 2)) - (10 * _df_.M2 * pow(_df_.M1, 3))
                       - (5 * _df_.M4 * _df_.M1)) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 2.5))
        _df_['omega'] = (_df_.M6 - (6 * _df_.M5 * _df_.M1) + (15 * _df_.M4 * pow(_df_.M1, 2)) - (20 * _df_.M3 * pow(_df_.M1, 3))
                         + (15 * _df_.M2 * pow(_df_.M1, 4)) - (5 * pow(_df_.M1, 6))) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 3))
        _df_.drop(columns=['N', 'M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6'], inplace=True)
        return _df_

    # Method for calculating the bulk variables assuming a "plaw" or a "exponential" relation for the drop velocity
    def bulk_variables(self, _speed_='plaw'):
        """
        Purpose: calculate for each record the Number of drops (N), the number of drops per unit volume N_V,
        rainfall rate (R), reflectivity  (Z), liquid water content (W), rainfall rate per unit drop (R/N=r),
        reflectivity per unit drop (Z/N_V=z), liquid water content per unit drop (W/N_V=w)

        Return: data frame with all the above quantities as columns

        _speed_: specifiy the law to use for the drop velocity deafult values is "plaw"
        """

        PI = 3.141592653589793  # approximate value of pi
        seconds_in_hour = 3600  # for converting rainfall rate in mm/h

        # initialize data frame for bulk variables
        NNvbulk = pd.DataFrame(columns=['N', 'Nv', 'R', 'Z', 'W', 'r', 'z', 'w'])
        NNvbulk = NNvbulk.astype(
            dtype={'N': 'int64', 'Nv': 'float64', 'R': 'float64', 'Z': 'float64', 'W': 'float64', 'r': 'float64',
                   'z': 'float64', 'w': 'float64'})

        # number of drops trough the catchment area N
        NNvbulk['N'] = self.data.sum(axis=1)

        # Rainfall rate depends on the 3rd "flux" moment no matter what is the equation for v(D)
        R_df_ = self.flux_moment_calculator([3])
        NNvbulk['R'] = (PI / 6) * (1 / (self.instrument_area * self.time_interval)) * seconds_in_hour * NNvbulk.N * R_df_.M3

        list_moments_order = list([0, 3, 6])
        if (_speed_ == 'expo'):
            _df_ = self.cloud_moment_calculator(list_moments_order, _speed_='expo')
        else:
            _df_ = self.cloud_moment_calculator(list_moments_order)

        # the division by 1000000 is necessary to have catchment area in meter squared
        NNvbulk['Nv'] = (round(1 / ((self.instrument_area / 1000000) * self.time_interval) * NNvbulk.N * _df_.M0)).astype(int)
        NNvbulk['Nv_float'] = (1 / ((self.instrument_area / 1000000) * self.time_interval) * NNvbulk.N * _df_.M0)
        NNvbulk['Z'] = NNvbulk.Nv_float * _df_['M6']
        NNvbulk['W'] = NNvbulk.Nv_float * _df_['M3']
        NNvbulk['r'] = NNvbulk['R'] / NNvbulk['N']
        NNvbulk['z'] = NNvbulk['Z'] / NNvbulk['Nv_float']
        NNvbulk['w'] = NNvbulk['W'] / NNvbulk['Nv_float']

        NNvbulk.drop(columns=['Nv_float'], inplace=True)
        return NNvbulk

    # Method for calculating the spectrum (probability density) for each record (flux representation)
    def flux_drop_pdf(self):
        """
        Purpose: calculate the instantaneous drop pdf

        Return: data frame with pdf value per class (one row per record)
        """
        dropfreq = cp.deepcopy(self.data)
        dropfreq['ndrops'] = dropfreq.sum(axis=1)
        dropfreq = dropfreq.div(dropfreq['ndrops'].values, axis=0)
        dropfreq.drop(['ndrops'], axis=1, inplace=True)
        class_span = self.classlimits.loc['right', :] - self.classlimits.loc['left', :]
        droppdf = dropfreq/class_span

        return droppdf

    # Method for calculating additional paramters describing the pdf of drop diamters (flux representation)
    def flux_additional_pdf_paramters(self):
        """
        Purpose: Calculate "extra" parematers describing the pdf p(D)
                 See https://hess.copernicus.org/articles/16/329/2012/hess-16-329-2012.pdf
        Return: data frame with the values of additional parameters
        """
        pdf_df = self.flux_drop_pdf()
        _df_ = cp.deepcopy(self.flux_phase_space_parameters()[['N']])

        # mode variables
        idxmax = pdf_df.idxmax(axis=1)  # column -- class with the maximum value of the pdf
        pdfmax = pdf_df.max(axis=1)  # max value of pdf
        diammax = self.classlimits.loc['left', idxmax] + (self.classlimits.loc['right', idxmax] - self.classlimits.loc['left', idxmax]) / 2
        _df_.loc[:, 'class_mode'] = idxmax.values
        _df_.loc[:, 'modeclassn'] = (idxmax.str.extract(r'([0-9]+)')).values.astype(int)
        _df_.loc[:, 'pdf_mode'] = pdfmax
        _df_.loc[:, 'D_mode'] = diammax.values

        # span variables
        res = pdf_df[pdf_df != 0].stack()
        res.rename_axis(index=['record number', 'class'], inplace=True)
        res2 = res.reset_index().drop(columns=[0])
        firstclass = res2.groupby(['record number']).head(1)
        firstclass.set_index('record number', inplace=True)
        firstclass.columns = ['first_notzero_class']
        lastclass = res2.groupby(['record number']).tail(1)
        lastclass.set_index('record number', inplace=True)
        lastclass.columns = ['last_notzero_class']

        _df_ = _df_.join(firstclass)
        _df_ = _df_.join(lastclass)
        _df_['fclassn'] = _df_.first_notzero_class.str.extract(r'([0-9]+)').astype(int)
        _df_['lclassn'] = _df_.last_notzero_class.str.extract(r'([0-9]+)').astype(int)
        _df_['D_span'] = self.classlimits.loc['right', _df_.last_notzero_class].values - \
            self.classlimits.loc['left', _df_.first_notzero_class].values

        # gradient
        _df_.loc[:, 'pln'] = (_df_.modeclassn-4)  # .astype(str)
        _df_.loc[:, 'prn'] = (_df_.modeclassn+4)  # .astype(str)
        _df_.loc[_df_.modeclassn-4 < _df_.fclassn, 'pln'] = (_df_.fclassn)  # .astype(str)
        _df_.loc[_df_.modeclassn+4 > _df_.lclassn, 'prn'] = (_df_.lclassn)  # .astype(str)

        _df_.loc[:, 'pl'] = 'C' + (_df_.modeclassn-4).astype(str)
        _df_.loc[:, 'pr'] = 'C' + (_df_.modeclassn+4).astype(str)
        _df_.loc[_df_.modeclassn-4 < _df_.fclassn, 'pl'] = 'C' + (_df_.fclassn).astype(str)
        _df_.loc[_df_.modeclassn+4 > _df_.lclassn, 'pr'] = 'C' + (_df_.lclassn).astype(str)
        _df_.reset_index(inplace=True)

        pdf_df_matrix = pdf_df.values

        rn = _df_.index.values
        pdfpl = _df_.pln.values-1
        pdfpr = _df_.prn.values-1
        pdfmode = _df_.modeclassn.values-1

        _df_.loc[:, 'delta_pdf_left'] = pdf_df_matrix[rn, pdfpl] - pdf_df_matrix[rn, pdfmode]
        _df_.loc[:, 'delta_pdf_right'] = pdf_df_matrix[rn, pdfpr] - pdf_df_matrix[rn, pdfmode]
        _df_.loc[:, 'delta_D_left'] = self.classlimits.loc['left', _df_.class_mode].values + \
            (self.classlimits.loc['right', _df_.class_mode].values - self.classlimits.loc['left', _df_.class_mode].values)/2 - \
            self.classlimits.loc['left', _df_.pl].values
        _df_.loc[:, 'delta_D_right'] = self.classlimits.loc['right', _df_.pr].values - \
            (self.classlimits.loc['left', _df_.class_mode].values +
             (self.classlimits.loc['right', _df_.class_mode].values - self.classlimits.loc['left', _df_.class_mode].values)/2)
        _df_.loc[:, 'grad_left'] = _df_.delta_pdf_left / _df_.delta_D_left
        _df_.loc[:, 'grad_right'] = _df_.delta_pdf_right / _df_.delta_D_right

        _df_.set_index('record number', inplace=True)
        _df_.drop(columns=['class_mode', 'first_notzero_class', 'last_notzero_class', 'pl', 'pr', 'pln', 'prn', 'delta_pdf_left',
                           'delta_pdf_right', 'delta_D_left', 'delta_D_right'], inplace=True)
        _df_.columns = ['N', 'class_mode', 'pdf_mode', 'D_mode', 'first_notzero_class', 'last_notzero_class', 'D_span',
                        'grad_left', 'grad_right']

        return _df_

    # Method for calculating the renormalized spectrum (probability density) for each record (flux representation)
    def flux_renormalized_spectrum(self, _bin_=0.2, Dr_left=-5, Dr_right=25):
        """
        Purpose: Calculate the renormalized probability density function for each record (flux representation)
                 D -> Dr = (D-mu)/sigma

        Return: data frame with the renormalized spectrum for each record

        _bin_: bin size used for calcularing the renormalized spectrum
        Dr_left: minimum vale of Dr used when calculating the renormalized spectrum
        Dr_right: maximum vale of Dr used when calculating the renormalized spectrum
        """
        _spectrum_ = self.flux_drop_pdf()
        pdf = _spectrum_.values

        cll = cp.deepcopy(self.classlimits)

        _N_ = self.flux_phase_space_parameters()[['N']].values
        _mu_ = self.flux_phase_space_parameters()[['mu']].values
        _sigma_ = self.flux_phase_space_parameters()[['sigma']].values

        _b_ = np.append(cll.loc[['left']].values, cll.loc[['right']].values[0, cll.shape[1]-1])
        _borders_ = np.tile(_b_, (_spectrum_. shape[0], 1))
        _borders_renorm = (_borders_[:, :]-_mu_[:])/_sigma_[:]

        _nbins_ = np.floor((Dr_right - Dr_left)/_bin_).astype(int)
        res = np.zeros((_spectrum_.shape[0], _nbins_))

        # print(res.shape)
        for _row_ in range(0, _spectrum_.shape[0]):
            # print(_row_)
            for i in range(0, pdf.shape[1]):
                lb = (_borders_renorm[_row_][i]+5)/_bin_
                lb_int = lb.astype(int)
                lb_rem = lb-lb_int
                rb = (_borders_renorm[_row_][i+1]+5)/_bin_
                rb_int = rb.astype(int)
                rb_rem = rb-rb_int
                # print(i,lb,lb_int,lb_rem,rb,rb_int,rb_rem)
                if(pdf[_row_][i] > 0):
                    # print(_row_,i,lb,lb_int,lb_rem,rb,rb_int,rb_rem)
                    for k in range(lb_int+1, rb_int):
                        res[_row_][k] = pdf[_row_][i]*_sigma_[_row_][0]
                        # print(i,k,res[_row_][k])
                    res[_row_][lb_int] = res[_row_][lb_int] + (1-lb_rem)*pdf[_row_][i]*_sigma_[_row_][0]
                    res[_row_][rb_int] = res[_row_][rb_int] + rb_rem*pdf[_row_][i]*_sigma_[_row_][0]

            res[_row_][:] = res[_row_][:]*_N_[_row_][0]

        xr = np.arange(Dr_left, Dr_right, _bin_)
        renpdf = res.sum(axis=0)/_N_.sum()

        _pd_ = pd.DataFrame({'Dr': xr, 'pdf': renpdf})
        _pd_ = _pd_.loc[_pd_.pdf > 0, :]
        return _pd_


# class to process 2DVD disdrometer data
class disdrorain_2dvd(object):
    def __init__(self, datapath=None, dataframe=None, fieldsep=' ',
                 instrument_area=5000, time_interval=60):
        """
        didrorain class for 2DVD data
        # datapath = path to file
        # dataframe = name of data frame
        # fieldsep = string separating field data
        # instrument_area = catchement area of the 2dvd disrometer in mm^2.
            Default value is 10000 = 100cm^2
        # time_interval = time interval of disdrometer recording in seconds.
            Default value is 60 = 1 minute
        """

        if datapath is not None:
            self.data = pd.read_csv(datapath, sep=fieldsep, header=None)
        if dataframe is not None:
            self.data = dataframe.copy()
        self.data.rename(columns={0: 'timestamp', 1: 'diameter', 2: 'speed'}, inplace=True)
        self.data.index.name = 'record number'
        self.instrument_area = instrument_area
        self.time_interval = time_interval

    # this attribute is not initialized at time of object class creation
    # but only if requested
    @LazyProperty
    # Bulk Variable attribute: obtain bulk variables values using the "plaw" drop velocity
    def bulkvar(self):
        return self.bulk_variables()

    # this attribute is not initialized at time of object class creation
    # but only if requested
    @LazyProperty
    # phase space parameters attribute for flux representation:
    # obtain the satistical moments (phase space parameteres) of the flux pdf
    def psp_flux(self):
        return self.flux_phase_space_parameters()

    # this attribute is not initialized at time of object class creation
    # but only if requested
    @LazyProperty
    # phase space parameters attribute for cloud representation under the "plaw" drop velocity:
    # obtain the satistical moments (phase space parameteres) of the cloud-plaw pdf
    def psp_cloud(self):
        return self.cloud_phase_space_parameters()

    # Method for calculating the bulk variables
    def bulk_variables(self):
        """
        Purpose: calculate for each record the Number of drops (N), the number of drops per unit volume N_V,
        rainfall rate (R), reflectivity  (Z), liquid water content (W), rainfall rate per unit drop (R/N=r),
        reflectivity per unit drop (Z/N_V=z), liquid water content per unit drop (W/N_V=w)

        Return: data frame with all the above quantities as columns
        """

        PI = 3.141592653589793  # approximate value of pi
        seconds_in_hour = 3600  # for converting rainfall rate in mm/h

        # initialize data frame for bulk variables
        NNvbulk = pd.DataFrame(columns=['N', 'Nv', 'R', 'Z', 'W', 'r', 'z', 'w'])
        NNvbulk = NNvbulk.astype(
            dtype={'N': 'int64', 'Nv': 'float64', 'R': 'float64', 'Z': 'float64', 'W': 'float64', 'r': 'float64',
                   'z': 'float64', 'w': 'float64'})

        # number of drops trough the catchment area N
        drops_per_record = self.data.groupby('timestamp')['diameter'].count()
        NNvbulk['N'] = drops_per_record.values

        # Rainfall rate depends on the 3rd "flux" moment no matter what is the equation for v(D)
        R_df_ = self.flux_moment_calculator([3])
        NNvbulk['R'] = (PI / 6) * (1 / (self.instrument_area * self.time_interval)) * seconds_in_hour * NNvbulk.N * R_df_.M3

        list_moments_order = list([0, 3, 6])
        _df_ = self.cloud_moment_calculator(list_moments_order)

        # the division by 1000000 is necessary to have catchment area in meter squared
        NNvbulk['Nv'] = (round(1 / ((self.instrument_area / 1000000) * self.time_interval) * NNvbulk.N * _df_.M0)).astype(int)
        NNvbulk['Nv_float'] = (1 / ((self.instrument_area / 1000000) * self.time_interval) * NNvbulk.N * _df_.M0)
        NNvbulk['Z'] = NNvbulk.Nv_float * _df_['M6']
        NNvbulk['W'] = NNvbulk.Nv_float * _df_['M3']
        NNvbulk['r'] = NNvbulk['R'] / NNvbulk['N']
        NNvbulk['z'] = NNvbulk['Z'] / NNvbulk['Nv_float']
        NNvbulk['w'] = NNvbulk['W'] / NNvbulk['Nv_float']

        NNvbulk.drop(columns=['Nv_float'], inplace=True)
        return NNvbulk

    # Moment Calculator method for the flux pdf
    def flux_moment_calculator(self, _alphalist_):
        """
        Purpose: calculate for each record all the moments in the list of moments _alphalist_ (flux pdf)

        Return: data frame with moments as columns

        _alphalist_: list of moments order (needs to be a python list)
        """

        moments_dict = dict()
        drops_per_record = self.data.groupby('timestamp')['diameter'].count()
        _tempdf_ = self.data.copy()
        for elem in _alphalist_:
            _tempdf_.loc[:, 'd_alpha'] = pow(_tempdf_.diameter, elem)
            _sum_ = _tempdf_.groupby('timestamp')['d_alpha'].sum()
            _malpha_ = _sum_/drops_per_record
            moments_dict[f"M{elem}"] = _malpha_.values

        _df_ = pd.DataFrame(moments_dict)
        return _df_

    # Moment Calculator method for the cloud pdf
    def cloud_moment_calculator(self, _alphalist_):
        """
        Purpose: calculate for each record all the moments in the list of moments _alphalist_ (cloud pdf)

        Return: data frame with moments as columns

        _alphalist_: list of moments order (needs to be a python list)
        """
        moments_dict = dict()
        drops_per_record = self.data.groupby('timestamp')['diameter'].count()
        _tempdf_ = self.data.copy()

        if 0 not in _alphalist_:
            _alphalist_.append(0)

        for elem in _alphalist_:
            _tempdf_.loc[:, 'd_alpha'] = pow(_tempdf_.diameter, elem)/_tempdf_.speed
            _sum_ = _tempdf_.groupby('timestamp')['d_alpha'].sum()
            _malpha_ = _sum_/drops_per_record
            moments_dict[f"M{elem}"] = _malpha_.values

        _df_ = pd.DataFrame(moments_dict)

        columns_drop = list()
        columns_name = list()

        if 0 not in _alphalist_:
            for elem in _alphalist_:
                columns_drop.append(f"M{elem}")
                if elem != 0:
                    columns_name.append(f"M{elem}")
                    _df_[f"X{elem}"] = _df_[f"M{elem}"]/_df_['M0']
        else:
            for elem in _alphalist_:
                columns_name.append(f"M{elem}")
                if elem != 0:
                    columns_drop.append(f"M{elem}")
                    _df_[f"X{elem}"] = _df_[f"M{elem}"]/_df_['M0']

        _df_.drop(columns=columns_drop, inplace=True)
        _df_.columns = columns_name
        return _df_

    # Method for calculating the phase space parameters (statistical moments) of the pdf in the flux representation
    def flux_phase_space_parameters(self):
        """
        Purpose: calculate for each record the following "flux" quantity:
            drop count (N), mean (mu), standard deviation (sigma), skewness (gamma), kurtosis (kappa),
            fifth central moment (eta), sixth central moment (omega)

        Return: data frame with all the above quantities as columns
        """
        list_moments_order = list([1, 2, 3, 4, 5, 6])
        _df_ = self.flux_moment_calculator(list_moments_order)
        drops_per_record = self.data.groupby('timestamp')['diameter'].count()
        _df_['N'] = drops_per_record.values
        _df_['mu'] = _df_.M1
        _df_['sigma'] = pow((_df_.M2 - pow(_df_.M1, 2)), 0.5)
        _df_['gamma'] = (_df_.M3 + (2 * pow(_df_.M1, 3)) - (3 * _df_.M1 * _df_.M2)) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 1.5))
        _df_['kappa'] = (_df_.M4 - (3 * pow(_df_.M1, 4)) + (6 * _df_.M2 * pow(_df_.M1, 2))
                         - (4 * _df_.M1 * _df_.M3)) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 2))
        _df_['eta'] = (_df_.M5 + (4 * pow(_df_.M1, 5)) + (10 * _df_.M3 * pow(_df_.M1, 2)) - (10 * _df_.M2 * pow(_df_.M1, 3))
                       - (5 * _df_.M4 * _df_.M1)) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 2.5))
        _df_['omega'] = (_df_.M6 - (6 * _df_.M5 * _df_.M1) + (15 * _df_.M4 * pow(_df_.M1, 2)) - (20 * _df_.M3 * pow(_df_.M1, 3))
                         + (15 * _df_.M2 * pow(_df_.M1, 4)) - (5 * pow(_df_.M1, 6))) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 3))
        _df_.drop(columns=['M1', 'M2', 'M3', 'M4', 'M5', 'M6'], inplace=True)
        return _df_

    # Method for calculating the phase space parameters (statistical moments) of the pdf in the cloud representation
    def cloud_phase_space_parameters(self):
        """
        Purpose: calculate for each record the following "cloud" quantity:
            drop count (Nv), mean (mu), standard deviation (sigma), skewness (gamma), kurtosis (kappa),
            fifth central moment (eta), sixth central moment (omega)

        Return: data frame with all the above quantities as columns
        """
        list_moments_order = list([0, 1, 2, 3, 4, 5, 6])
        _df_ = self.cloud_moment_calculator(list_moments_order)
        drops_per_record = self.data.groupby('timestamp')['diameter'].count()
        _df_['N'] = drops_per_record.values
        _df_['Nv'] = (round(1 / ((self.instrument_area / 1000000) * self.time_interval) * _df_.N * _df_.M0)).astype(int)
        _df_['mu'] = _df_.M1
        _df_['sigma'] = pow((_df_.M2 - pow(_df_.M1, 2)), 0.5)
        _df_['gamma'] = (_df_.M3 + (2 * pow(_df_.M1, 3)) - (3 * _df_.M1 * _df_.M2)) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 1.5))
        _df_['kappa'] = (_df_.M4 - (3 * pow(_df_.M1, 4)) + (6 * _df_.M2 * pow(_df_.M1, 2))
                         - (4 * _df_.M1 * _df_.M3)) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 2))
        _df_['eta'] = (_df_.M5 + (4 * pow(_df_.M1, 5)) + (10 * _df_.M3 * pow(_df_.M1, 2)) - (10 * _df_.M2 * pow(_df_.M1, 3))
                       - (5 * _df_.M4 * _df_.M1)) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 2.5))
        _df_['omega'] = (_df_.M6 - (6 * _df_.M5 * _df_.M1) + (15 * _df_.M4 * pow(_df_.M1, 2)) - (20 * _df_.M3 * pow(_df_.M1, 3))
                         + (15 * _df_.M2 * pow(_df_.M1, 4)) - (5 * pow(_df_.M1, 6))) / (pow((_df_.M2 - (_df_.M1 * _df_.M1)), 3))
        _df_.drop(columns=['N', 'M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6'], inplace=True)
        return _df_

# Remove "too narrow" counts: we keep records where at least _nclmin_=4 classes where occupied
    def remove_counts_below_threshold(self, _countth_=60):
        """
        Purporse: Remove records with a too small drop count

        Return: object of class disdrorain_2dvd without the records with too small count

        _countth_: the count threshold
        """
        disdroclass_good = cp.deepcopy(self)
        _count_ = self.data.groupby('timestamp')['diameter'].count().to_frame()
        _count_.rename(columns=({'diameter': 'ndrops'}), inplace=True)
        disdroclass_good.data.set_index('timestamp', inplace=True)
        disdroclass_good.data = disdroclass_good.data.join(_count_)
        disdroclass_good.data = disdroclass_good.data[disdroclass_good.data['ndrops'] >= _countth_]
        disdroclass_good.data.reset_index(inplace=True)
        disdroclass_good.data.drop(columns=['ndrops'], inplace=True)
        disdroclass_good.data.index.names = ['record number']

        return (disdroclass_good)

    # Method for removing drops with off-bound speed
    def remove_drops_offbound_speed(self):
        """
        Purporse: Remove drops with "suspect" speed.
                  See Kruger, A., and W. F. Krajewski, 2002: "Two-Dimensional Video Disdrometer: A Description"
                  Journal of Atmospheric and Oceanic Technology,19 (5), 602–617
        Return: object of class disdrorain_2dvd without the drops with offbound speed
        """
        disdroclass_good = cp.deepcopy(self)
        disdroclass_good.data['lower_speed_bound'] = (9.65-10.3*np.exp(-0.6*disdroclass_good.data.diameter))*0.6
        disdroclass_good.data['upper_speed_bound'] = (9.65-10.3*np.exp(-0.6*disdroclass_good.data.diameter))*1.4
        disdroclass_good.data.loc[(disdroclass_good.data.speed <= disdroclass_good.data.upper_speed_bound)
                                  & (disdroclass_good.data.speed >= disdroclass_good.data.lower_speed_bound), 'keep_flag'] = '1'
        disdroclass_good.data.loc[(disdroclass_good.data.speed > disdroclass_good.data.upper_speed_bound)
                                  | (disdroclass_good.data.speed < disdroclass_good.data.lower_speed_bound), 'keep_flag'] = '0'

        disdroclass_good.data = disdroclass_good.data.loc[disdroclass_good.data.keep_flag == '1']
        disdroclass_good.data.drop(columns=['lower_speed_bound', 'upper_speed_bound', 'keep_flag'], inplace=True)

        return (disdroclass_good)

    def twodvd_to_parsivel(self):
        """
        Purporse: cast 2DVD data into PARSIVEL counts
        Return: data frame with PARSIVEL counts 
        """

        def findclass_alaparsivel(row):
            """
            Purporse: find class to which a drop of agaive ndiameter belongs
            Return: interger = diameter class 
            """
            if row['diameter'] <= 1.25:
                return int(row['diameter']/0.125)+1
            if ((row['diameter'] > 1.25) & (row['diameter'] <= 2.5)):
                return int((row['diameter']-1.25)/0.25)+11
            if ((row['diameter'] > 2.5) & (row['diameter'] <= 5)):
                return int((row['diameter']-2.5)/0.5)+16
            if ((row['diameter'] > 5) & (row['diameter'] <= 10)):
                return int((row['diameter']-5)/1)+21
            if ((row['diameter'] > 10) & (row['diameter'] <= 20)):
                return int((row['diameter']-10)/2)+26
            else:
                return int((row['diameter']-20)/3)+31

        disdroclass_good = cp.deepcopy(self)
        disdroclass_good.data['bin'] = disdroclass_good.data.apply(lambda row: findclass_alaparsivel(row), axis=1)
        # print(disdroclass_good.data.head(5))
        _bins_ = disdroclass_good.data.groupby(['timestamp', 'bin']).agg({'diameter': 'count'})
        _bins_ = _bins_.rename(columns={'diameter': 'bincount'})
        _mydata_ = pd.DataFrame(data=_bins_.index.get_level_values('timestamp').unique())
        # print(_mydata_)

        for i in range(0, 32):
            _mydata_[i] = 0

        for _minid_ in _bins_.index.get_level_values('timestamp').unique():
            _aa_ = _bins_.loc[_minid_]
            for _bin_ in _aa_.index.get_level_values('bin'):
                _cc_ = _aa_.loc[_bin_].values
                _mydata_.loc[_mydata_.timestamp == _minid_, (_bin_-1)] = _cc_[0]
        _mydata_.drop(columns=['timestamp'], inplace=True)

        return (_mydata_)
