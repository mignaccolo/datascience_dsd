# <font color=6495ED>datascience_dsd</font>
### <font color=87CEFA>__A data science approach to the rainfall drop size distribution__</font>

This repository contains the code necessary to analyze disdrometer data, and explore the "data scientist" parametrization of the rain drop size distribution. 

The data scientist approach to the rain drop size distribution is discussed in detail in "A worldwide data science investigation of rainfall" by M. Ignaccolo & C. DeMichele: [[manuscript]](https:\\duckduckgo.com).

In addition to the code 4 representative data sets (one for each type of disdrometer: 2dvd, Parsivel, RD69, and RD80) are included together with 4 walk through documents (in markdown and html format) that showcase how to perform the analysis  conducted in [[manuscript]](https:\\duckduckgo.com).

<font color=FFA500>__Note__</font>: All code is Python 3.8, and GNU GPLv3 licensed. Code has been tested on machine with Linux Mint 20 Cinnamon v:4.6.7, Linux Kernel v:5.4.0-97-generic.



### <font color=87CEFA>__Table of contents__</font>

##### <font color=FF4500>Folders</font>

- __disdrometer_data_raw__: 4 raw disdrometers data (part of the 166 datasets used in [[manuscript]](https:\\duckduckgo.com)) . One dataset for each type of disdrometer considered in.
  - ALE sub folder: 2DVD disdrometer data
  - BBY sub folder: Joss-Waldvogel RD80  disdrometer data
  - DRW sub folder: Joss-Waldvogel RD69  disdrometer data
  - PES sub folder: Parsivel disdrometer data 
- __disdrometer_data_clean__:  the cleaned=ready for analysis version of the ALE, BBY, DRW, PES datasets
- __python_code__: all the python programs necessary to implement the data science approach to rain drop size distributions.

<font color=FF4500>__Files__</font>

- __LICENSE.txt__: License file

- __README.md__: This very file in markdown format.

- __requirements.txt__: The list of python3 packages required to run all python code. Run the below command to install all the packages

  ```shell
  pip3 install -r requirements.txt
  ```

- __data_catalog.csv__: The "database" catalog of all the dataset used in [[manuscript]](https:\\duckduckgo.com). All python code used for the analysis expect the disdrometer datasets to be organized in a catalog. 

  - GEO_ID: Country or Earth's Region where data were collected 
  - DSN: Data Set Number (catalog does not need to be ordered by DSN)
  - ID1: 3/4 letters acronym indicating name of the instrument's location. E. g. 'DRW'  for Darwin, AUS
  - ID2: full acronym= filename associated with the dataset
  - INSTRUMENT: Type of disdrometer used to collect the data
  - CELLLIMITS: if applicable the filename of disdrometer class limits. 
    - = NONE for 2DVD disdormeter data
    - = standard (no actual file) for Joss-WaldVogel RD80 data 
  - AREA_INSTRUMENT: the catchment area of the instrument
  - TIME_RESOLUTION: time resolution of data
    - 2DVD are always binned into 60s time intervals for purposes of statistical analysis 
  - ORIGIN: Data Origin. E.g. NASA campaign data
  - ADD_INFO: Any other additional information. E.g. 'Composite' indicates data set is the composition of other database datasets 
  - NREC: After the cleaning procedure is applied, the number of records (TIME_RESOLUTION time intervals) in the dataset. 

- __statistical_moments_renormalization_parameters.csv__: The minimum and maximum values used to rescale the statistical moments in a MIN-MAX renormalization.

  - renorm_type: a keyword referencing a set of values. User is free to adopt any keyword and renormalization parameters. We adopted keywords that relate to the methodology used to derive proper renormalization parameters. 

    - flux $\rightarrow$ min and max value are obtained analyzing the 5% and 95% percentiles of the statistical moments of the flux probability density function of  $p(D)$: see online material for detailed explanation.

    - cloudexpo $\rightarrow$ min and max value are obtained analyzing the 5% and 95% percentiles of the statistical moments of the cloud probability density function of  $f(D)$ when drops of non-2DVD data have an "exponential" diameter-speed velocity relation: see online material for detailed explanation.

    - cloudplaw $\rightarrow$  min and max value are obtained analyzing the 5% and 95% percentiles of the statistical moments of the cloud probability density function of  $f(D)$ when drops of non-2DVD data have an "power law" diameter-speed velocity relation: see online material for detailed explanation.

    - renormalization equation is 
      $$
      x \rightarrow x_{r}=\frac{x-x_{min}}{x_{max}-x_{min}}
      $$
      

  - statistical_moment: one of the 5 statistical moments considered ($[\mu, \sigma, \gamma,\kappa,\eta]$)

  - xmin: minimum value of the parameter

  - xmax: maximum value of the parameter

- __Walkthrough_Analysis_ALE.html__: Step by step guided analysis reproducing the analysis in [[manuscript]](https:\\duckduckgo.com) for the ALE dataset. 

- __Walkthrough_Analysis_BBY.html__: Step by step guided analysis reproducing the analysis in [[manuscript]](https:\\duckduckgo.com) for the BBY dataset.

- __Walkthrough_Analysis_DRW.html__: Step by step guided analysis reproducing the analysis in [[manuscript]](https:\\duckduckgo.com) for the DRW dataset.

- __Walkthrough_Analysis_PES.html__: Step by step guided analysis reproducing the analysis in [[manuscript]](https:\\duckduckgo.com) for the PES dataset.



### <font color=87CEFA>__Acknowledgements__</font>

We wish to sincerely thank colleagues and institutions who have generously provided disdrometer data.

1. Dr. C. R. Williams, the National Oceanic and Atmospheric Administration (NOAA) and Physical Sciences Laboratory (PSL) for the public availability of the BBY, and DRW datasets.
2. Dr. Ali Tokay and National Aeronautics and Space Administration (NASA) for the public availability of the the  Global Precipitation Measurement (GPM) mission Ground Validation HYMEX campain: ALE and PES datasets.

Mutual Information analysis was done adopting the Local Non-uniformity Correction via the python package [NPEET\_LNC](https://github.com/BiuBiuBiLL/NPEET_LNC), while nearest neighbors search in the statistical moments space associated with a rain drop size distribution were done using the python package [GriSpy](https://grispy.readthedocs.io/en/latest/).

