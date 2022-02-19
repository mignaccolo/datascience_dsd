# <font color=#ADFF2F>Analysis ALE data set</font>

- Location: Ales, France
- Latitude: 44.13N
- Longitude: 4.09E
- instrument: 2DVD disdrometer
- Resolution: 60 seconds (after integration)
- Catchment Area: 10,000 mm$^{2}$
- Time interval: 30SEP2012 to 23OCT2012 non contiguous





## <font color=#6495ED>From Raw data to Clean data</font>

Execute the following linux commands from "/datascience_dsd/disdrometer_data_raw/PES"

1. unzip drop counts

   ```shell
   unzip hymex_france.zip;
   cd hymex_france;
   ```

2. expand all archives

   ```shell
   ls *.tar.gz | awk '{print "tar -zxvf "$1}' | sh;
   ```

3. Extract drop data from each .txt file

   ```sh
   ls hymex_2dvd_sn35*_drops.txt | awk '{print "awk \047BEGIN{getline d;getline d}{print substr(FILENAME,17,8),\$1,\$2,\$4}\047 "$1 " > __sn35_"substr($1,17,8)}' | sh;
   ```

4. Concatenate all drop data files

   ```shell
   ls __sn35_2012* | sort | awk '{print "cat "$a " >> ale_raw"}'  | sh;
   ```

5. clean data

   - set negative counts to zero and remove time stamp

     ```shell
     awk '{printf"%s-%s %s %s\n",$1,substr($2,1,5),$3,$4}' ale_all > _temp_
     
     # mv _temp_ and pes_raw to parent directory
     mv _temp_ .././;
     mv ale_raw .././;
     # switch to parent directory
     cd ../;
     ```

   - prepare data for statistical analysis

     ```shell
     python3 ../../python_code/clean_2DVD_data.py --disdroparsivellimits celllimits_PARSIVEL --disdrodata _temp_ --disdroout ale_2dvd;
     rm _temp_; # remove temp files
     rm -rf hymex_italy_pescara; # remove unzipped folder
     mv ale_2dvd ../../disdrometer_data_clean/ALE/./ # move clean data set in clean data folder
     ```

6. Examine file "summary_ale_2dvd" to see results of cleaning procedure

7. "ale_2dvd", the clean data set, is in /datascience_dsd/disdrometer_data_clean/ALE





## <font color=#6495ED>Analysis</font>

1. <font color=#191970>__Calculate drop size distribution parameters (including renormalized values) for the flux representation  $(N,p(D))$__ </font> (

   execute from "/datascience_dsd/disdrometer_data_clean/ALE"

   ```shell
   python3 ../../python_code/drop_size_distribution_paramters.py --disdrodatapath . --disdrocatalog ../../data_catalog.csv --disdroacronym ale_2dvd --renorm ../../statistical_moments_renormalization_parameters.csv --renorm_type flux --output analysis_results/dsd_parameters_ale_2dvd_flux_r
   ```

   Output is a space separated file "dsd_parameters_ale_2dvd_flux_r" containing following columns:

   - $N$, $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$, $\omega$, $\mu_r$, $\gamma_r$, $\sigma_r$, $\kappa_r$, $\eta_r$
   - $\omega$ is the 6-th central moment: not used in the analysis
   - $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$ are normalized using the "flux ranges

   

2. <font color=#191970>__Calculate drop size distribution parameters   for cloud representation $(N_{V},f(D))$ (renormalized values will be calculated using the cloudexpo renormalization values)__ </font>

   execute from "/datascience_dsd/disdrometer_data_clean/ALE"

   ```shell
   python3 ../../python_code/drop_size_distribution_parameters.py --disdrodatapath . --disdrocatalog ../../data_catalog.csv --disdroacronym ale_2dvd --representation cloud2dvd --renorm ../../statistical_moments_renormalization_parameters.csv  --renorm_type cloudexpo --output analysis_results/dsd_parameters_ale_2dvd_cloud_rexpo
   ```

   Output is a space separated file "dsd_parameters_ale_2dvd_cloud_rexpo" containing following columns:

   - $N_{V}$, $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$, $\omega$, $\mu_r$, $\gamma_r$, $\sigma_r$, $\kappa_r$, $\eta_r$
   - $\omega$ is the 6-th central moment: not used in the analysis
   - $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$ are normalized using the "cloud expo" ranges

   

3. <font color=#191970>__Calculate drop size distribution parameters   for cloud representation $(N_{V},f(D))$ (renormalized values will be calculated using the cloudplaw renormalization values)__ </font>

   execute from "/datascience_dsd/disdrometer_data_clean/ALE"

   ```shell
   python3 ../../python_code/drop_size_distribution_parameters.py --disdrodatapath . --disdrocatalog ../../data_catalog.csv --disdroacronym ale_2dvd --representation cloud2dvd --renorm ../../statistical_moments_renormalization_parameters.csv  --renorm_type cloudplaw --output analysis_results/dsd_parameters_ale_2dvd_cloud_rplaw
   ```

   Output is a space separated CSV file "dsd_parameters_ale_2dvd_cloud_rplaw" containing following columns:

   - $N_{V}$, $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$, $\omega$, $\mu_r$, $\gamma_r$, $\sigma_r$, $\kappa_r$, $\eta_r$
   - $\omega$ is the 6-th central moment: not used in the analysis
   - $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$ are normalized using the "cloud plaw" ranges

   

4. <font color=#191970>__Calculate the cross correlation of the drop size distribution parameters for all drop size distribution representations__ </font> (execute from "/datascience_dsd/disdrometer_data_clean/ALE")

   ```shell
   python3 ../../python_code/xcorr_drop_size_distribution_parameters.py --disdrodatapath . --disdrocatalog ../../data_catalog.csv --disdroacronym ale_2dvd --xcorroutmatrix analysis_results/xcorr_dsd_paramters_ale_2dvd_matrix --xcorroutplot analysis_results/xcorr_dsd_paramters_ale_2dvd_4plot
   ```

   Output is two space separated files "xcorr_dsd_paramters_ale_2dvd_matrix" and "xcorr_dsd_paramters_ale_2dvd_4plot"

   - xcorr_dsd_paramters_drw_r1min_matrix
     - For each representation: (flux  and cloud): cross correlation values arranged in a symmetric matrix for all possible couples within the 6-tuple ($N$, $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$) or  or ($N_V$, $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$) 
   - xcorr_dsd_paramters_drw_r1min_4plot
     - as previous file, but cross correlation values arranged for create "bokeh" visualization 

   

5. <font color=#191970>__Calculate the rescaled mutual information of the drop size distribution parameters for all drop size distribution representations__ </font>

   execute from "/datascience_dsd/disdrometer_data_clean/ALE"

   ```shell
   python3 ../../python_code/mutualinfo_drop_size_distribution_parameters.py --disdrodatapath . --disdrocatalog ../../data_catalog.csv --disdroacronym ale_2dvd --rmioutmatrix analysis_results/rmi_dsd_paramters_ale_2dvd_matrix --rmioutplot analysis_results/rmi_dsd_paramters_ale_2dvd_4plot
   ```

   Output is two space separated files "rmi_dsd_paramters_ale_2dvd_matrix" and "rmi_dsd_paramters_ale_2dvd_4plot"

   - rmi_dsd_paramters_ale_2dvd_matrix
     - For each representation: (flux and cloud): rescaled mutual information values arranged in a symmetric matrix for all possible couples within the 6-tuple ($N$, $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$) or  or ($N_V$, $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$)
   - rmi_dsd_paramters_ale_2dvd_4plot
     - as previous file, but rescaled mutual information values arranged for create "bokeh" visualization 

   

6. <font color=#191970>__Perform the PCA analysis of the drop size distribution parameters for all drop size distribution representations__ </font>

   execute from "/datascience_dsd/disdrometer_data_clean/ALE"

   ```shell
   python3  ../../python_code/pca_dropsize_distribution_parameters.py --disdrodatapath . --disdrocatalog ../../data_catalog.csv --disdroacronym ale_2dvd --pcavarianceexp analysis_results/pca_dsd_paramters_ale_2dvd_varexpl --pcacomponents analysis_results/pca_dsd_paramters_ale_2dvd_compo
   ```

   Output is two space separated files "pca_dsd_paramters_ale_2dvd_varexpl" and "pca_dsd_paramters_ale_2dvd_compo"

   - pca_dsd_paramters_ale_2dvd_varexpl
     - For each representation: (flux, and cloud): variance explained by each PCS component for the 6-tuple ($N$, $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$), varused=N\_pdf, and  for 5-tuple ( $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$), varused=pdf 
   - pca_dsd_paramters_ale_2dvd_compo
     - For each representation: (flux and cloud): PCA coefficients for the 6-tuple ($N$, $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$), varused=N\_pdf, and  for 5-tuple ( $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$), varused=pdf

   

7. <font color=#191970>__Perform the LAF fit (flux representation / parameter $\sigma$)__ </font>

   execute from "/datascience_dsd/disdrometer_data_clean/ALE"

   ```shell
   python3 ../../python_code/adaptive_fitting.py --disdrodatapath . --disdrocatalog ../../data_catalog.csv --disdroacronym ale_2dvd --pdfparameter sigma --representation flux --renorm_table ../../statistical_moments_renormalization_parameters.csv --renorm_type flux --radiusseq ~/work_vivo/pioggia/DATABASE/ETLed_data/_radius_  --occupancy 20 --output analysis_results/adaptfit_flux#ale_2dvd#sigma
   ```

   Output is the file "adaptfit_flux#ale_2dvd#sigma": Fields are

   - __mu_r__: the value of the renormalized mean $\mu$ 
   - __gamma_r__: the value of the renormalized skewness $\gamma$ 
   - __predicted_r__: the value of the predicted renormalized parameter ($\sigma_{r}$ in this case): the predicted value is the median of all rescaled values of the parameter ($\sigma_{r}$ in this case) inside a ball of radius "__radius__" centered in $(\mu_{r}, \gamma_{r})$
   - __radius__: the radius of the ball centered in $(\mu_{r}, \gamma_{r})$
   - __predicted_5%__: the 5% percentile of all the renormalized values of the parameter ($\sigma$ in this case) inside a ball of radius "__radius__" centered in $(\mu_{r}, \gamma_{r})$ 
   - __predicted_q1__: the 25% percentile of all the renormalized values of the parameter ($\sigma$ in this case) inside a ball of radius "__radius__" centered in $(\mu_{r}, \gamma_{r})$
   - __predicted_q3__: the 75% percentile of all the renormalized values of the parameter ($\sigma$ in this case) inside a ball of radius "__radius__" centered in $(\mu_{r}, \gamma_{r})$
   - __predicted_95%__: the 95% percentile of all the renormalized values of the parameter ($\sigma$ in this case) inside a ball of radius "__radius__" centered in $(\mu_{r}, \gamma_{r})$
   - __mu__: the value of $\mu$ corresponding the renormalized value $\mu_{r}$
   - __gamma__: the value of $\gamma$ corresponding the renormalized value $\gamma_{r}$
   - __sigma_median__: the value of the predicted parameter ($\sigma$ in this case): this is the un-renormalized version of __predicted_r__
   - __sigma_5%__: the 5% percentile of all values of the parameter ($\sigma$ in this case) inside a ball of radius "__radius__" centered in $(\mu, \gamma)$: this is the un-renormalized version of __predicted_5%__
   - __sigma_q1__: the 25% percentile of all values of the parameter ($\sigma$ in this case) inside a ball of radius "__radius__" centered in $(\mu, \gamma)$: this is the un-renormalized version of __predicted_q1%__
   - __sigma_q3__: the 75% percentile of all values of the parameter ($\sigma$ in this case) inside a ball of radius "__radius__" centered in $(\mu, \gamma)$: this is the un-renormalized version of __predicted_q3%__
   - __sigma_95%__: the 95% percentile of all values of the parameter ($\sigma$ in this case) inside a ball of radius "__radius__" centered in $(\mu, \gamma)$: this is the un-renormalized version of __predicted_95%__

    To perform the fit for other statistical moments ($\kappa$ and $\eta$ ) substitute the word "sigma" with "kappa" or "eta" in the "--pdfparameter" and "--output" option. To perform the fit for the cloud2dvd  representation, substitute "cloud2dvd" in the "--representation" and use "cloudexpo" or "cloudplaw"  in the "--renorm_type" option.

   

8. <font color=#191970>__Calculate the in situ accuracy of the LAF fit (flux representation / parameter $\sigma$)__ </font>

   execute from "/datascience_dsd/disdrometer_data_clean/ALE"

   ```shell
   python3 ../../python_code/accuracy_adaptive_fitting_insitu.py --input analysis_results/adaptfit_flux#ale_2dvd#sigma --output analysis_results/insitu_fitaccuracy_flux#ale_2dvd#sigma --pdfparameter sigma
   ```

   Output is the space separated file "insitu_fitaccuracy_flux#drw_r1min#sigma": Fields are

   - __analysis_type__: we either calculate the probability [pdf] associated with a given __value__ of the coeeficient of variation, or the cumulative probability [cdf]
   - __variable__: [coeffvar] indicates the coefficient of variation (relative dispersion) relative to the semi interquartile range. [coeffvar595] indicates the coefficient of variation (relative dispersion) relative to the semi [5%,95%] percentile range
   - __value__: the value of [coeffvar] or [coefffvar595]
   - __prob__: the value of the probability associated to __value__ 
   - __cumprob__: the value of the cumulative probability associate to __value__. __cumprob__ is missing [NaN] when __analysis_type__ is [pdf]
   - __avg_radius__: the average radius of the ball used for estimate the parameter (statistical moments) for a given __value__ of [coeffvar] or [coefffvar595]  
   - __site__: the site for which the analysis is performed. This information is contained in the input file

   

9. <font color=#191970>__Compare LAF fitted curves at different sites (flux representation / parameter $\sigma$)__ </font>

   execute from "/datascience_dsd/disdrometer_data_clean/ALE"

   ```shell
   python3 ../../python_code/compare_insitu_fitting_curves.py --fitsiteA analysis_results/adaptfit_flux#ale_2dvd#sigma --fitsiteB ../DRW/analysis_results/adaptfit_flux#drw_r1min#sigma  --pdfparameter sigma --output analysis_results/comparison_curves#ale_2dvd#drw_r1min#sigma
   ```

   Output are two space separated files:

   1. "comparison_curves#ale_2dvd#drw_r1min#sigma#couples"
   2. "comparison_curves#ale_2dvd#drw_r1min#sigma#both"

   Both files contains the L2RD, wL2RD, MRD, and wMRD coefficients of relative dispersion formatted in different ways (for plotting purposes).

   Note: it is not important which dataset is considered as "A" and which one is considered as "B".
   
   
   
10. <font color=#191970>__Calculate the rescaled mutual information of the drop size distribution moments for all drop size distribution representations__ </font> (execute from "/datascience_dsd/disdrometer_data_clean/ALE")

    ```shell
    python3 ../../python_code/mutualinfo_drop_size_distribution_moments.py --disdrodatapath . --disdrocatalog ../../data_catalog.csv --disdroacronym ale_2dvd --rmioutmatrix analysis_results/rmi_dsd_moments_ale_2dvd_matrix --rmioutplot analysis_results/rmi_dsd_moments_ale_2dvd_4plot
    ```

    Output is two space separated files "rmi_dsd_moments_ale_2dvd_matrix" and "rmi_dsd_moments_ale_2dvd_4plot"

    - rmi_dsd_moments_ale_2dvd_matrix
      - For each representation: (flux, cloudplaw, and cloudexpo): rescaled mutual information values arranged in a symmetric matrix for all possible couples within the 6-tuple ($N$, $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$) or  or ($N_V$, $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$)
    - rmi_dsd_moments_ale_2dvd_4plot
      - as previous file, but rescaled mutual information values arranged for create "bokeh" visualization 

