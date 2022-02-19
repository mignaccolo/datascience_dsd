

# <font color=#ADFF2F>Analysis DRW data set</font>

- Location: Darwin, Australia
- Latitude: 12.45S
- Longitude: 130.83E
- instrument: Joss-Waldvogel RD69
- Resolution: 60 seconds
- Catchment Area: 5,000 mm$^{2}$
- Time interval: Julian 307 2005 --> Julian 041 2006 / no gaps





## <font color=#6495ED>From Raw data to Clean data</font>

Execute the following linux commands from "/datascience_dsd/disdrometer_data_raw/DRW"

1. unzip drop counts

   ```shell
unzip counts.zip;
   ```
   
2. concatenate daily counts

   ```shell
   cd counts/;
   cat dat_200* > drw_raw;
   mv drw_raw .././;
   cd ..;
   ```

3. clean data

   execute from "/datascience_dsd/disdrometer_data_raw/DRW"

   - set negative counts to zero and remove time stamp

     ```shell
     awk '{for (i=1;i<NF-1;i++) {if($i>=0) printf"%s ",$i; else printf"0 "}; if ($(NF-1)>=0) printf "%s\n",$(NF-1); else printf "0\n"}'  drw_raw > _temp_;
     ```

   - clean data

     ```shell
     python3 ../../python_code/clean_data.py  --disdrolimits celllimits_RD69_20cl_darwin_horiz --disdrodata _temp_ --disdroarea 5000 --disdroout drw_r1min;
     rm _temp_; # remove temp file
     rm -rf counts; # remove counts folder
     mv drw_r1min ../../disdrometer_data_clean/DRW/./ # move clean data set in clean data folder
     ```
   
4. Examine file "summary_drw_r1min" to see results of cleaning procedure

5. "drw_r1min", the clean data set, is in /datascience_dsd/disdrometer_data_clean/DRW





## <font color=#6495ED>Analysis</font>

1. <font color=#191970>__Calculate drop size distribution parameters (including renormalized values) for the flux representation  $(N,p(D))$__ </font>

   execute from "/datascience_dsd/disdrometer_data_clean/DRW"

   ```shell
   python3 ../../python_code/drop_size_distribution_parameters.py --disdrodatapath . --disdrocatalog ../../data_catalog.csv --disdroacronym drw_r1min --renorm ../../statistical_moments_renormalization_parameters.csv --renorm_type flux --output analysis_results/dsd_parameters_drw_r1min_flux_r
   ```

   ​	Output is a space separated file "dsd_parameters_drw_r1min_flux_r" containing following columns:

   - $N$, $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$, $\omega$, $\mu_r$, $\gamma_r$, $\sigma_r$, $\kappa_r$, $\eta_r$
   - $\omega$ is the 6-th central moment: not used in the analysis

   

2. <font color=#191970>__Calculate drop size distribution parameters (including renormalized values)  for cloudexpo representation $(N_{V},f(D))$__ </font>

   execute from "/datascience_dsd/disdrometer_data_clean/DRW"

   ```shell
   python3 ../../python_code/drop_size_distribution_parameters.py --disdrodatapath . --disdrocatalog ../../data_catalog.csv --disdroacronym drw_r1min --representation cloudexpo --renorm ../../statistical_moments_renormalization_parameters.csv  --renorm_type cloudexpo --output analysis_results/dsd_parameters_drw_r1min_cloudexpo_r
   ```

   Output is a space separated file "dsd_parameters_drw_r1min_cloudexpo_r" containing following columns:

   - $N_{V}$, $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$, $\omega$, $\mu_r$, $\gamma_r$, $\sigma_r$, $\kappa_r$, $\eta_r$
   - $\omega$ is the 6-th central moment: not used in the analysis

   

3. <font color=#191970>__Calculate drop size distribution parameters (including renormalized values)  for cloudplaw representation $(N_{V},f(D))$___ </font>

   execute from "/datascience_dsd/disdrometer_data_clean/DRW"

   ```shell
   python3 ../../python_code/drop_size_distribution_parameters.py --disdrodatapath . --disdrocatalog ../../data_catalog.csv --disdroacronym drw_r1min --representation cloudplaw --renorm ../../statistical_moments_renormalization_parameters.csv --renorm_type cloudplaw --output analysis_results/dsd_parameters_drw_r1min_cloudplaw_r
   ```

   Output is a space separated file "dsd_parameters_drw_r1min_cloudplaw_r" containing following columns:

   - $N_{V}$, $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$, $\omega$, $\mu_r$, $\gamma_r$, $\sigma_r$, $\kappa_r$, $\eta_r$
   - $\omega$ is the 6-th central moment: not used in the analysis

   

4. <font color=#191970>__Calculate the cross correlation of the drop size distribution parameters for all drop size distribution representations__ </font>

   execute from "/datascience_dsd/disdrometer_data_clean/DRW"

   ```shell
   python3 ../../python_code/xcorr_drop_size_distribution_parameters.py --disdrodatapath . --disdrocatalog ../../data_catalog.csv --disdroacronym drw_r1min --xcorroutmatrix analysis_results/xcorr_dsd_paramters_drw_r1min_matrix --xcorroutplot analysis_results/xcorr_dsd_paramters_drw_r1min_4plot
   ```

   Output is two space separated files "xcorr_dsd_paramters_drw_r1min_matrix" and "xcorr_dsd_paramters_drw_r1min_4plot"

   - xcorr_dsd_paramters_drw_r1min_matrix
     - For each representation: (flux, cloudplaw, and cloudexpo): cross correlation values arranged in a symmetric matrix for all possible couples within the 6-tuple ($N$, $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$) or  or ($N_V$, $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$)
   - xcorr_dsd_paramters_drw_r1min_4plot
     - as previous file, but cross correlation values arranged for create "bokeh" visualization 

   

5. <font color=#191970>__Calculate the rescaled mutual information of the drop size distribution parameters for all drop size distribution representations__ </font> (execute from "/datascience_dsd/disdrometer_data_clean/DRW")

   ```shell
   python3 ../../python_code/mutualinfo_drop_size_distribution_parameters.py --disdrodatapath . --disdrocatalog ../../data_catalog.csv --disdroacronym drw_r1min --rmioutmatrix analysis_results/rmi_dsd_paramters_drw_r1min_matrix --rmioutplot analysis_results/rmi_dsd_paramters_drw_r1min_4plot
   ```

   Output is two space separated files "rmi_dsd_paramters_drw_r1min_matrix" and "rmi_dsd_paramters_drw_r1min_4plot"

   - rmi_dsd_paramters_drw_r1min_matrix
     - For each representation: (flux, cloudplaw, and cloudexpo): rescaled mutual information values arranged in a symmetric matrix for all possible couples within the 6-tuple ($N$, $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$) or  or ($N_V$, $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$)
   - rmi_dsd_paramters_drw_r1min_4plot
     - as previous file, but rescaled mutual information values arranged for create "bokeh" visualization 

   

6. <font color=#191970>__Perform the PCA analysis of the drop size distribution parameters for all drop size distribution representations__ </font>

   execute from "/datascience_dsd/disdrometer_data_clean/DRW"

   ```shell
   python3  ../../python_code/pca_dropsize_distribution_parameters.py --disdrodatapath . --disdrocatalog ../../data_catalog.csv --disdroacronym drw_r1min --pcavarianceexp analysis_results/pca_dsd_paramters_drw_r1min_varexpl --pcacomponents analysis_results/pca_dsd_paramters_drw_r1min_compo
   ```

   Output is two space separated files "pca_dsd_paramters_drw_r1min_varexpl" and "pca_dsd_paramters_drw_r1min_compo"

   - pca_dsd_paramters_drw_r1min_varexpl
     - For each representation: (flux, cloudplaw, and cloudexpo): variance explained by each PCA component for the 6-tuple ($N$, $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$), varused=N\_pdf, and  for 5-tuple ( $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$), varused=pdf 
   - pca_dsd_paramters_drw_r1min_compo
     - For each representation: (flux, cloudplaw, and cloudexpo): PCA coefficients for the 6-tuple ($N$, $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$), varused=N\_pdf, and  for 5-tuple ( $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$), varused=pdf

   

7. <font color=#191970>__Perform the LAF fit  (flux representation / parameter $\sigma$)__ </font>

   execute from "/datascience_dsd/disdrometer_data_clean/DRW"

   ```shell
   python3 ../../python_code/adaptive_fitting.py --disdrodatapath . --disdrocatalog ../../data_catalog.csv --disdroacronym drw_r1min --pdfparameter sigma --representation flux --renorm_table ../../statistical_moments_renormalization_parameters.csv --renorm_type flux --radiusseq ~/work_vivo/pioggia/DATABASE/ETLed_data/_radius_  --occupancy 20 --output analysis_results/adaptfit_flux#drw_r1min#sigma
   ```

   Output is a space separated file "adaptfit_flux#drw_r1min#sigma": Fields are

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

    

   To perform the fit for other statistical moments ($\kappa$ and $\eta$ ) substitute the word "sigma" with "kappa" or "eta" in the "--pdfparameter" and "--output" option. To perform the fit for the cloudexpo of cloudplaw representation, substitute "cloudexpo" or "cloudplaw" in the "--representation" and "--renorm_type" option. 

   

8. <font color=#191970>__Calculate the in situ accuracy of the LAF fit (flux representation / parameter $\sigma$)__ </font>

   execute from "/datascience_dsd/disdrometer_data_clean/DRW"

   ```shell
   python3 ../../python_code/accuracy_adaptive_fitting_insitu.py --input analysis_results/adaptfit_flux#drw_r1min#sigma --output analysis_results/insitu_fitaccuracy_flux#drw_r1min#sigma --pdfparameter sigma
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

   execute from "/datascience_dsd/disdrometer_data_clean/DRW"

   ```shell
   python3 ../../python_code/compare_insitu_fitting_curves.py --fitsiteA analysis_results/adaptfit_flux#drw_r1min#sigma --fitsiteB ../ALE/analysis_results/adaptfit_flux#ale_2dvd#sigma  --pdfparameter sigma --output analysis_results/comparison_curves#drw_r1min#ale_2dvd#sigma
   ```

   Output are two space separated files:

   1. "comparison_curves#ale_2dvd#drw_r1min#sigma#couples"
   2. "comparison_curves#ale_2dvd#drw_r1min#sigma#both"

   Both files contains the L2RD, wL2RD, MRD, and wMRD coefficients of relative dispersion formatted in different ways (for plotting purposes).

   Note: it is not important which dataset is considered as "A" and which one is considered as "B".

   

10. <font color=#191970>__Calculate the rescaled mutual information of the drop size distribution moments for all drop size distribution representations__ </font> (execute from "/datascience_dsd/disdrometer_data_clean/DRW")

    ```shell
    python3 ../../python_code/mutualinfo_drop_size_distribution_moments.py --disdrodatapath . --disdrocatalog ../../data_catalog.csv --disdroacronym drw_r1min --rmioutmatrix analysis_results/rmi_dsd_moments_drw_r1min_matrix --rmioutplot analysis_results/rmi_dsd_moments_drw_r1min_4plot
    ```

    Output is two space separated files "rmi_dsd_moments_drw_r1min_matrix" and "rmi_dsd_moments_drw_r1min_4plot"

    - rmi_dsd_moments_drw_r1min_matrix
      - For each representation: (flux, cloudplaw, and cloudexpo): rescaled mutual information values arranged in a symmetric matrix for all possible couples within the 6-tuple ($N$, $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$) or  or ($N_V$, $\mu$, $\sigma$, $\gamma$, $\kappa$, $\eta$)
    - rmi_dsd_moments_drw_r1min_4plot
      - as previous file, but rescaled mutual information values arranged for create "bokeh" visualization 

11. 



