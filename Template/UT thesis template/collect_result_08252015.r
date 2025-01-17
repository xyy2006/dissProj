1, 5000 Type I and power we have for CV simulated method
# 1, alpha table is the same as I stored before
  #  Test_Simulation_result_collect_andPlot_forMoreTestPlayers.r -f Result_analyzeWith_independenceOn_RVsimulatedU_scoUpdated_simulated_data_AR1_simuDatasets5000_nulls_h0_totalSampleSize -o alpha_ind_5000replicates -S SNOW -j F
\begin{table}[ht]
\centering
\begin{tabular}{rrrrrrrrrrrrrrrrrrrr}
  \hline
 & pSSU & pSSUw & pScore & pSum & pUminP & geescoreP & UminP & aveP & aveP\_weighted & KM\_P & geescor
eP\_withDiagSigma & spu\_Inf\_P & spu\_weighted\_Inf\_P & aspu\_P & aspu\_weighted\_P & aspu\_sco\_P &
 aspu\_weighted\_sco\_P & aspu\_aspuw.sco\_P & aspu\_aspuw\_P \\
  \hline
independenceOn\_RVsimulatedU\_scoUpdated\_AR1\_h0\_500 & 0.047 & 0.048 & 0.047 & 0.052 & 0.048 & 0.047
 & 0.053 & 0.052 & 0.051 & 0.058 & 0.059 & 0.054 & 0.053 & 0.060 & 0.058 & 0.058 & 0.055 & 0.058 & 0.0
60 \\
  independenceOn\_RVsimulatedU\_scoUpdated\_AR1\_h0\_1000 & 0.047 & 0.046 & 0.044 & 0.058 & 0.048 & 0.
044 & 0.051 & 0.057 & 0.057 & 0.059 & 0.056 & 0.051 & 0.051 & 0.057 & 0.057 & 0.057 & 0.054 & 0.057 &
0.057 \\
  independenceOn\_RVsimulatedU\_scoUpdated\_AR1\_h0\_2000 & 0.049 & 0.047 & 0.051 & 0.048 & 0.048 & 0.
051 & 0.052 & 0.048 & 0.049 & 0.058 & 0.055 & 0.052 & 0.052 & 0.061 & 0.058 & 0.058 & 0.058 & 0.058 &
0.058 \\
  independenceOn\_RVsimulatedU\_scoUpdated\_AR1\_h0\_3000 & 0.051 & 0.052 & 0.052 & 0.052 & 0.050 & 0.
052 & 0.053 & 0.051 & 0.052 & 0.061 & 0.060 & 0.054 & 0.053 & 0.063 & 0.060 & 0.059 & 0.057 & 0.059 &
0.062 \\
   \hline
\end{tabular}
\end{table}
# 2, power as well, the plot.
  # Test_Simulation_result_collect_andPlot_forMoreTestPlayers.r -f Result_analyzeWith_independenceOn_RVsimulatedU_scoUpdated_simulated_data_AR1_simuDatasets5000_signals_h0.001_totalSampleSize -o power_ind_5000replicates -S SNOW -j F

#==========================================================================================================#  

2, 5000 Type I and power we have for CV permuted method, only has sample size 1000 and 3000
# 1, alpha table, new table
  # Test_Simulation_result_collect_andPlot_forMoreTestPlayers.r -f Result_analyzeWith_independenceOn_RVpermutedU_scoUpdated_simulated_data_AR1_simuDatasets5000_nulls_h0_totalSampleSize -o alpha_ind_5000replicates_permuted -S SNOW -j F
\begin{table}[ht]
\centering
\begin{tabular}{rrrrrrrrrrrrrrrrrrrr}
  \hline
 & pSSU & pSSUw & pScore & pSum & pUminP & geescoreP & UminP & aveP & aveP\_weighted & KM\_P & geescor
eP\_withDiagSigma & spu\_Inf\_P & spu\_weighted\_Inf\_P & aspu\_P & aspu\_weighted\_P & aspu\_sco\_P &
 aspu\_weighted\_sco\_P & aspu\_aspuw.sco\_P & aspu\_aspuw\_P \\
  \hline
independenceOn\_RVpermutedU\_scoUpdated\_AR1\_h0\_1000 & 0.047 & 0.046 & 0.044 & 0.058 & 0.048 & 0.044
 & 0.049 & 0.057 & 0.058 & 0.050 & 0.050 & 0.048 & 0.049 & 0.052 & 0.050 & 0.053 & 0.051 & 0.053 & 0.0
52 \\
  independenceOn\_RVpermutedU\_scoUpdated\_AR1\_h0\_3000 & 0.051 & 0.052 & 0.052 & 0.052 & 0.050 & 0.0
52 & 0.051 & 0.051 & 0.050 & 0.053 & 0.053 & 0.053 & 0.051 & 0.054 & 0.051 & 0.052 & 0.052 & 0.051 & 0
.054 \\
   \hline
\end{tabular}
\end{table}
# 2, power plot fails, due to no 500 and 2000.

#==========================================================================================================#  deal with:
  # Result_longiPointIncreaseTest_1timepoints_analyzeWith_independenceOn_RVsimulatedU_scoUpdated_simulated_data_AR1_simuDatasets5000_signals_h0.001_totalSampleSize1000_commonVariant
# enter R.
setwd("/work/02040/yyang/dissertation_proj/SimulationData") 
results <- foreach (i = 1:4) %do%
{
  filename <- sprintf("Result_longiPointIncreaseTest_%dtimepoints_analyzeWith_independenceOn_RVsimulatedU_scoUpdated_simulated_data_AR1_simuDatasets5000_signals_h0.001_totalSampleSize1000_commonVariant.RDS", i)
  result <- readRDS(filename)

}
names(results) <- sprintf("%d_measurements", 1:4)
save(results, file = "powerData_for_Result_longiPointIncreaseTest_1to4timepoints_analyzeWith_independenceOn_RVsimulatedU_scoUpdated_simulated_data_AR1_simuDatasets5000_signals_h0.001_totalSampleSize1000_commonVariant.Rdata")
# conclusion: because the variance is not tuned correctly, the power increase is trivial.

#==========================================================================================================#  deal with:
  # Result_longiPointIncreaseTest_1timepoints_analyzeWith_independenceOn_RVsimulatedU_scoUpdated_simulated_data_AR1_simuDatasets5000_signals_h0.001_totalSampleSize3000_commonVariant_withFixedTimeEffect
# enter R.
setwd("/work/02040/yyang/dissertation_proj/SimulationData") 
results <- foreach (i = 1:4) %do%
{
  filename <- sprintf("Result_longiPointIncreaseTest_%dtimepoints_analyzeWith_independenceOn_RVsimulatedU_scoUpdated_simulated_data_AR1_simuDatasets5000_signals_h0.001_totalSampleSize3000_commonVariant_withFixedTimeEffect.RDS", i)
  result <- readRDS(filename)

}
names(results) <- sprintf("%d_measurements", 1:4)
# conclusion: lack file. At that time, I didnt' finish it. Besides, the power increase is also very minor.


#==========================================================================================================#  deal with:
  #Result_longiPointIncreaseTest_1timepoints_analyzeWith_independenceOn_RVsimulatedU_scoUpdated_simulated_data_AR1_simuDatasets2000_signals_h0.001_totalSampleSize500_commonVariant_withFixedTimeEffect.RDS
  #Result_longiPointIncreaseTest_1timepoints_analyzeWith_independenceOn_RVsimulatedU_scoUpdated_simulated_data_AR1_simuDatasets2000_signals_h0.001_totalSampleSize2000_commonVariant_withFixedTimeEffect.RDS
setwd("/work/02040/yyang/dissertation_proj/SimulationData") 
results <- foreach (i = 1:4) %do%
{
  filename <- sprintf("Result_longiPointIncreaseTest_%dtimepoints_analyzeWith_independenceOn_RVsimulatedU_scoUpdated_simulated_data_AR1_simuDatasets2000_signals_h0.001_totalSampleSize3000_commonVariant_withFixedTimeEffect.RDS", i)
  result <- readRDS(filename)

}
names(results) <- sprintf("%d_measurements", 1:4)
# conclusion: 500 sample is minor increase. 2000 better, 3000 minor.

#========================================================================================================#
###::: re-simulate with true case in ARIC analyses and use LMM simulation method.skeleton
# then analyze the data
# -----------type I error--------------------.
# > temp
                   # pSSU pSSUw pScore  pSum pUminP geescoreP UminP  aveP aveP_weighted  KM_P geescoreP_withDiagSigma
# ar1_1000          0.056 0.052  0.063 0.081  0.048     0.063 0.053 0.082         0.078 0.062                   0.063
# ar1_2000          0.067 0.070  0.053 0.085  0.057     0.053 0.059 0.083         0.078 0.077                   0.081
# ar1_3000          0.049 0.051  0.043 0.073  0.044     0.043 0.045 0.071         0.066 0.057                   0.062
# ar1_500           0.073 0.072  0.049 0.087  0.060     0.049 0.064 0.088         0.082 0.081                   0.086
# exchangeable_1000 0.043 0.038  0.053 0.058  0.041     0.053 0.048 0.058         0.053 0.050                   0.048
# exchangeable_2000 0.054 0.051  0.046 0.040  0.050     0.046 0.055 0.042         0.040 0.059                   0.067
# exchangeable_3000 0.042 0.039  0.048 0.045  0.047     0.048 0.048 0.044         0.044 0.047                   0.049
# exchangeable_500  0.060 0.055  0.045 0.056  0.049     0.045 0.052 0.053         0.056 0.070                   0.064
# independence_1000 0.043 0.038  0.053 0.058  0.041     0.053 0.047 0.057         0.052 0.050                   0.047
# independence_2000 0.054 0.051  0.046 0.040  0.050     0.046 0.056 0.045         0.041 0.059                   0.066
# independence_3000 0.042 0.039  0.048 0.045  0.047     0.048 0.047 0.044         0.044 0.047                   0.048
# independence_500  0.060 0.055  0.045 0.056  0.049     0.045 0.052 0.055         0.055 0.070                   0.062
# unstructured_1000 0.042 0.034  0.055 0.056  0.043     0.055 0.046 0.056         0.052 0.047                   0.051
# unstructured_2000 0.054 0.053  0.049 0.041  0.048     0.049 0.054 0.043         0.041 0.060                   0.067
# unstructured_3000 0.039 0.039  0.048 0.046  0.047     0.048 0.049 0.046         0.042 0.047                   0.050
# unstructured_500  0.062 0.057  0.042 0.059  0.050     0.042 0.055 0.059         0.059 0.071                   0.065
                  # spu_Inf_P spu_weighted_Inf_P aspu_P aspu_weighted_P aspu_sco_P aspu_weighted_sco_P aspu_aspuw.sco_P
# ar1_1000              0.051              0.053  0.074           0.072      0.074               0.065            0.073
# ar1_2000              0.062              0.059  0.088           0.083      0.084               0.085            0.085
# ar1_3000              0.062              0.045  0.059           0.059      0.058               0.058            0.059
# ar1_500               0.066              0.064  0.089           0.083      0.081               0.076            0.079
# exchangeable_1000     0.047              0.048  0.053           0.050      0.053               0.055            0.052
# exchangeable_2000     0.055              0.055  0.064           0.064      0.062               0.061            0.064
# exchangeable_3000     0.052              0.048  0.059           0.057      0.054               0.050            0.051
# exchangeable_500      0.058              0.052  0.069           0.065      0.063               0.062            0.060
# independence_1000     0.047              0.047  0.051           0.050      0.052               0.050            0.049
# independence_2000     0.053              0.056  0.065           0.064      0.061               0.065            0.065
# independence_3000     0.053              0.047  0.058           0.056      0.053               0.047            0.051
# independence_500      0.063              0.052  0.066           0.059      0.063               0.059            0.060
# unstructured_1000     0.049              0.046  0.049           0.044      0.053               0.049            0.051
# unstructured_2000     0.059              0.054  0.062           0.062      0.064               0.064            0.064
# unstructured_3000     0.052              0.049  0.057           0.054      0.049               0.050            0.050
# unstructured_500      0.058              0.055  0.071           0.066      0.066               0.060            0.064
                  # aspu_aspuw_P
# ar1_1000                 0.072
# ar1_2000                 0.086
# ar1_3000                 0.059
# ar1_500                  0.086
# exchangeable_1000        0.050
# exchangeable_2000        0.064
# exchangeable_3000        0.057
# exchangeable_500         0.065
# independence_1000        0.047
# independence_2000        0.065
# independence_3000        0.057
# independence_500         0.061
# unstructured_1000        0.047
# unstructured_2000        0.067
# unstructured_3000        0.056
# unstructured_500         0.069

# -----------power of wk:ind varying sample sizes--------------------.
###::: cmd used
# Test_Simulation_result_collect_andPlot_forMoreTestPlayers.r -f Result_longiPointIncreaseTest_4timepoints_analyzeWith_independenceOn_RVsimulatedU_scoUpdated_simulated_data_AR1_simuDatasets1000_signals_h0.001_totalSampleSize\\d+?_commonVariant_withFixedTimeEffectInLMM.RDS -o PowerCurve_independenceWkCor_on_AR1data_h0.001_CVs_simulatedInLMM -S MC -j F
# Test_Simulation_result_collect_andPlot_forMoreTestPlayers.r -f Result_longiPointIncreaseTest_4timepoints_analyzeWith_independenceOn_RVsimulatedU_scoUpdated_simulated_data_AR1_simuDatasets1000_signals_h0.001_totalSampleSize\\d+?_commonVariant_withFixedTimeEffectInLMM.RDS -o PowerCurve_reducedTests_independenceWkCor_on_AR1data_h0.001_CVs_simulatedInLMM -S MC -j F
###::: Conclusion is the graphs are similar to before simulated with correlated residual model.

#-------------------power gain from increase the longi points of LMM simulation data---------------------------#
#"Result_longiPointIncreaseTest_\\d+?timepoints_analyzeWith_independenceOn_RVsimulatedU_scoUpdated_simulated_data_AR1_simuDatasets1000_signals_h0.001_totalSampleSize\\d+?_commonVariant_withFixedTimeEffectInLMM"
###::: cmd used
# Test_Simulation_result_collect_andPlot_forMoreTestPlayers_longiIncreasePower.r -f Result_longiPointIncreaseTest_\\d+?timepoints_analyzeWith_independenceOn_RVsimulatedU_scoUpdated_simulated_data_AR1_simuDatasets1000_signals_h0.001_totalSampleSize\\d+?_commonVariant_withFixedTimeEffectInLMM -S MC -j F -o PowerIncreaseCurve_independenceWkCor_on_AR1data_h0.001_CVs_simulatedInLMM
###:::  conclusion: power gain is very trivial. I plan to use the correlated residual simulation model previously generated.

#----------------------------------power gain from increase the longi points of correlated residual simulation data-------------------#
###::: cmd used:
# Test_Simulation_result_collect_andPlot_forMoreTestPlayers_longiIncreasePower.r -f Result_longiPointIncreaseTest_\\d+?timepoints_analyzeWith_independenceOn_RVsimulatedU_scoUpdated_simulated_data_AR1_simuDatasets2000_signals_h0.001_totalSampleSize\\d+?_commonVariant_withFixedTimeEffect\\.RDS -S MC -j F -o PowerIncreaseCurve_independenceWkCor_on_AR1data_h0.001_CVs_simulatedInCorrelatedResidual
###::: Conclusion: power gain is satisfactory.

#-------------------------------analyze the df_sig_collection to compare performance of LaSPU with aSPU on baseline. and compare with SSU------------#
c("APOC3","ANGPTL8","C19orf80","PAFAH1B2","ANGPTL4","LIPG","LPL","FAM65A") %in% df_sig_collected$aspu_P$df_gw_marginal_sig$Gene
subset(df_sig_collected$pSSU$df_gw_marginal_sig, Gene == "PAFAH1B2")  # use this in the df with all genes to fetch the result even not marginally sig (p < 0.001).