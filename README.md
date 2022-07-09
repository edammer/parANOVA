# parANOVA(.dex) and plotVolc
Fast parallel ANOVA+Tukey statistics with 0 fallback option and volcano plot function

Compatible with the Seyfried Systems Biology Pipeline.

Sample pipeline input which benefits from fallback when Tukey p values are 0 (<1e-10) is provided as .RData ; optional 0 fallback is to the Bonferroni FDR on a two-sided unequal variance T test.

Outputs of the parANOVA.dex function include a data frame and CSV of one-way ANOVA statistics for overall model and Tukey-corrected p value for each pairwise comparison in the supplied abundance data.

A second plotVolc function is provided which benefits from the absence of 0 values when fallback p value estimation is enabled, and can leverage information from earlier in the pipeline about module membership (colors). Outputs of plotVolc include PDF and interactive (plotly) HTML plots which rely on ggplot2 package as well.

These functions generate variables expected by <a href="https://github.com/edammer/GOparallel/">the GOparallel gene ontology enrichment code</a>.
