# ANOVA, Volcano, and DEx Barplot of Coexpression Module Members
Fast parallel ANOVA+Tukey statistics with 0 fallback option, volcano plot, and module-aware stacked barplot

Compatible with the Seyfried Systems Biology Pipeline.

Sample pipeline input which benefits from fallback p value calculation when unreliable ANOVA+Tukey p values compute to <1e-10 is provided as .RData
The fallback is to the Bonferroni FDR on a two-sided unequal variance T test, which approximates the Tukey calculation when significance is higher and p values are smaller.

Outputs of the parANOVA.dex function include a data frame and CSV of one-way ANOVA statistics for overall model and Tukey-corrected p value for each pairwise comparison in the supplied abundance data.

A second plotVolc function is provided which benefits from the fallback p value estimation, and can leverage information from earlier in the pipeline about module membership (colors). Outputs of plotVolc include PDF and interactive (plotly) HTML plots which rely on ggplot2 package as well.

These functions generate data stored in variables expected by <a href="https://github.com/edammer/GOparallel/">the GOparallel gene ontology enrichment code</a>.

A third function, DEXpercentStacked, outputs the fraction of each module occupied by DEx proteins for each pairwise comparison, using a heat scale to also show the mean log2(fold change) for the up and down proteins in each module, meeting the thresholding criteria set during use of the plotVolc function.
_______________

See notes the .R wrapper for these functions, and use therein of sample RData with pipeline inputs for the above functions available from this repository.
Sample outputs the wrapper is expected to generate are provided in <a href="https://github.com/edammer/parANOVA/raw/main/parANOVA.dex-sampleOutput.zip">this .ZIP file</a>.
