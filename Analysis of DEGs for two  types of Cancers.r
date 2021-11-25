Project Description:

You are required to conduct a study to analyze multi-platform molecular data for different types of
cancer. The cancer types are:
1. Lung Squamous Cell Carcinoma (LUSC),
2. Kidney Renal Clear Cell Carcinoma (KIRC).
The data (sent with this statement) for each cancer type:
1. Gene expression (GE),
◦ Two GE files for each cancer type:
1. “type-rsem-fpkm-tcga-t_paired.txt”, where type {lusc, kirc}: GE ∈ data for tissues with
cancer,
2. “type-rsem-fpkm-tcga_paired.txt”: GE data for tissues in a healthy case.
◦ Data are paired: each GE file will have the same number of cases (patients) and in the same
order.
◦ Files are tab-separated.
2. Copy number alterations (CNAs).
◦ A single file per cancer having the CNAs data for chromosome segments: arm-level and
focal-level copy numbers.
◦ Files are tab-separated.
◦ The set of cases may differ than that of the GE data.

Requirements:

1. Hypothesis Testing. For each cancer type, infer the differentially expressed genes (DEGs).
◦ Apply the appropriate test statistic for the following two cases:
1. Samples are paired,
2. Samples are independent.
◦ Report the set of DEGs in the above two pairing cases, and report how different these two
sets of genes.
▪ Does the difference increase with increasing the number of samples per cancer?
◦ Use the set of DEGs of the paired case and perform Gene Set Enrichment Analysis (GSEA)
on this set of genes (Bonus 5% per cancer type).
▪ Suggestion: you can use this GSEA Software.
2. Regression. For each cancer type, perform a multivariable regression analysis.
◦ The dependent variable is the gene expression for one gene, and the independent variables
are the CNAs.
◦ Pick up the five most differentially expressed genes in Requirement 1 and perform the
multivariable regression analysis on each one of these five genes separately.
Page 1
◦ The set of cases in the CNAs profiles and the GE profiles are not identical. You need to
make the regression analysis on the intersection of these cases.
◦ Report the significant copy number predictors for each gene.
• Support your findings/results/conclusions with figures.
• You have to deliver the following:
◦ All the code scripts you used for your analysis,
▪ Comments are a must.

◦ Project report:

▪ It should look like a research paper. It should have the following sections:
• Introduction,
• Methods: describe all the steps carefully and include all the used software packages,
• Results and Discussion: report your results in details and discuss them,
◦ You can augment your results with textual files or spreadsheets.
• Conclusion: list the overall findings of your analysis.
• Members Contribution: list in details what each member in your group did in
this project. Each member in the group may receive a different grade based on
the contribution weight.
