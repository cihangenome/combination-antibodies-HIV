# combination-antibodies-HIV
This repository contains the R scripts that were used to analyze the sequencing data derived from HIV patients who received anti-HIV antibodies. 

There are three patient groups, two of which have received anti_HIV antibodies over the course of several weeks, whereas the third one is the placebo group. The R scripts deposited here process the bulk TCR repertoire sequencing data from each patient at three time points in each group. The TCR repertoire analysis results are shown in Figure 5B of the manuscript titled "Combination anti-HIV antibodies provide sustained virologic remission". These figure has four panels for each patient group. 

Changes in the HIV-specific breadth and depth of CD8+ T cells of study participants are shown in the two upper panels. Violin plots show the Gaussian kernel probability density of the breadth/depth values over time. The mean and standard deviation values of the time point-specific distribution are shown as circles and vertical lines, respectively. Principal component analysis (PCA) of the changes in the TCR repertoire characteristics is shown (the two lower panels). Each ellipse shows the 95% confidence interval in the PCA space and the center of each ellipse is indicated by larger sized symbols that represent specific time points. Lower left panels depict PCA results with the frequencies of the HIV-specific clonotypes ranked among the top 25 with respect to their P-values associated with the pairwise comparisons between the three time points. Lower right panels depict PCA results with the gene usage profiles derived from the top TRBV-TRBJ gene pairs present in the above clonotypes. Principal component (PC) 1 and PC2 represent a lower-dimensional representation of the input data consisting of the frequencies of the HIV-specific clonotypes (lower left panel) and the usage levels of the TRBV-TRBJ gene pairs (lower right panel) for each patient group. The P values were determined using the Wilcoxon signed-rank test.
