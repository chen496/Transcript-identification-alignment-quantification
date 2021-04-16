cDNA-seq analysis:
--------

Transcript identification/alignment, quantification and differential gene expression analysis

Usually, we map the sequence to NCBI genome databases, like mm10 or CRch38. 
For this specific project, we need to map the mRNA sequence to NCBI CCDS database. 
There is a big difference for alignment against mouse genome databases vs. NCBI CCDS database because of the following reason. 
Our project is phage display of mouse proteins encoded by mouse cDNA derived from mouse tissue mRNA. 
One major problem of phage display cDNA library is the reading frame, in which many clones may be non-coding cDNA sequences, such as 5' and 3' UTR. 
All these non-coding region sequences should not be eliminated by data analysis because they do not express any proteins. 
Otherwise, false positives will be very high. Therefore, it is important to align the NGS data against mouse CCDS database at
https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi.     

Data:
NGS128 mRNA sequence data is from a normal group, NGS130 mRNA sequence data  is from an experiment group.

Quantification for comparative ligandomics to identify disease-associated targets in Dr. WeiLi's paper

![alt text](https://github.com/chen496/RNA-seq-anaysis/blob/28d493bc2e814932049a4c8ae263901b19a1d1ea/1.mRNA-seq%20alignment-Bowtie2/phase%20display.png
)

Overview
--------
Read mapping and transcript identification strategies:

![alt text](https://github.com/chen496/RNA-seq-anaysis/blob/2515339f5474c5059b8501804d8ed2ab4792458d/1.mRNA-seq%20alignment-Bowtie2/Read%20mapping%20and%20transcript%20identification%20strategies.png
)


The pipeline has three steps:

(1) Alignment
Bowtie2: Use an ungapped aligner Bowtie2 map reads to the reference transcriptome

(2) Transcript quantification
Algorithms that quantify expression from transcriptome mappings used in this project inlcuded
1) RSEM (RNA-Seq by Expectation Maximization)
2) Salmon
3) kallisto 

(3) Differential gene expression analysis
Deseq2




Results:
----------

The RSEM gives us two levels expression estimates in isoforms and genes, respectively. 
RSEM outputs the expected read count in each transcript and estimated the effective transcript's sequence length. 
Divide the expected read counts by the length of expected transcript's length of each gene and then multiply 1000, so that get RPK values.  
The tables only contain the genes with raw count sum in two samples(NGS128 and NGS130) are greater than 10. 
Salmon, Kallisto, and RSEM find 453, 977, and 333 genes with mapped reads. 
Different method do give different results, but three method largely agree on those genes with largest number of reads. 



References
----------

Experiment:

.. [1] Li, W., Pang, I.H., Pacheco, M.T.F. and Tian, H., 2018. Ligandomics: a paradigm shift in biological drug discovery. Drug discovery today, 23(3), pp.636-643.

.. [2] Caberoy, N.B., Zhou, Y., Jiang, X., Alvarado, G. and Li, W., 2010. Efficient identification of tubby‐binding proteins by an improved system of T7 phage display. Journal of Molecular Recognition: An Interdisciplinary Journal, 23(1), pp.74-83.

Algorithms:

.. [3] Conesa, A., Madrigal, P., Tarazona, S., Gomez-Cabrero, D., Cervera, A., McPherson, A., Szcześniak, M.W., Gaffney, D.J., Elo, L.L., Zhang, X. and Mortazavi, A., 2016. A survey of best practices for RNA-seq data analysis. Genome biology, 17(1), pp.1-19.

.. [4] Everaert, C., Luypaert, M., Maag, J.L., Cheng, Q.X., Dinger, M.E., Hellemans, J. and Mestdagh, P., 2017. Benchmarking of RNA-sequencing analysis workflows using whole-transcriptome RT-qPCR expression data. Scientific reports, 7(1), pp.1-11.

.. [5] Patro, R., Duggal, G., Love, M.I., Irizarry, R.A. and Kingsford, C., 2017. Salmon provides fast and bias-aware quantification of transcript expression. Nature methods, 14(4), pp.417-419.

.. [6] Bray, N.L., Pimentel, H., Melsted, P. and Pachter, L., 2016. Near-optimal probabilistic RNA-seq quantification. Nature biotechnology, 34(5), pp.525-527.

.. [7] Li, B. and Dewey, C.N., 2011. RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC bioinformatics, 12(1), pp.1-16.
