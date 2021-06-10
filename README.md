# RNASeq_Data_Analysis_DESeq2_Test-3
Differential Gene Expression Analysis using DESeq2 tool to find the top 1000 differentially expressed genes with p value, adjusted p value and log2 fold change values. In addition, the identification of top 20 significant up and down regulated genes &amp; their functions are also carried out.
Logic and Justification

1.	Top 1000 DEGs:

        a.	The adjusted p value(q) column was sorted from largest to smallest.
        b.	The rows corresponding to q=NA were removed.
        c.	The adjusted p value(q) column was given a threshold of q<0.05.
        d.	The top 1000 rows were extracted to another csv file.
        e.	The top 1000 DEGs are hereby obtained according to log2foldchange, p-value and adjusted p value


2.	Top 20 significant upregulated and downregulated genes:

        a.	The adjusted p value(q) column was sorted from largest to smallest.
        b.	The rows corresponding to q=NA were removed.
        c.	The adjusted p value(q) column was given a threshold of q<0.05
        d.	Log2foldchange column was sorted from largest to lowest.
        e.	Top 20 rows will be the top 20 significant upregulated genes.
        f.	Last 20 rows will be the top 20 significant downregulated genes.

3.	Function of genes:

       
        a.	The function of the significant upregulated and downregulated genes
                was identified using pantherdb.org and attached as excel sheets respectively.

_______________________________________

