# Manhattan Plot
An R function for colorful style of Manhattan plot of genome-wide association studies results, or any association results with chromosome coordinates and p-values.

```R
source("https://raw.githubusercontent.com/jianggl2000/Manhattan-Plot/refs/heads/main/manhattan.R")

#Simulate data for Manhattan plot
nchr=26
nsnps=1000
testdata=data.frame(
    SNP=sapply(1:(nchr*nsnps), function(x) paste("rs",x,sep='')),
    CHR=rep(1:nchr,each=nsnps), 
    BP=rep(1:nsnps,nchr), 
	P=runif(nchr*nsnps)
)

#display the header of the GWAS data
head(testdata)

#plot the manhattan figure
manhattan(testdata)
```
The header of the testdata. Make sure your data frame contains columns of CHR, BP, and P. 
<pre>
  SNP CHR BP          P
1 rs1   1  1 0.03589639
2 rs2   1  2 0.96449950
3 rs3   1  3 0.45762041
4 rs4   1  4 0.68747009
5 rs5   1  5 0.96663348
6 rs6   1  6 0.79437829
</pre>

The names of CHR is number 1 to 26, where 23 for chrX, 24 for chrY, 25 for chrom XY and 26 for chrom MT.

The figure,

![Manhattan plot](manhattan.png)
