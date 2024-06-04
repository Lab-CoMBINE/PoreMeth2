# PoreMeth2
PoreMeth2 is an R package for the identification of Differentially Methylated Regions from Nanopore methylation data (inferred by methcallers such as Nanopolish, DeepSignal, Dorado or Guppy) of paired samples and for their functional interpretation.

## Installation
### Install devtools from CRAN (Required)
	install.packages("devtools")

### install PoreMeth2
	options(timeout=9999999)

	devtools::install_github("Lab-CoMBINE/PoreMeth2")


## Features
The functions provided by PoreMeth2 allow to:

- Detect shifts in CpGs' methylation levels and segment such shifts to obtain precise DMRs.
- Predict DMRs' mechanism of origin (random or selection) by simultaneously segmenting entropy levels.
- Functionally interpret DMRs by annotating them on genic and regulatory elements such as CpG Islands, Transcription Factor Binding Sites and Enhancers.

## Usage

### 1. Data Preparation
In order to use PoreMeth2 to identify DMRs, methylation files from both samples need to be elaborated by a Perl script which is also provided with this package.
The script `parse_nanopolish_entropy.pl` takes as input read-level methylation calls  and uses them to obtain $\beta$ (methylation levels across reads) and entropy (a measure of the relative proportion of possible epialleles) values for each CpG site in both samples.

#### Input 
The input file for the Perl script should contain read-level methylation calls obtained with Nanopolish or Modkit (on Guppy/Dorado output data) and should be structured as the output of `nanopolish call-methylation` (see [documentation](https://nanopolish.readthedocs.io/en/latest/quickstart_call_methylation.html)) and `modkit extract` (see [documentation](https://github.com/nanoporetech/modkit/blob/master/book/src/intro_extract.md)). Note that The tool used to call methylation needs to be specified to the Perl parser as option. 

> Note that DeepSignal2 methylation calls are also suitable to this tool, but since it requires using Guppy and Tombo (which are deprecated) we don't provide methods to elaborate them. If desired they can however be adapted to nanopolish format 

#### Usage
The Perl script is used as follows:

		cat /my/methylation/calls.tsv | perl parse_nanopolish_entropy.pl > /my/parsed/reads.txt	

#### Output

The output from the parser is a tab separated file that follows this structure:

| column | description |
| --- | --- |
| chrom | chromosome | 
| position | genomic position of the CpG |
| entropy | entropy ($S$) value | 
| entropy_cov | number of reads used to calculate $S$ |
| beta | methylation level ($\beta$) of the CpG |
| beta_cov | number of reads used to calculate $\beta$ |

### 2. DMR identification

#### PoreMeth2DMR
The R function `PoreMeth2DMR` performs the double segmentation of $\beta$ and $S$ to identify DMRs and calculate $\Delta\beta$ and $\Delta S$ between test and control sample. 

	   dmrs  <-  PoreMeth2DMR(TableTest, TableControl, omega = 0.1, eta = 1e-5, FW = 3)	

Where:
-  `TableTest` and `TableControl` are the tables output by  `parse_nanopolish_entropy.pl`  on the Test and Control samples respectively. 
	> Note that we strongly suggest to filter input data based on beta_cov (a threshold of 10 is appropriate for 25-30x sequencing data, while it should be lowered for experiments with less coverage).
	
- `omega` is an optional parameter that modulates the relative weight between the experimental and the biological variance. When omega is close to 1, the biological variance is much larger than the experimental one, while for values of omega close to 0 the experimental noise gives the leading contribution to the total variance. We suggest to use `omega` in the range 0.1-0.5.
- `eta` is an optional parameter that represents the baseline probability the mean process (m_i) changes its value for the HSLM algorithm. Suggested values are inside $10^{-7}$-$10^{−3}$ range.
- `FW` is the minimum number of datapoints for a DMR to be called (DMRs made of a number of CpGs smaller than `FW` are discarded).

The output of `PoreMeth2DMR` is a data table containing the following fields:


| column | description |
| --- | --- |
| chr | chromosome |
| start | DMR's start position|
| end | DMR's end position |
| DeltaBeta | $\Delta\beta$ value between Test and Control samples |
| DeltaEntropy | $\Delta S$ value between Test and Control samples |
| BetaTest | mean $\beta$ value inside the DMR's coordinates for the Test sample|
| BetaControl | mean $\beta$ value inside the DMR's coordinates for the Control sample|
| EntropyTest | mean $S$ value inside the DMR's coordinates for the Test sample |
| EntropyControl | mean $S$ value inside the DMR's coordinates for the Control sample |
| NumCpG | number of CpG dinucleotides included in the region |
|  p | pvalue associated to the DMR | 


> Note that the function reports DMRs with any $\Delta\beta$ absolute value, but it is recommended to filter out DMRs with |$\Delta\beta$| < 0.2 to exclude unreliable results.

#### PoreMeth1DMR 
In case read-level methylation calls are not available and entropy cannot be calculated, it is possible to identify DMRs based strictly on $\beta$ shifts between samples by using the function `PoreMeth1DMR` from the previous version of this tool. 

The output of the function is identical to the one described for `PoreMeth2DMR`, without the fields related to Entropy. 

### 3. DMR interpretation
#### PoreMethAnnotate
DMRs obtained with `PoreMeth2DMR` can be annotated to genic and regulatory elements with the following command:
	
		PoreMethAnnotate(TableDMR, FileOut, NumProc = 5, AnnotationType = "Genes", Assembly = "hg19")	

Where:
- `TableDMR` is the output table of `PoreMeth2DMR`
- `NumProc` is an optional argument for the number of cores to use in parallel
- `AnnotationType` is an optional argument that specifies whether to annotate DMRs on genic elements only (`"Genes"`) or genic elements and regulatory features (`"GenesReg"`). 
- `Assembly` is an optional parameter that specifies the reference version to use for annotation (`"hg19"` or `"hg38"`).

The output file for `"Genes"` mode will contain the following additional fields:


| column | description |
| --- | ---|
| chr.GenCode | chromosome | 
| start.GenCode | start position of the annotated genic feature | 
| end.GenCode | end position of the annotated genic feature |
| feature.GenCode | genic feature ( either `Promoter`, `FirstExon`, `LastExon` or numbered `Exon`/`Intron` ) |
| strand.GenCode | either `+` or `-` |
| symbol.GenCode | gene name |
| type.GenCode | feature type ( e.g. `protein_coding`, `processed_transcript`, `lincRNA`... ) |
| chr.GenCode.overlap | chromosome |
| start.GenCode.overlap | start position of the overlap with the annotated genic feature |
| end.GenCode.overlap | end position of the overlap with the annotated genic feature |
| ratio1.GenCode.overlap | ratio of the overlap's length with respect to genic feature's length |
| ratio2.GenCode.overlap | ratio of the overlap's length with respect to DMR's length |
	

The output from "GenesReg" mode will also output, for CGIs, Enhancers, DNAse and TFBSs: 
- name
- chr
- start
- end
- chr overlap
- start overlap 
- end overlap 

### 4. Plots and statistics
#### PoreMeth2DMRStatistics
Given a DMR table (annotated or not) `PoreMeth2DMRStatistics` returns a summary the number of hyper/hypo-methylated regions with different values of $\Delta S$ across different genic features, CpG Islands and enhancers. 
The function can be used with:
	
		PoreMeth2DMRStatistics(TableDMR, Assembly = "hg19", BetaThr = 0.2, EntropyThr = 0.1, PValueThr = 0.05)

Where:
- `TableDMR` is the output table of `PoreMeth2DMR`
- `Assembly` is an optional parameter that specifies the reference version to use for statistics (`"hg19"` or `"hg38"`).
- `BetaThr` is the $\Delta \beta$ threshold applied for DMRs' classification 
- `EntropyThr` is the $\Delta S$ threshold applied for DMRs' classification
- `PValueThr` is the value threshold to consider a DMR


The output table presents the following rows for each category: 
	
| | | 
| ---| ---|
| HyperHigh | hyper-methylated ($\Delta \beta$ > `\|BetaThr\|`) and with high entropy ($\Delta S$ > `\|EntropyThr\|`)
| HyperMid | hyper-methylated ($\Delta \beta$ > `\|BetaThr\|`) and with mid entropy ( `-\|EntropyThr\|` < $\Delta S$ < `\|EntropyThr\|`)
| HyperLow | hyper-methylated ($\Delta \beta$ > `\|BetaThr\|`) and with low entropy ($\Delta S$ < `-\|EntropyThr\|`)
| HypoHigh | hypo-methylated ($\Delta \beta$ < `-\|BetaThr\|`) and with high entropy ($\Delta S$ > `\|EntropyThr\|`)
| HypoMid | hypo-methylated ($\Delta \beta$ < `-\|BetaThr\|`) and with mid entropy ( `-\|EntropyThr\|` < $\Delta S$ < `\|EntropyThr\|`)
| HypoLow | hypo-methylated ($\Delta \beta$ < `-\|BetaThr\|`) and with low entropy ($\Delta S$ < `-\|EntropyThr\|`)


#### PoreMeth2SingleExpQualityPlot
This function allows to automatically print plots for the evaluation of input data quality by displaying stats about `beta_cov` and `entropy_cov` (see Data Preparation) and $\beta$ and $S$ density functions.
Run with:
		
		PoreMeth2SingleExpQualityPlot(TableIn)

Where: 
- `TableIn` is the output from `parse_nanopolish_entropy.pl`.

#### PoreMeth2PairedExpQualityPlot
With this function it is possible to display stats about `beta_cov` and `entropy_cov` for common CpG dineucleotides in the Test/Control samples pair.

		PoreMeth2PairedExpQualityPlot(TableTest,TableControl)

- `TableTest` is the output from `parse_nanopolish_entropy.pl`.
- `TableControl` is the output from `parse_nanopolish_entropy.pl`.

#### PoreMeth2Plot
This function permits to plot DMRs $\Delta \beta$ and $\Delta S$ levels with genomic and regulatory annotations.
It can be used with:
		
		PoreMeth2Plot <- function(Input = NA, AnnotatedRes, PoreMeth2DMRResults, Meth1, Meth2)

Where:
- `Input` Are the coordinates to plot (chr:start-end) or gene symbol.
- `AnnotatedRes` Results from PoreMethAnnotate.
- `PoreMeth2DMRResults` Results from PoreMeth2DMR.
- `Meth1` Results from parse_nanopolish_entropy.
- `Meth2` Results from parse_nanopolish_entropy.

By default the resulting panel shows 3 plots a) the shifting levels of the DMRs, the $\Delta \beta$ of the CpGs on which the DMRs are calculated, b) the shifting levels of the Entropy, the $\Delta S$ of the CpGs on which the DMRs are calculated, c) the following genomic regions: the gene features (promoters, exons, introns and direction), the CGIs, the DNAse sites, the enhancers and the TFBS
Other the aesthetic parameters are available for a better readibility and understanding of the plots.
