# PICIfinder

PICIfinder is a tool that provide utilities to detect PICI directly from the Gram-positive bacteria genome. It is mainly based on the conserved pattern features of PICI, and uses hmmsearch to compare the input genome sequences to find the fragments that match the PICI pattern. The detection algorithm consists of four parts:
1. Prediction of candidate fragments based on co-localization of gene modules.
2. Distinguishing false positives of PICI candidate fragments.
3. Range determination of PICI by random forest model.
4. Determination of the grade of PICI.

## Quick Start

### Download and install
git clone [here](https://github.com/ZhangAnYunLab-SCU/PICIfinder)

### Requirements
Before running PICIfinder, you need to install python packages and some tools:
```
Bio.SeqIO 
FastaIO 
sklearn.ensemble 
sklearn.model_selection 
prokka 
prodigal 
hmmsearch 
diamond
```

### PICIfinder Usage
```
Usage: PICI_finder.py [OPTIONS]

Options: 
-i, --seqfile PATH input fasta file [required] 
-t, --cpu INTEGER number of parallel PICI_run runs] 
-l, --length INTEGER length in basepairs to limit input sequences [default=1000, can increase but not decrease] 
-o, --out PATH path to deposit output folder and temporary files, will create if doesn't exist [default= working directory] 
--help Show this message and exit.
```
### Example
The sample files are available in `example_data`. Run them according to the following instructions.<br/>
`python PICIfinder.py -i SaPI1.fasta -l 3000 -t 1 -o outpath/SaPI1`
