# TopDownCrawl

Top-Down-Crawl is able to rapidly align quantitative binding data from experiments such as SELEX-seq and SMiLE-seq. Methods like these are able to provide highly reprodible k-mer level measurements of binding affinity, but understanding where an enriched k-mer falls along the bidning site depends on alignment. k-mer level data is more informative than simplistic PWM-based models which assume independent contributions between positions of the binding site.

The python script (TDC.py) accepts a table (.txt, .tsv, .csv, .xlsx, .xls) as input including the sequences in column 1, and the associated binding measurements in column 2. Values for reverse complements will be averaged if provided, or assumed to be equivalent if not. The first row is expected to contain header names of your choosing. An example input file is included (example_MEF2B_R1_E_k7.tsv).


### STEP 1: Install Required Packages

`pip install pandas matplotlib logomaker fonts font-roboto`

### STEP 2: Run Top-Down Crawl

`python TDC.py example_MEF2B_R1_E_k7.tsv`

### OUTPUT
```
<input>_aligned.tsv    Return aligned and reverse-complement-averaged input (tab-delimited)
<input>_summary.tsv    Show how many sequences align to each shift
<input>_PWM.png        Position weight matrix (PWM) of the aligned input, weighting each sequence by its associated binding metric, removing the 2 most extreme shifts 
<input>_PWM_rc.png     Reverse complement of the provided PWM
```
