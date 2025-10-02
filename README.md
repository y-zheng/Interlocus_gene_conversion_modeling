# How to Run the Simulation Script

This R script is a forward-time population-genetic simulator for a duplicated neutral locus with point mutation, reciprocal recombination, and interlocus gene conversion (IGC), plus optional selection on copy number. It requires R (with the **abind** package) and an initial population file. 

Run it from the command line as:

```bash
Rscript Neutral_locus_2copy.r.txt <sel> <gcr> <gcml> <rec.rate> <rep>
```

- **<sel>**: per-extra-copy selection coefficient (applied multiplicatively on diploids).  
- **<gcr>**: per-nucleotide IGC rate.  
- **<gcml>**: mean IGC track length (bp).  
- **<rec.rate>**: reciprocal recombination rate.  
- **<rep>**: replicate label used in input/output file names.  

The script expects an input file at:  
```
Neutral_Initial/Neutral_inieq_rep_<rep>.txt
```
This should contain 2N haplotypes (rows), each line starting with **g** followed by a binary string of alleles.  

Output files are written to `Round1.5/` and include:  
- `*_stats.txt`: summary statistics every 100 generations (θπ, θW, divergence, FST, etc.).  
- Population snapshots at the end of 200,000 generations (haplotype states per copy, SNP positions, and copy numbers).  
- If the duplicate is lost during simulation, a `*_failrecord_copy2.txt` file is written and the simulation restarts.  

---

# Generating the Initial Population File

The script requires an initial population file in the following format:  

- Path: `Neutral_Initial/Neutral_inieq_rep_<rep>.txt`  
- Exactly **2 × Ne** lines (here, 4000 if Ne=2000).  
- Each line: the character **g** immediately followed by a fixed-length string of **0/1** alleles representing the same set of segregating sites across haplotypes.  
  Example line:  
  ```
  g010101100110...
  ```

### Using `ms`

One common way to generate such a file is with Hudson’s **ms** simulator:

```bash
mkdir -p Neutral_Initial
ms 4000 1 -t 200 -s 500 |   awk '/^positions:/{ for(i=0;i<4000;i++){ getline; print "g"$0 } }'   > Neutral_Initial/Neutral_inieq_rep_1.txt
```

- `4000` = 2 × Ne haplotypes.  
- `-t` sets θ; `-s` forces the number of segregating sites (optional).  
- Prefixes each haplotype string with `g` for compatibility with the R script.  

### Alternatives

- Use `scrm` (ms-compatible).  
- Convert haplotypes from other simulators (e.g. SLiM, coalescent tools, or VCF output) into the same **g + 0/1 string** format with a custom script (Perl, Python, awk, etc.).  

**Requirements for the initial file:**  
1. All haplotypes must have the same number of sites.  
2. Sites must be binary (0/1).  
3. One haplotype per line, prefixed with `g`.  
