# STATISTICS OF THE GENOME ASSEMBLY
### Haptolina_ericina_var_UIO028.mainGenome.fasta

| Statistic    | Value       |
| ------------ | ----------- |
| total_length | 428,890,065 |
| scaffolds    | 388         |
| mean_length  | 1,105,387   |
| longest      | 6,770,620   |
| shortest     | 11,071      |
| N_count      | -           |
| Gaps         | -           |
| N50          | 1,772,137   |
| N50n         | 72          |
| N70          | 1,152,416   |
| N70n         | 130         |
| N90          | 574,375     |
| N90n         | 230         |

----
### Results from Repeatmodeler
**file name:** Haptolina_ericina_var_UIO028.mainGenome.fasta			
**sequences**:	388		
**total length:**  428,890,065 bp  (428890065 bp excl N/X-runs)			
**GC level:  **      	59.78%		
**bases masked:**  	 137,480,010 	32.05%	

| Type                                 | Number of Elements | Length Occupied (bp) | Percentage of Sequence |
| ------------------------------------ | ------------------ | -------------------- | ---------------------- |
| **Retroelements**                    | 15,465             | 23,570,249           | 5.50%                  |
| - SINEs:                             | 0                  | -                    | 0%                     |
| - Penelope                           | 0                  | -                    | 0%                     |
| - LINEs:                             | 6,869              | 7,965,013            | 1.86%                  |
| -- CRE/SLACS                         | 1,622              | 947,072              | 0.22%                  |
| -- L2/CR1/Rex                        | 0                  | -                    | 0%                     |
| -- R1/LOA/Jockey                     | 335                | 26,832               | 0.01%                  |
| -- R2/R4/NeSL                        | 0                  | -                    | 0%                     |
| -- RTE/Bov-B                         | 0                  | -                    | 0%                     |
| -- L1/CIN4                           | 2,333              | 3,880,544            | 0.90%                  |
| - LTR elements:                      | 8,596              | 15,605,236           | 3.64%                  |
| -- BEL/Pao                           | 0                  | -                    | 0%                     |
| -- Ty1/Copia                         | 5,079              | 7,471,044            | 1.74%                  |
| -- Gypsy/DIRS1                       | 2,238              | 7,059,670            | 1.65%                  |
| -- Retroviral                        | 300                | 133,771              | 0.03%                  |
| **DNA transposons**                  | 3,358              | 3,100,485            | 0.72%                  |
| - hobo-Activator                     | 236                | 25,016               | 0.01%                  |
| - Tc1-IS630-Pogo                     | 0                  | -                    | 0%                     |
| - En-Spm                             | 0                  | -                    | 0%                     |
| - MuDR-IS905                         | 0                  | -                    | 0%                     |
| - PiggyBac                           | 0                  | -                    | 0%                     |
| - Tourist/Harbinger                  | 983                | 506,109              | 0.12%                  |
| - Other (Mirage, P-element, Transib) | 0                  | -                    | 0%                     |
| Rolling-circles                      | 671                | 154,078              | 0.04%                  |
| Unclassified:                        | 214,792            | 90,652,630           | 21.14%                 |
| Total interspersed repeats:          | -                  | 117,323,364          | 27.36%                 |
| Small RNA:                           | 0                  | -                    | 0%                     |
| Satellites:                          | 1,762              | 1,393,428            | 0.32%                  |
| Simple repeats:                      | 178,984            | 16,533,747           | 3.86%                  |
| Low complexity:                      | 20,573             | 2,075,393            | 0.48%                  |
* most repeats fragmented by insertions or deletions			
  have been counted as one element			
RepeatMasker version 4.1.2-p1 , rushjob mode			         
run with rmblastn version 2.11.0+			
The query was compared to classified sequences in "HeUIO028_scaff3-families.fa"	

### Results from BRAKER2
Theses statisics are from the gtf file after running BRAKER2 on the genome.
(including simple bash commnad to get the stats)

The format of the gtf file is as follows:

| Scaffold     | Source    | Feature     | Start   | End     | Score | Strand | Frame | Attributes                             |
|--------------|-----------|-------------|---------|---------|-------|--------|-------|----------------------------------------|
| scaffold_103 | AUGUSTUS  | start_codon | 598011  | 598013  | .     | +      | 0     | transcript_id "g105827.t1"; gene_id "g105827"; |
| scaffold_103 | AUGUSTUS  | CDS         | 598011  | 598595  | 1     | +      | 0     | transcript_id "g105827.t1"; gene_id "g105827"; |
| scaffold_103 | AUGUSTUS  | exon        | 598011  | 598595  | .     | +      | .     | transcript_id "g105827.t1"; gene_id "g105827"; |
| scaffold_103 | AUGUSTUS  | gene        | 598011  | 598595  | 1     | +      | .     | g105827 |
| scaffold_103 | AUGUSTUS  | transcript  | 598011  | 598595  | 1     | +      | .     | g105827.t1 |
| scaffold_103 | AUGUSTUS  | stop_codon  | 598593  | 598595  | .     | +      | 0     | transcript_id "g105827.t1"; gene_id "g105827"; |
| scaffold_1   | AUGUSTUS  | stop_codon  | 2460597 | 2460599 | .     | -      | 0     | transcript_id "g71632.t2"; gene_id "g71632"; |
| scaffold_1   | AUGUSTUS  | CDS         | 2460597 | 2460879 | 0.72  | -      | 1     | transcript_id "g71632.t2"; gene_id "g71632"; |
| scaffold_1   | AUGUSTUS  | exon        | 2460597 | 2460879 | .     | -      | .     | transcript_id "g71632.t2"; gene_id "g71632"; |
| scaffold_1   | AUGUSTUS  | intron      | 2460880 | 2462310 | 0.49  | -      | .     | transcript_id "g71632.t2"; gene_id "g71632"; |
| scaffold_1   | AUGUSTUS  | CDS         | 2462311 | 2462603 | 0.68  | -      | 0     | transcript_id "g71632.t2"; gene_id "g71632"; |
| scaffold_1   | AUGUSTUS  | exon        | 2462311 | 2462603 | .     | -      | .     | transcript_id "g71632.t2"; gene_id "g71632"; |
| scaffold_1   | AUGUSTUS  | intron      | 2462604 | 2462664 | 0.55  | -      | .     | transcript_id "g71632.t2"; gene_id "g71632"; |
| scaffold_1   | AUGUSTUS  | CDS         | 2462665 | 2463132 | 0.55  | -      | 0     | transcript_id "g71632.t2"; gene_id "g71632"; |
| scaffold_1   | AUGUSTUS  | exon        | 2462665 | 2463132 | .     | -      | .     | transcript_id "g71632.t2"; gene_id "g71632"; |
| scaffold_1   | AUGUSTUS  | start_codon | 2463130 | 2463132 | .     | -      | 0     | transcript_id "g71632.t2"; gene_id "g71632"; |
| scaffold_1   | AUGUSTUS  | gene        | 2460597 | 2463666 | 0.78  | -      | .     | g71632 |
| scaffold_1   | AUGUSTUS  | transcript  | 2460597 | 2463666 | 0.37  | -      | .     | g71632.t1 |
_etc._


```bash
# Count the number of unique genes
awk '$3 == "gene"' your_output.gtf | cut -f9 | sort -u | wc -l
110010
 # Count the number of unique transcripts
awk '$3 == "transcript"' your_output.gtf | cut -f9 | sort -u | wc -l
113522
```
The average length is 1884.28 base pairs, and the length of the largest sequence is 157,716 base pairs(!?!).
**NB. the longest gene seems too long? Have to check in the mapping, if this is reasonable or not**



