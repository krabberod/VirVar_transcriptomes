I renamed the gtf file to PkV-RF01.gtf
Changed the gene_id from "1," to "PkV-RF01_g1" "2," to "PkV-RF01_g2" etc. with this sed command:
```bash
sed -i.bak '/^#/! s/gene_id "\([0-9]*\),"/gene_id "PkV-RF01_g\1"/g;
     /^#/! s/transcript_id "\([0-9]*\),"/transcript_id "PkV-RF01_g\1"/g;
     /^#/! s/ID "\([0-9]*\),"/ID "PkV-RF01_g\1"/g' PkV-RF01_final.gtf
```

Original file:
```
##gtf-version 3
##source-version GeneMark.hmm_PROKARYOTIC 3.25
##date: Wed May  2 12:14:54 2018
# Sequence file name: PkV-RF01_final.fasta
# Model file name: PkV-RF01_hmm_combined.mod
# RBS: true
# Model information: GeneMarkS_gcode_1
PkV-RF01	GeneMark.hmm	gene	734	1762	-1485.666949	-	0	gene_id "1,"; ID "nbis-gene-1";
PkV-RF01	GeneMark.hmm	transcript	734	1762	-1485.666949	-	0	gene_id "1,"; transcript_id "nbisL2-cds-1"; ID "nbisL2-cds-1"; Parent "nbis-gene-1"; original_biotype "mrna";
PkV-RF01	GeneMark.hmm	exon	734	1762	-1485.666949	-	.	gene_id "1,"; transcript_id "nbisL2-cds-1"; ID "nbis-exon-1"; Parent "nbisL2-cds-1";
PkV-RF01	GeneMark.hmm	CDS	734	1762	-1485.666949	-	0	gene_id "1,"; transcript_id "nbisL2-cds-1"; ID "cds-1"; Parent "nbisL2-cds-1";
PkV-RF01	GeneMark.hmm	gene	3862	4602	-1075.440697	-	0	gene_id "2,"; ID "nbis-gene-2";
PkV-RF01	GeneMark.hmm	transcript	3862	4602	-1075.440697	-	0	gene_id "2,"; transcript_id "nbisL2-cds-2"; ID "nbisL2-cds-2"; Parent "nbis-gene-2"; original_biotype "mrna";
PkV-RF01	GeneMark.hmm	exon	3862	4602	-1075.440697	-	.	gene_id "2,"; transcript_id "nbisL2-cds-2"; ID "nbis-exon-234"; Parent "nbisL2-cds-2";
PkV-RF01	GeneMark.hmm	CDS	3862	4602	-1075.440697	-	0	gene_id "2,"; transcript_id "nbisL2-cds-2"; ID "cds-2"; Parent "nbisL2-cds-2";
```
Renamed file:
```
##gtf-version 3
##source-version GeneMark.hmm_PROKARYOTIC 3.25
##date: Wed May  2 12:14:54 2018
# Sequence file name: PkV-RF01_final.fasta
# Model file name: PkV-RF01_hmm_combined.mod
# RBS: true
# Model information: GeneMarkS_gcode_1
PkV-RF01	GeneMark.hmm	gene	734	1762	-1485.666949	-	0	gene_id "PkV-RF01_g1"; ID "nbis-gene-1";
PkV-RF01	GeneMark.hmm	transcript	734	1762	-1485.666949	-	0	gene_id "PkV-RF01_g1"; transcript_id "nbisL2-cds-1"; ID "nbisL2-cds-1"; Parent "nbis-gene-1"; original_biotype "mrna";
PkV-RF01	GeneMark.hmm	exon	734	1762	-1485.666949	-	.	gene_id "PkV-RF01_g1"; transcript_id "nbisL2-cds-1"; ID "nbis-exon-1"; Parent "nbisL2-cds-1";
PkV-RF01	GeneMark.hmm	CDS	734	1762	-1485.666949	-	0	gene_id "PkV-RF01_g1"; transcript_id "nbisL2-cds-1"; ID "cds-1"; Parent "nbisL2-cds-1";
PkV-RF01	GeneMark.hmm	gene	3862	4602	-1075.440697	-	0	gene_id "PkV-RF01_g2"; ID "nbis-gene-2";
PkV-RF01	GeneMark.hmm	transcript	3862	4602	-1075.440697	-	0	gene_id "PkV-RF01_g2"; transcript_id "nbisL2-cds-2"; ID "nbisL2-cds-2"; Parent "nbis-gene-2"; original_biotype "mrna";
PkV-RF01	GeneMark.hmm	exon	3862	4602	-1075.440697	-	.	gene_id "PkV-RF01_g2"; transcript_id "nbisL2-cds-2"; ID "nbis-exon-234"; Parent "nbisL2-cds-2";
PkV-RF01	GeneMark.hmm	CDS	3862	4602	-1075.440697	-	0	gene_id "PkV-RF01_g2"; transcript_id "nbisL2-cds-2"; ID "cds-2"; Parent "nbisL2-cds-2";
```
