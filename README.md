# GoBi-2

Implement a simple paired-end read simulator considering the following parameters:

• -length <int> : read length  
• -frlength <int> -SD <int>: fragment length distribution. The input is the mean and standard  
deviation of a normal distribution. This must be used to draw fragment lengths (you
have to round the drawn numbers). For simplicity re-draw the fragment length if the value is
smaller than length or larger than the transcript length.  
• -readcounts: table of gene_id, transcript_id, count tuples: For every transcript simulate
count many reads. For each read pair, select a fragment length (see frlength), and draw a
start position from a uniform distribution of all possible start positions on the transcript.  
• -mutationrate: mutation rate in percent (between 0.0 and 100.0): percent of the simulated  
bases to be changed from the original.  
• -fasta: genome FASTA file  
• -fidx: genome FASTA file index  
• -gtf: annotation file for the transcript locations  
• -od: output directory (where the files fw.fastq, rw.fastq, read.mappinginfo will be written)  
  
The simulator should write three files, two for the simulated paired-end sequences in **FASTQ format**
(fw.fastq, rw.fastq, one FASTQ file for the first fragment, one for the second, set the quality score
to the maximum for all bases), and a tab separated file (read.mappinginfo) with the following
headers:

• readid: integer (starting from 0)  
• chr id: chromosome  
• gene id  
• transcript id  
• fw regvec: genomic region vector for the forward read (1-based end exclusive)  
• rw regvec: genomic region vector for the reverse read (1-based end exclusive)  
• t fw regvec: region vector for the forward read in transcript coordinates (0-based end exclusive)  
• t rw regvec: region vector for the reverse read in transcript coordinates (0-based end exclusive)  
• fw mut: mutated positions in the forward read (comma separated integer list)  
• rw mut: mutated positions in the reverse read (comma separated integer list)  

The format for genomic region vectors is:  
start1-end1(|startx-endx)+

Use your implementation and simulate reads with the following settings:  
**length=75, frlength=200, SD=80, mutationrate=1.0%** on the input file  
**readcounts=readcounts.simulation** containing the number of reads to simulate per transcript.  
