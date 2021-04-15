# iCLIP analysis notes 
example on Grot et al 2017 data 
based on the following papers:
A.Bush et al 2020 


## General Quality check FASTQC
The starting point of the pipeline is a file of sequencing reads infastqfile format (or in a compacted version asfastq.gz), which includes the reads of all samples that were multiplexed in the sequencing run. The quality of the sequencing run can be checked usingFastQC:
```
fastqc ––extract ––nogroup ––outdir <outdir><data.fastq.gz>
```
fastqc ––extract ––nogroup ––outdir pre-processing_fastqc data.fastq.gz ERR1530360.fastq.gz


Among the quality measures reported by FastQC, we mainly focus on 

- Per Base Sequence Quality:

Especially for longer reads, the Per Base Sequence Quality (Phred score) can drop towards the end of the reads. In case of extremely low qualities, we recommend trimming the 3 ends of reads.
Post demultiplexing and UMI extraction processing the 5' Phred score should improve, while at 3' you might expect a drop in quality (possibly due to 3'adapter removal)


- Per Base Sequence Content: 
The Per Base Sequence Content over all iCLIP reads shows the UMI (here positions 1-3 and 8-9) and the experimental barcode (here positions 4-7; Fig. 3). The pattern from position 10 onward reflects the RNA binding preference of the RBP, and consequently varies between the studied proteins. Keep in mind that because this picture looks different from standard RNA-sequencing data, FastQC might red-flag some of its checks, but these failed check are not necessarily meaningful in the context of iCLIP data.


## Demultiplexing 
Following initial quality control and quality filtering, demultiplexingand adapter trimming are performed on the quality-filtered data usingFlexbar.In addition to separating the reads of the different samples,this step extracts the barcode region from the 5end of the reads(––barcode-trim-end LTAIL), and trims adapters from the3end ifpresent(––adapter-trim-end RIGHT).We recommend not allowing any mismatches for barcode matching (set by––barcode-error-rate 0),while an error rate of 0.1 is acceptable for trimming adapter%T%C%Aposition in read (nt)fraction per base (%) sequence content across all bases Fig.3.Quality control of the iCLIP sequencing reads.  (––adapter-seq adapter.seq ––adapter-error-rate 0.1).
 Inorder to remove all adapter traces, even if only the very first nucleo-tides of them were sequenced, we commonly require just 1 nt ofoverlap between the read)end and the beginning of the adapter(––adapter-min-overlap 1). Since very short sequences are likely to overlap by chance,this stringent setting means that many reads will be trimmed by a few extra bases, but ensures that even the shortest remain in g adapter fragments are trimmed off.Only trimmed read swith are maining length of at least min Read Length are kept for further analysis. Demultiplexing,adapter trimming and barcode removal can either be done separately or in one step. For the latter, we commonly use Flexbar with the following parameters: 
 flexbar -r <data.filtered.fastq.gz>––zip-output GZ––barcodes barcodes.fasta––barcode-unassigned––barcode-trim-end LTAIL––barcode-error-rate 0––adapter-seq adapter.seq––adapter-trim-end RIGHT––adapter-error-rate 0.1––adapter-min-overlap 1––min-read-length minReadLength––umi-tags Busch, et al.Methods 178 (2020)