1/02/2020    
Giulia Manferrari 

# CLIP Metanalysis 

## CLIP 
### Experimental Design and barcodes 


## Original CLIP-seq   

In the original CLIP approach, reverse transcription needs to proceed from a universal 30 ligated adapter to a universal 50 ligated adapter, since both adapters are required for PCR ampliﬁcation. However, in over 80% of cases, the reverse transcriptase stalls at the short polypeptide left at the UV-induced crosslink site, resulting in truncated cDNAs that lack the 50 adapter, and are therefore not ampliﬁed in CLI




## iCLIP
Hupperts et al. 2014 
 
Experimental steps iCLIP 
- UV-C irradiation: covalent crosslink/bonds at sites of protein-RNA interactions (preserved in vivo RNA-RBP binding pattern) 
- Cell lysis and partial RNA digestion to obtain RNA fragments in an optimal size range
- RBP specific immunoprecipitation of the protein–RNA complexe
- RNA de-phosphorylation: this steps removes the cyclic phosphate that remains at the 3' after RNAseI digestion. If not removed the phosphate prevents adapter ligation and may favour selft-circularization
- L3 3' adapter ligation (pre-adenylated adaptor): this is an Oligo/DNA adapter pre-adenylated **(rAPP)**. It is same for every sample IDT (rAppAGATCGGAAGAGCGGTTCAG/ddC/)  

rApp (5'): pre-adenylation acts as a substrate for T4-ligase in absence of ATP 
/ddC/ (3'): to prevent the 5'-App oligo from self-ligating, the 3'-end is capped with a 3'-terminal blocking group, such as a dideoxy nucleotide (ddC) Blocking the 3'-end in this manner prevents the oligo from either circularization (by self-ligation) or concatemerization to other 5'-App oligos


- 5' RNA radioactively labellING (Y32-P-ATP)
- SSD-PAGE separation of RBP-RNA complexes 
- RNA is recovered from the nitrocellulose membrane by digesting the protein with proteinase K which leaves a polypeptide at the crosslink site
- RNA is then reverse transcribed into cDNA, which most often truncates at the polypeptide remaining at the crosslink site
- Barcoded RT primers: example 5Phos/NNAACCNNNAGATCGGAAGAGCGTCGTGgatcCTGAACCGC
- Size selection of the cDNA (Gel purification) removes free reverse transcription (RT) primers, and is followed by cDNA circularization, which attaches the second adapter to the 3' end of cDNA
- Restriction enzyme digestion linearizes the cDNA before PCR ampliﬁcation 
- PCR amplification: P5/P3 Solexa primers
     P5 Solexa primer 58 nt: AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
     P3 Solexa primer 61 nt: CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
- PCR product: product size + 58/61 nt solexa primers + barcode 
- High-troughput sequencing: after high-troughput sequencing the barcode sequencing precede the cDNA and importantly just precede the sequence at the truncation/peptide site (or cross linking site). There fore it's essential to know the nember of random nucleotide that precede and especially follow the barcode  

Sequencing Protocol: 
Standard Illumina sequencing protocol 
Reccomended: 50 nt single end run


Note on pre-adenylated adapters 
pre-adenylated (App)

irCLIP: 
in the latest protocol of irCLIP there are listed two possible L3 strategies. 
1. Standard Infrared adaptor (order from IDT as HPLC purified DNA with 5’ phosphate)  
   L3-IR-phos: /5Phos/AG ATC GGA AGA GCG GTT CAG AAA AAA AAA AAA /iAzideN/AA AAA AAA AAA A/3Bio/ 
   *The IR adaptor conjugates to IR dye 800 via the internal Azide modification, enabling the visualization of the adaptor and ligation products.*
2. Barcoded L3-???-phos 3′ adaptors (order from IDT as HPLC purified DNA with 5’ phosphate, ??? stands for the barcode, which differs between the 8 adaptors)  
  L3-ATT-phos	/5Phos/WN ATT AGA TCG GAA GAG CGG TTC AG/3Bio/


Truncated cDNA accounts for 80% of cDNA library 

RT primers : All reverse transcription primers are equippedwith 5-nt random barcodes to enable removal of PCR artifacts, plus a 4-nt barcode sequences that allow multiplexing. The 4-nt bar-codes differ by at least two nucleotides, insuring that a single pointmutation is not sufficient to convert one barcode into another.
*The quality of the library depends on the efficiency of theindividual primer (different primers might yield different library quality but they are all comparable in terms of efficiency)*

L3-App Adapter: 

### Meta analysis iCLIP, eCLIP and PAR-CLIP data
The same data analysis method applies to iCLIP and its more recent variants, including infrared CLIP (irCLIP) and enhanced CLIP (eCLIP) which also amplify truncated cDNAs.

- CLIP/PAR-CLIP: Loss of truncating cross-linking site (90% of all cases) 
- PAR-CLIP: UV-A different sensitivity (RBP-dependent)
- iCLIP: higher sensitivity as it can amplify truncated cDNA. 
  iclip and eclip have similar sensitivity - irCLIP is an order of magnitude greater 
  -eCLIP: reater care needs to be taken when analyzing data from methods that neither denature nor visualize the complexes, such as eCLIP, since it cannot be assumed that the sequenced reads represent only RNAs in contact with the protein of interest
*iclip sensitivity can be assessed at a basic level (number of unique cDNAs) and on a fucntional level(RNA maps)
crosslinking signal around regulated events should be used, whenever possible, as a more appropriate measure of data sensitivity - especially when comparing different methods


## DATA EXPLORATION/ANALYSIS
- Sensitivity: number of unique cDNAs 
- RNA-map plots 
- distribution of raw cross-link site (check for know xlink site and silenced exon for binding validation)


# SAMPLE PRE-PROCESSING

# Data processing: example with sample 1 Grot.et al data. 

The nf-core pipeling does NOT include a DE-MULTIPLEXING STEP.   
Therefore fastq files needs to be pre-processed before running the nf-core pipeline. Especially because the 5' start of the raw fastq file will contain the experimental barcode and UMI(RGB) just before the cross-linking site *(in a standard sequencing you trim everything that it's at the 5' in iclip data however it's important to preserve the 5' site as this correspond the x-linking site)

i.e this is an example of an iclip fastq 

i.e sample 1 Grot et. al. this command show the head of the fastq data. 
Here '@ERR1530360.1 SNL125:350:HWMM3ADXX:1:1108:2082:9996/1' is the header of the fastq file and does not contain any information on the UMI in the raw data
```
zcat sample1_r1.fastq.gz | head  
@ERR1530360.1 SNL125:350:HWMM3ADXX:1:1108:2082:9996/1
TGCTTTCTGATTATGGTATATTCATATAATGCAATATTATGCAGTGTAAAAAAAAGAAATTGGGGTGTGCTCTAA

```
The 5' start of the sequence will start with the randomer(or UMI or RGB) in between the Experimental barcode. 
#i.e `NNNXXXXNN = TGC*TTTC*TG` NNN: UMI/RGB
                               XXX: experimental barcode. 

The nf-core pipeline has an UMI DE-DUPLICATION step that will identify PCR duplication artifacts by identifying sequences with the same UMI. 
By definition, unique reads are reads with the same experimental barcode but different UMI. Reads with the same barcode and same UMI are considered to be PCR duplicates and needs to be de-duplicated.  
The UMI tool in the nf-core pipeline will take care of the de-duplication step if we provide the number of NN or UMI in the header. 
However it won't de-multiplex, or remove the experimental barcode.  
Therefore fastq files need to be pre-process to: 
-de-multiplex /remove sample barcode 
-add the UMI in the header 

*Note: in theory in a standard analyisis the header should contain both the UMI and the barcode. Because sometimes single points mutation can occur (?)* 

i.e This UMI tool command extract the UMI tool and place it in the header, although it leaves the experimental barcode in the reads, so that's wrong.
```
umi_tools extract --stdin=sample1_r1.fastq.gz --bc-pattern=NNNXXXXNN --log=processed.log --stdout sample1_r1.processed.fastq.gz 
#NNNXXXXNN = TGC*TTTC*TG  the umi tool removed the random barcode (TGC and TG) and inserrted this UMI in the header of the processed reads 

@ERR1530360.1_TGCTG SNL125:350:HWMM3ADXX:1:1108:2082:9996/1
TTTCATTATGGTATATTCATATAATGCAATATTATGCAGTGTAAAAAAAAGAAATTGGGGTGTGCTCTAA

```
This results in following: 

```
Read:          TAGCCGGCTTTGCCCAATTGCCAAATTTTGGGGCCCCTATGAGCTAG 
  Barcode:       NNNXXXXNN
                     |
                     v
                 TAGCCGGCT
                     |
                     V
       random-> TAG CCGG CT <- random
                     ^
                     |
                  library


 Processed read: CCGGTTGCCCAATTGCCAAATTTTGGGGCCCCTATGAGCTAG
                 ^^^^  
```
The random NN are RGB to account for PCR artifacts
The barcode (sample assignment) it's used to de-multiplex samples. 

## De-multiplexing 
To demultiplex samples use iCount. iCount is a Phyton package developed by the Ule lab. 
The package has a demultiplexing function
Each experiment is marked with a unique barcode sequence at the very beginning of the sequencing reads. Part of the barcode are also so-called randomer nucleotides that are used to identify unique cDNA molecules after mapping.

Once iCount it's installed: 
activate the icount environment where iCount commands can be run
```
conda activate icount 
```
then generate a demultiplex folder where the demultiplex files will be output
`mkdir demultiplexed`

and then extract the sample assignment and randomer sequence with the command demultiplex.
(within a wrap)
```
sbatch -c 8 --mem=32G -t 12:00:00 --wrap="iCount demultiplex sample1_r1.fastq.gz AGATCGGAAGAGCGGTTCAG NNNTTTCNN --out_dir "demultiplexed""
```

 The command expects the adapter sequence AGATCGGAAGAGCGGTTCAG, followed by the sample barcodes, in this example, expected to be present in the sequencing file. 
 #sample1 barcode TTTC 
 NNNXXXXNN = NNNTTTCNN

 *Reads that cannot be assigned to any of the specified sample barcodes (for the given number of allowed mismatches) are stored in a separate file named demux_nomatch.fastq.gz. You should have a look at such reads and try to understand why they do not conform to expectations*

Further Documentation on iCount can be found here  https://icount.readthedocs.io/en/latest/tutorial.html#preparing-iclip-data-for-mapping

Other programs can be used to extract the barcode and place the UMI in the header
- Ultraplex https://github.com/ulelab/ultraplex
- Cut&Adapt https://cutadapt.readthedocs.io/en/stable/guide.html#demultiplexing