Tag-based RNA-seq reads processing pipeline, version December 29, 2017
#Modified for the SCC by Davies in 2021

#More info on the packages we will be using
#fastx_toolkit: http://hannonlab.cshl.edu/fastx_toolkit/download.html  --- trimming and quality control 
#bowtie2: http://bowtie-bio.sourceforge.net/index.shtml ---- mapping software 
#pipeline came from https://github.com/z0on/tag-based_RNAseq - check it out for additional analyses etc


#-----------------------GRABBING DATA AND SCRIPTS
# go in and grab all of the tag-seq scripts we will use along with the data files
#1. login to the server either through SCC on demand or Terminal/Mobaxterm 
cd /projectnb/bi594/YOURUSERNAME
mkdir gene_expression
cd gene_expression
ln -s /projectnb/bi594/daviessw/gene_expression/data_scripts/* .
ls -l

#let's look at everything
head -n 100 7.5A.fastq  ##these are our raw sequence files --- # Looking at first 100 lines of a sample #Grep sequencer name on fastq file -> number of sequences in the file 
head -n 100 Crep454_500.fasta #this is our standard transcriptome file
grep '>' Crep454_500.fasta | wc -l ##counts the number of contigs in your fasta file -- '>' indicates the beginning of sequence information

head Crep454_seq2iso.tab #this is a file that assigns a unique 'isogroup' for sequences depending on how it was assembled, multiple sequences can be inferred to be the same gene

head Crep454_iso2gene.tab #this is the gene annotation file
wc -l Crep454_seq2iso.tab

head Crep454_iso2go.tab #this is the Gene Ontology (GO) annotation file

#------------------------------TRIMMING FILES
# (Assuming we have many files with extension fastq, and we have fastx_toolkit installed and working)
# adaptor trimming, deduplicating, and quality filtering:

module load fastx-toolkit

# creating and launching the cleaning process for all files in the same time:
tagseq_trim_launch.pl '\.fastq$' > clean
#uses taqseq_clipper.pl function to pipe into fastxclipper -> clip off polyA tails (because mRNA), any read less than 20 bp, and based on Quality score (?), then fastxclipper again to clip illumina primers (same length and quality again), 
#then pipes into fastq_quality_filter to filter any reads that do not have 90% of read above quality score of 20 and output it into a new file called .trim for each .fastq file

##look at what is being done in clean
nano clean
tagseq_clipper.pl
#you'll see that this is where we are removing PCR duplicates that were generated during library preparation. 

nano clean
#in the next steps you'll see that these commands are clipping poly A tails, Illumina adapters, short reads (<20bp) and low quality (to keep 90% of a read must be quality score >20)

# now execute all commands written to file 'clean', preferably in array format
scc6_qsub_launcher.py -N trim -P bi594 -jobsfile clean
#this should create a trim.array.qsub and a trim_array_commands.txt files

qsub trim_array.qsub #submitting the file as a job on the terminal 

#let's see what is happening on the cluster
qstat -u daviessw (change to your username)
#To halt job 
qdel JOBID

##when the job is done, have a look in the trim.e* file
cat trim.(tab complete)
##this has all of the info for trimming. You'll see many sequences are PCR duplicates because this is TagSeq data and remember that we incorporated the degenerate bases into the cDNA synthesis. 
#goods -> reads retained ; dups -> duplicate reads removed ; noheader, N and quality were also removed 
#reads should be different lengths? 

#---------------------------------------Making the mapping database for your reference transcriptome
# download and format reference transcriptome:
#I download my reference database into the same folder, lots of people prefer having a database folder where they keep all of their databases. Up to you!
#you should start an interactive module for this for most transcriptome, but this Crepidula one is really small 

module load bowtie2
# creating bowtie2 index for your transcriptome:
bowtie2-build Crep454_500.fasta Crep454_500.fasta #prepping transcriptome so that it is mappable 
ls -l

#---------------------------------------Mapping reads to reference transcriptome

# cd where the trimmed read files are (extension "trim")
tagseq_bowtie2map.pl "trim$" Crep454_500.fasta  > maps #creating maps file 
nano maps
#utilizing bowtie2 --local function -x(mapping database) -U(unpaired reads from the trim file) -S(creating a .sam file as output) -K (number of maps to look for before it stops)

scc6_qsub_launcher.py -N maps -P bi594 -jobsfile maps
#this should create a maps.array.qsub and a maps_array_commands.txt files

qsub maps_array.qsub


###now you have individual sam files for each trimmed files

# alignment rates: #REPORT IN TABLE (raw reads, trimmed reads, mapped reads, efficiency rate, alignment rate = maps/trims) because we only tried to map trimmed reads 
nano maps.o(tab complete)

#---------------------------------------Generating read-counts-per gene 

#In newer Illumina sequencing multiple genes can map to one isogroup 

# NOTE: Must have a tab-delimited file giving correspondence between contigs in the transcriptome fasta file and genes. Typically, each gene is represented by several contigs in the transcriptome. 
head Crep454_seq2iso.tab

module load samtools
# counting hits per isogroup:
samcount_launch_bt2.pl '\.sam' Crep454_seq2iso.tab > sc
nano sc
#using samcount.pl script; count the number of sequences that map to an isogroup in the specified sam file using the bowtie alligner and creating a counts file as an output


scc6_qsub_launcher.py -N sc -P bi594 -jobsfile sc
#this should create a sc.array.qsub and a sc_array_commands.txt files

qsub sc_array.qsub

#nano sc.o(tab complete)
#you will see this: disregarding reads mapping to multiple isogroups
#we do not count reads that map to multiple places in this script, conservative approach.

#now you have individual counts files for each of your samples. Let's compile them into a single table!

####At this point we have trimmed every fastq, mapped it(.sam) and counted the reads in every .sam file per isogroup (.counts)


# assembling all counts into a single table:
#.pl (pearlscript)
expression_compiler.pl *.sam.counts > Crep_counts.txt

head Crep_counts.txt

# DONE! use your favorite R packaged (DESeq2, WGCNA) to make sense of the counts.