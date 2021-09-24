# backup





## PacBio - IsoSeq3 




# SQANTI3 - Version 3.0



## Setup

A)  Building from source , creating an environment, and activating it

 1)  > $ git clone https://github.com/ConesaLab/SQANTI3.git
 2)  > $ cd SQANTI3
 3)  > $ conda env create -f SQANTI3.conda_env.yml
 4)  > $ source activate SQANTI3.env

B) Installing dependencies

 1)  GTF to Gene Prediction is necessary for SQANTI to work; downloading it directly into /SQANTI3/utilities
 2)  wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred -P <PATH:TO>/SQANTI3/utilities/ 
     chmod +x <PATH:TO>/SQANTI3/utilities/gtfToGenePred
 4)  cDNA_Cupcake also has to be installed
 5)  > $ (SQANTI3.env)$ git clone https://github.com/Magdoll/cDNA_Cupcake.git
 6)  > $ (SQANTI3.env)$ cd cDNA_Cupcake
 7)  > $ (SQANTI3.env)$ python setup.py build
 8)  > $ (SQANTI3.env)$ python setup.py install
 9)  Before running SQANTI3 /SQANTI3/cDNA_Cupcake/sequence and /SQANTI3/cDNA_Cupcake/ have to be added to $PYTHONPATH
 10) > $ export PYTHONPATH=$PYTHONPATH:/PATH:TO/SQANTI3/cDNA_Cupcake/sequence/
 11) > $ export PYTHONPATH=$PYTHONPATH:/PATH:TO/MYO/SQANTI3/cDNA_Cupcake/



## Running it SQANTI3_qc.py
### Necessary Input

A) Assembled Transcriptome as GTF/GFF3 (with the --gtf option)
B) Reference Annotation as GTF/GFF3 
C) Reference Genome as .fasta 

> $ python sqanti3_qc.py --gtf [TranscriptomeOfInterest.Sample.gtf/gff3] [ReferenceAnnotation.gtf/gff3] [ReferemceGenome.Fasta]"


### Output of Interest for us 

1) SampleName_classification.txt
2) SampleName_sqanti_report.pdf

- The second file is a visualized pdf of the data found in the classification.txt


## Creating the Stacked Bar Plot with Fractions of Transcripts with MatlabR2021a

### Reading in the data from SQANTI's Classification.txt

- > Data=importdata('Sample_classification.txt');
- > opts=detectImportOptions('Sample_classification.txt');
- > opts.SelectedVariableNames={'structural_category'};
- > StructuralCategories=readtable('Sample_classification.txt',opts);


### Converting the data to a Numerical Array

- > c=table2cell(StructuralCategories); 
- > v=vertcat( c(:) ); 
- > Categories=categorical(v);
- > categories(Categories)
- > Cats=countcats(Categories);
- > CatsTrans=Cats(:)';

### Calculating the Fractions

- > SumOfUniqueIsoforms=sum(CatsTrans);
- > Fractions=[CatsTrans./SumOfUniqueIsoforms];

### Generating an x Vector for plotting the data correctly; replace 1 with sample numbers

- > x=[1 2 ... n]; 


### Plotting the Data
- > bar(x,Fractions,0.5,'stacked')


### Labeling the plot and making it pretty

- > set(gca,'XTickLabel',{'Flair Simple MT1','Flair Simple MT2','Flair Simple MB1','Flair Simple MB2'},                               'FontSize',16,'FontName','Times')
- > ax=gca;
- > ax.XAxis.FontSize = 20;
- > ylabel('Fraction of Transcripts','fontweight','bold','FontSize',24,'FontName','Times')
- > xlabel('SQANTI3 Quality Control','fontweight','bold','FontSize',24,'FontName','Times')
- > lgd=legend({'Full-Splice Match','Incomplete-Splice Match','Novel in Catalog','Novel not in Catalog','Genic                         Genomic','Antisense','Fusion','Intergenic'});
- > title(lgd,'Overlap Type (structural category','Fontname','Times')
- > set(lgd, 'location', 'northeastoutside','fontsize',15)
- > box off
- > ax.XTick=[1 2 3 4 5 6];
- > ax.XTickLabels={'WT1-LadderSeq', 'WT2-LadderSeq', 'WT3-LadderSeq', 'WT4-LadderSeq','m6acDNAFlair','m6aWTFlair'};
- > xtickangle(45);
- > title('Comparing Fractions of Transcripts','FontSize',26);






# GFFCompare - Visualizing structural categories

## Building from Source

1) > $ git clone https://github.com/gpertea/gffcompare
2) > $ cd gffcompare
3) > $ make release

## Executing it

### Necessary input

1) Query gff file(s) of assembled transcriptomes provided as GTF/GFF file 
2) Reference annotation provided as GTF/GFF file

### Usage

- > $ gffcompare [options]* {-i <input_gtf_list> | <input1.gtf> [<input2.gtf> .. <inputN.gtf>]}

### Extracting the structural categories from the generated .tmap file

Complete, exact match of intron chain (Classcode =)
- > $ cat sample.gtf.tmap | awk '$3=="="{print $0}' | cut -f3 | sort | wc -l 

Contained in reference; intron compatible (Classcode c)
- > $ cat sample.gtf.tmap | awk '$3=="c"{print $0}' | cut -f3 | sort | wc -l

Containment of reference; reverse containment (Classcode k)
- > $ cat sample.gtf.tmap | awk '$3=="k"{print $0}' | cut -f3 | sort | wc -l

Retained intron(s), all introns matched or retained (Classcode m)
- > $ cat sample.gtf.tmap | awk '$3=="m"{print $0}' | cut -f3 | sort | wc -l

Retained intron(s), not all introns matched/covered (Classcode n)
- > $ cat sample.gtf.tmap | awk '$3=="n"{print $0}' | cut -f3 | sort | wc -l

Multi-exon with at least one junction match (Classcode j)
- > $ cat sample.gtf.tmap | awk '$3=="j"{print $0}' | cut -f3 | sort | wc -l

Single exon transfrag partially covering an intron, possible pre-mRNA fragment(Classcode e)
- > $ cat sample.gtf.tmap | awk '$3=="e"{print $0}' | cut -f3 | sort | wc -l

Other same strand overlap with reference exons (Classcode o)
- > $ cat sample.gtf.tmap | awk '$3=="o"{print $0}' | cut -f3 | sort | wc -l

Intron match on the opposite strand (likely a mapping error) (Classcode s)
- > $ cat sample.gtf.tmap | awk '$3=="s"{print $0}' | cut -f3 | sort | wc -l

Exonic overlap on the opposite strand (like o or e but on the opposite strand) (Classcode x)
- > $ cat sample.gtf.tmap | awk '$3=="x"{print $0}' | cut -f3 | sort | wc -l

fully contained within a reference intron(Classcode i)
- > $ cat sample.gtf.tmap | awk '$3=="i"{print $0}' | cut -f3 | sort | wc -l

Contains a reference within its introns (Classcode y)
- > $ cat sample.gtf.tmap | awk '$3=="y"{print $0}' | cut -f3 | sort | wc -l

Possible polymerase run-on; no actual overlap (Classcode p)
- > $ cat sample.gtf.tmap | awk '$3=="p"{print $0}' | cut -f3 | sort | wc -l

Repeat. at least 50% bases soft-masked (Classcode r)
- > $ cat sample.gtf.tmap | awk '$3=="r"{print $0}' | cut -f3 | sort | wc -l

None of the above; unknown or intergenic (Classcode u)
- > $ cat sample.gtf.tmap | awk '$3=="u"{print $0}' | cut -f3 | sort | wc -l

### Repeat for each sample 



## Creating the Stacked Bar Plot with Fractions of Transcripts with MatlabR2021a

### Importing the data and preparing it for the plot 

First, add all the numerical outputs from the structural categories
- > SampleCategories = Classcode = + Classcode c + ... + Classcode ü

Then, create a 1x14 vektor and store the numerical outputs from the classcode command in it  
- > SampleVektor = [ Classcode = , Classcode c , ... , Classcode u]

Generate a vektor with the fractions of the structural categories
- > SampleFractionsVektor = SampleVektor ./ SampleCategories


### Repeat the three steps above for each sample


Create a matrix with the fractions for each sample 
- > FractionMatrixForPlot = [ Sample1; Sample2; ...; SampleN]


### Plotting the data and adjusting the colors 

- > bar1 = bar(FractionMatrixForPlot, 'stacked')
- > bar1(1).FaceColor=[0.282352941176471 0.352941176470588 0.4];
- > bar1(2).FaceColor=[0.647058823529412 0.929411764705882 0.576470588235294];
- > bar1(3).FaceColor=[0.541176470588235 0.709803921568627 0.419607843137255];
- > bar1(4).FaceColor=[0.27843137254902 0.490196078431373 0.305882352941176];
- > bar1(5).FaceColor=[0.882352941176471 0.92156862745098 0.356862745098039];
- > bar1(6).FaceColor=[0.929411764705882 0.658823529411765 0.250980392156863];
- > bar1(7).FaceColor=[0.968627450980392 0.407843137254902 0.125490196078431];
- > bar1(8).FaceColor=[1 0.47843137254902 0.8];
- > bar1(9).FaceColor=[0.84 0.36 0.95];
- > bar1(10).FaceColor=[0.650980392156863 0.129411764705882 1];
- > bar1(11).FaceColor=[0 0.447058823529412 0.741176470588235];
- > bar1(12).FaceColor=[0.42 0.42 0.90];
- > bar1(13).FaceColor=[0.650980392156863 0.650980392156863 0.650980392156863];
- > bar1(14).FaceColor=[0.635294117647059 0.0784313725490196 0.184313725490196];


### Labelling the plot 

- > set(gca,'XTickLabel',{'Flair Simple MT1','Flair Simple MT2','Flair Simple MB1','Flair Simple MB2'},                               'FontSize',18,'FontName','Times')
- > ax=gca;
- > ylim([0 1])
- > ax.XTick=[1 2 3 4 5 6];
- > ax.XAxis.FontSize = 20;
- > ax.XTickLabels={'LadderSeq-WT1','LadderSeq-WT2','LadderSeq-WT3','LadderSeq-WT4','m6acDNAFlair','m6aWTFlair'};
- > lgd=legend({'complete, exact match of intron chain (=)','contained in reference (c)','containment of reference (k)','retained     intron(s), all introns matched or retained (m)','retained intron(s), not all introns matched or retained (n)', ...
    'multi-exon with at least one junction match (j)','single exon transfrag covering an intron (e), possibly pre-mRNA fragment       (e)','other same strand overlap with reference exons','intron matching on opposite strand (s)', 'exonic overlap on opposite       strand (x)', 'fully contained within a reference intron (i)', ... 
    'contains a reference within its introns (y)', 'possible polymerase run-on (no overlap) (p)', 'other (u)'});
- > set(lgd, 'location', 'northeastoutside','fontsize',15)
- > title(lgd,'Overlap Type (structural category','Fontname','Times')
- > title('Gffcompare Classification of Transcripts','fontweight','bold','FontSize',26,'FontName','Times')
- > ylabel('Fraction of Transcripts','fontweight','bold','FontSize',24,'FontName','Times')
- > xlabel('Gffcompare Class Distribution','fontweight','bold','FontSize',24,'FontName','Times')
- > box off
























#### Backup ####  
  
 
  
  -Creating the Index with Hisat2
  
  > $ hisat2_extract_splice_sites.py genome.gtf > genome.ss
  
  > $ hisat2_extract_exons.py genome.gtf > genome.exon
  
  > $ hisat2-build -p 16 ReferenceGenome.fasta Output.Hisat2.Index

  > $ hisat2-build -p 16 --exon genome.exon --ss genome.ss genome.fa genome_tran

  -Aligning the sample reads against the index
  
  > $  hisat2 -p 10 -x GRCm39.Hisat2.Index --rna-strandness F -q -1 Read1_Sample.fastq -2 Read2_Sample.fastq -S Hisat2.MB1.Short

  -The reads were sorted by coordinates with SamTools
  
  > $ samtools sort {sample}.bam -o {sorted_sample}.bam

  -The sorted and aligned bam was also assembled without annotation by StringTie
  
  > $ stringtie sorted_sample.bam -o MB1.Short.StringTie.gtf -p 10
  
## Kallisto Pipeline


  -Setup
  
  > $ conda install -c bioconda kallisto

  -Indexing the sample fasta 
  
  > $ kallisto index -i {IndexFile.Output.idx} {ReferenceGenome.fa}

  -Running Kallistos' quantifying algorithm
  
  > $ kallisto quant -i {IndexFile.Output.idx} --single -l 200 -s 20 -t 10 {SampleReads.Fasta} -o /path/to/output --pseudobam

  -Extract the abundances with TPM > 0 from the SampleAbundance.tsv 

  > $ awk '$5 > 0 {print $1}' SampleAbundance.tsv > transcript_id.txt
  
  -Change the format from .txt to gtf with GtfFilter
  
  > $ gtfFilter -m whiteT -l transcript_id.txt ReferenceAnnotation.gtf CreatedSampleGTF.gtf

  -Compare assembly to reference with gffcompare
  
  > $ gffcompare -r Reference.gtf -o gffcmp Sample.gtf 



## Cufflinks Pipeline

-Creating the Index with Hisat2
  
  > $ hisat2-build -p 16 ReferenceGenome.fasta Output.Hisat2.Index

  -Aligning the sample reads against the index
  
  > $  hisat2 -p 10 -x GRCm39.Hisat2.Index --rna-strandness F -q Sample.fastq -S Hisat2.MB1.Short

  -The reads were sorted by coordinates with SamTools
  
  > $ samtools sort {sample}.bam -o {sorted_sample}.bam

  -The sorted and aligned bam was then assembled with annotation by Cufflinks

  > $ cufflinks -g {reference_annotation}.gtf {sorted_sample}.bam -p 10


