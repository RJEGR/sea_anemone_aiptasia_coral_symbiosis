# Assembly microRNAs

## Preprocessing libs

Test wich adapter w/ trim_galone

```bash
for i in $(ls -x /mnt/nfs/home/francesco.cicala/Corals/Raw_Seqs/Aiptasia_*/*
.fastq); do trim_galore $i --small_rna --cores 12; done

mkdir -p REPORTS
mv *.txt REPORTS
```

## Mirtrace

```bash
# Using well directory 

#cd /qc_passed_reads.all.uncollapsed

#  grep "^>" SRR2716058_trimmed.fasta  | awk '{print $2}' |  sort | uniq -c
# for i in $(ls *.fasta); do echo $i;  grep "^>" $i  | awk '{print $2}' |  sort | uniq -c; done

# for i in $(ls *.fasta); do seqkit fx2tab $i -l -g -H > ${i%.fasta}.profiling; done &

# upgrade version of profiling
 # seqkit fx2tab SRR2716058_trimmed.fasta -l -g -H -j 6 | awk '{$3 = substr($3, 1,2)} 1' | head

for i in $(ls *.fasta); do seqkit fx2tab $i -l -g -H -j 6 | awk '{$3 = substr($3, 1,2)} 1' > ${i%.fasta}.profiling; done &


for i  in $(ls *.fasta); do cat $i | seqkit grep -n -r -p "rnatype:mirna" -p "rnatype:unknown" >  ${i%.fasta}.mirna.unknown.fa; done

mkdir -p PROFILING_BY_READ_LENGTH

for i in $(ls *.profiling); do tar -czvf ${i}.tar.gz $i; done

mv *.tar.gz PROFILING_BY_READ_LENGTH
scp -r francesco.cicala@home.bca.unipd.it://mnt/nfs/home/francesco.cicala/Corals/Clean_Seqs/mirtrace.20240207-164934.646/qc_passed_reads.all.uncollapsed/PROFILING_BY_READ_LENGTH

# Process it in r

# or,
for i in $(ls *.profiling); do cat $i | awk '{print $2,$4}' | sort -k 2 | uniq -c > ${i%.profiling}.reads.length.summarized; done

```








## 1) Reference
# 1.1) decompress

```bash
zcat GCA_001417965.1_A
iptasia_genome_1.1_genomic.fna.gz > GCA_001417965.1_Aiptasia_genome_1.1_genomic.fna
```

# Concat
```bash
grep -c "^>" *.fna 

# GCA_001417965.1_Aiptasia_genome_1.1_genomic.fna:4312
# GCA_001939145.1_ASM193914v1_genomic.fna:9688
# GENOMES.fna:14000

#/mnt/nfs/home/francesco.cicala/Corals/Reference_Genome/Aipstasia_genome/GCA_001417965.1_Aiptasia_genome_1.1_genomic.fna.gz
#/mnt/nfs/home/francesco.cicala/Corals/Reference_Genome/Symbiodinium_microadriaticum/data/GCA_001939145.1

cat GCA_001417965.1_Aiptasia_genome_1.
1_genomic.fna GCA_001939145.1_ASM193914v1_genomic.fna > GENOMES.fna

```

## X) Run
```bash
conda env list
conda activate ShortStack4 
ShortStack --version


#readfile=`ls -x /mnt/nfs/home/francesco.cicala/Corals/Raw_Seqs/Aiptasia_*/*.fastq`

#readfile=`ls -x /mnt/nfs/home/francesco.cicala/Corals/Clean_Seqs/*.fq`

readfile=`ls -x /mnt/nfs/home/francesco.cicala/Corals/Clean_Seqs/mirtrace.20240207-164934.646/qc_passed_reads.rnatype_unknown.uncollapsed/*.mirna.unknown.fasta`

ShortStack --genomefile GENOMES.fa --known_miRNAs ALL-mat.fa --dn_mirna --outdir ShortStack_"$(date +%Y%m%d)"_out --threads 24 --dicermax 30 --mmap u --mincov 0.8 --readfile $readfile &>> "ShortStack_"$(date +%Y%m%d)".log" &

mkdir -p OUTPUTS

cp -r Counts.txt mir.fasta Results.txt strucVis/ alignment_details.tsv OUTPUTS/

 scp -r francesco.cicala@home.bca.unipd.it://mnt/nfs/home/francesco.cicala/Corals/RICARDO/ShortStack_20240203_out/OUTPUTS .

```



### MIRDEE2

```bash
bowtie-build GENOMES.fa GENOMES

for i in $(ls /mnt/nfs/home/francesco.cicala/Corals/Clean_Seqs/mirtrace.20240207-164934.646/qc_passed_reads.rnatype_unknown.uncollapsed/*.mirna.unknown.fasta); do remove_white_space_in_id.pl $i > ${i%.fasta}.whitespac.fa

ls * | grep whitespac.fa | sort > f2.tmp && cut -d '.' -f 1 f2.tmp > f1.tmp

paste f1.tmp f2.tmp | column -t > config.txt && rm *tmp


mapper.pl config.txt -d -c -j -l 18 -m -s reads_collapsed.fa -p GENOMES
```

continue w https://github.com/RJEGR/Small-RNASeq-data-analysis/blob/master/RAW_TUTORIAL_BKP/TUTORIAL.md#3-optional-mirdeep2pl-to-detect-novel-mirna-candidates