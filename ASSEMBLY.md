# Shorstacks

## Preprocessing libs

Test wich adapter w/ trim_galone


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

conda env list
conda activate ShortStack4 
ShortStack --version



#readfile=`ls -x /mnt/nfs/home/francesco.cicala/Corals/Raw_Seqs/Aiptasia_*/*.fastq`
readfile=`ls -x /mnt/nfs/home/francesco.cicala/Corals/Raw_Seqs/Aiptasia_Argonaute-CLIP/*.fastq`
ShortStack --genomefile GENOMES.fa --known_miRNAs ALL-mat.fa --dn_mirna --outdir ShortStack_"$(date +%Y%m%d)"_out --threads 24 --dicermax 30 --mmap u --mincov 0.8 --pad 1 --readfile $readfile --autotrim &>> "ShortStack_"$(date +%Y%m%d)".log" &
