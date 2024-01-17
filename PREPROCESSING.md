# Quality inspection
Fastqc
```bash
mkdir -p fastqc
fastqc Aiptasia_*/*.fastq -t 24 --nogroup -o ./fastqc &> fastqc.log &
multiqc ./fastqc/*zip -o multiqc
```

MirTrace
```bash
mirtrace qc -s meta_species_all Aiptasia_*/*.fastq -w --uncollapse-fasta --t
 20
```

Download
```bash
scp -r francesco.cicala@home.bca.unipd.it://mnt/nfs/home/francesco.cicala/Corals/Raw_Seqs/mirtrace.20240117-184928.103/OUTPUTS .

scp francesco.cicala@home.bca.unipd.it://mnt/nfs/home/francesco.cicala/Corals/Raw_Seqs/multiqc/multiqc_report.html
```