MirTrace
```bash
fastqc Aiptasia_*/*.fastq -t 24 --nogroup -o ./fastqc &> fastqc.log &


mirtrace qc -s meta_species_all Aiptasia_*/*.fastq -w --uncollapse-fasta --t
 20
```

Download
```bash
scp -r francesco.cicala@home.bca.unipd.it://mnt/nfs/home/francesco.cicala/Corals/Raw_Seqs/mirtrace.20240117-184928.103/OUTPUTS .

```