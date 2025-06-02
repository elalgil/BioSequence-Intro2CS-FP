# Bio Sequence / Elal Gilboa

## 🔧 Extensions:
- **EXTQUALITY** – Filters out low-quality reads, low-quality Kmers, and highly redundant Kmers during the mapping process.
- **EXTREVCOMP** – Determines the best orientation for each read (forward/reverse) based on the number of unique Kmers.
- **EXTCOVERAGE** – Calculates genome coverage based on unique and ambiguous read alignments.
- **EXTSIM** – Filters highly similar genomes from the reference Kmer DB and provides filtering statistics.

## 🚀 Usage
Run from the command line:
```bash
python3 main.py -t [chosen_task] [chosen_flags]
```

## ✅ Valid Flag Combinations:

### Reference Creation
```bash
-t reference -g genome.fa -k 21 -r output.kdb
```
Required:
- `-g`: genome file path (FASTA format)
- `-k`: kmer size (21–31)
- `-r`: output reference DB file path (`.kdb`)

Optional (for EXTENSIM):
- `--filter-similar`
- `--similarity-threshold THRESHOLD`

### Dump Reference
```bash
-t dumpref -r ref.kdb
```
Or:
```bash
-t dumpref -g genome.fa -k 21
```

### Align Reads
```bash
-t align -r ref.kdb -a output.aln --reads input.fq [options]
```
Or:
```bash
-t align -g genome.fa -k 21 -a output.aln --reads input.fq [options]
```
Options:
- `-m`: unique threshold
- `-p`: ambiguous threshold

### Dump Alignments
```bash
-t dumpalign -a output.aln
```
*Can only be used alone.*

Or:
```bash
-t dumpalign -r ref.kdb --reads input.fq [options]
```
Or:
```bash
-t dumpalign -g genome.fa -k 21 --reads input.fq [options]
```

## 🔌 Extension Parameters
- **EXTQUALITY**: `--min-read-quality MRQ`, `--min-kmer-quality MKQ`, `--max-genomes MG`
- **EXTREVCOMP**: `--reverse-complement`
- **EXTCOVERAGE**: `--coverage`, `--genomes g1,g2,...`, `--min-coverage MC`, `--full-coverage`

## 🧠 Design

### Data Structures:
- Reference Kmers DB is implemented as a `dict[str, List[int]]` for O(1) lookup time in pseudo-alignment.

### Classes:
- `Kmer`: holds sequence data, locations, uniqueness, and reference sources.
- `Read`: holds read data, quality, header, mapping status, and orientation.
- `Reference`: holds genome data, header, Kmers, and ID for similarity calculations.
- `Alignment`: manages mapping results, coverage, and quality statistics.
- `KmerDB`: contains all reference Kmers and optional similarity statistics.

### File Formats:
- `*.aln`: pickled `Alignment` object, compressed using `gzip`
- `*.kdb`: pickled `KmerDB` and genome lengths dictionary (from FASTA headers)

### FASTA/Q Loading:
- FASTQ: reads 4-line blocks and creates `Read` objects.
- FASTA: reads genome header + sequence to create `Reference` objects using generators.

---

Developed by Elal Gilboa  
Hebrew University – Bioinformatics Project
