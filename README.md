## Syn2Chr
Syn2Chr is a Python tools kit aiding _de novo_ assembly of pseudochromosomes from long contigs or scaffolds, especially contigs from high coverage whole genome PacBio _de novo_. Chromosome painting results and chromosome G-bands are also supported.
#### How to use
1. Search for genome synteny between the _de novo_ assembly and a reference genome with BLAST. Perfectly, the reference genome should come from same order or same familiy, and repeat masked. However, compare species with larger distance (like human vs. mouse) would still work. The BLAST output should be table format without alignment, e.g. format 6 for BLASTN. A typical command line of this step:
  - blastn -db /reference/GRCm38.repeatmask.fa -query /assembly/Mole_Falcon.fa -out blast.tbl -evalue 1e-40 -num_threads 8 -outfmt 6
2. Build a synteny map with SynBuild.py, output as .SVG figure
  - for cx in {1..22} X; do SynBuild.py ./blast.tbl -c Chr$cx -Q /assembly/Mole_Falcon.fa.fai -o ./Mole_chr$cx.svg; done
3. Based on the chromosome painting data, write a pseudochromosome synteny file (syntex of synteny file can be found in the manual).
  - ChrBuild.py ./scaffolds.fa ./synteny.txt -t -o Mole_toplevel
4. Output files
  - Mole_toplevel.fa <- Pseudochromosome fasta
  - Mole_toplevel.agp <- A Golden Path file
  - Mole_toplevel.fa <- Summary of Pseudochromosome assembly
