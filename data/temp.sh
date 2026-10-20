
## check 10x BAM


ls -l | wc -l

# samtools view pbmc3k_possorted_genome_bam_filterd.bam     View and convert BAM to SAM
# Extract the BAM header and write to header_filted_sam
samtools view -H pbmc3k_possorted_genome_bam.bam | head


samtools view pbmc3k_possorted_genome_bam.bam | cut -f 12- | tr "\t" "\n" | grep "CB:Z:" | uniq | head

# samtools view pbmc3k_possorted_genome_bam_filterd.bam     View and convert BAM to SAM
# grep 'xf:i:25'                                                Find counted UMIs
# cut -f 12-                                                    Remove first 12 columns SAM file
# tr "\t" "\n"                                                  Replace all tabs with a new line
# grep 'CB:Z:'                                                  Find lines that start with CB:Z



samtools view pbmc3k_possorted_genome_bam.bam | grep "xf:i:25" | cut -f 12- | grep "CB:Z:CGAGAAGTCAAGTTGC" | tr "\t" "\n" | grep "GX:Z" |awk "{count[$1]++} END {for(j in count) print count[j], j}" | sort -nr > feature_counts_per_barcode.txt

# samtools view pbmc_1k_v3_possorted_genome_bam_filterd.bam     View and convert BAM to SAM
# grep "xf:i:25"                                                Grep counted reads
# cut -f 12-                                                    Remove first 12 columns SAM file
# grep "CB:Z:CGAGAAGTCAAGTTGC"                                  Print lines with "CB:Z:CGAGAAGTCAAGTTGC". ***Replace this with your barcode of interest***
# tr "\t" "\n"                                                  Replace all tabs with a new line
# grep "CB:Z:""                                                 Grep lines that start with CB:Z
# awk "{count[$1]++} END {for(j in count) print count[j], j}"   Count each occurrence of a barcode
# sort -nr                                                      Sort numerically and in reverse order





# Filter BAM file for records with 'xf:i:25' and write that to
# body_filtered_sam. body_filtered_sam is on the sam format and will
# not include the BAM header.
samtools view pbmc3k_possorted_genome_bam.bam | LC_ALL=C grep "xf:i:25" > body_filtered_sam

# Extract the BAM header and write to header_filted_sam
samtools view -H pbmc3k_possorted_genome_bam.bam > header_filted_sam

# Combine the header and extracted records
cat header_filted_sam body_filtered_sam > pbmc_1k_v3_possorted_genome_bam_filterd.sam

# Convert SAM to BAM
samtools view -b pbmc_1k_v3_possorted_genome_bam_filterd.sam > pbmc_1k_v3_possorted_genome_bam_filterd.bam
