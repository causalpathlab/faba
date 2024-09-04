
https://kb.10xgenomics.com/hc/en-us/articles/360004396971-Are-10x-Single-Cell-gene-expression-libraries-strand-specific

> Question: Are 10x Single Cell gene expression libraries
> strand-specific?
>
> Answer: Yes, 10x Single Cell solutions are strand-specific. Cell
> Ranger 'Count' counts sense-strand reads only. There is a
> possibility of legitimate antisense transcripts but they are
> unlikely to follow the same exon structure as the sense transcripts.
> We do include a measurement of the levels of antisense transcripts
> detected in Cell Ranger, which typically hovers around 1% of mapped
> reads. Those antisense transcripts won't contribute any UMI counts
> since Cell Ranger does not count UMIs from antisense reads.




### What differentiates between methylated and unmethylated m6A sites?


