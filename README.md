# WGlab_small_RNA_analysis
This is the PERL script and the associated annotation for C. elegans small RNA analysis
The genome annotation used is WS215
The alignment software is Bowtie 0.12.7
A genome annotation index is processed by bowtie and then saved in the bowtie index folder
The working folder contains the PERL script and associated annotation, a folder 'bowtie0127' for temporary analysis file, and a 'fastaq' folder containing the sample sequences.
The sample sequence is in the fasta format with ">ID_xread-number\nsequence\n". For example, ">WG20222222_x323\nAGGGGGGCCAATT\n"
