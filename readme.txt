PARAMETERS SET:

- 1500 length filtering
- cons 100 graph (option; default)
- count node coverage 0.1 node length
- min code coverage 5 perc (option; default)
- filter mappings sum perc coverage
- filter mappings 0.9 of contig length
- incl
- no initial coverage filtering
- checkm2 (alternative - CheckM)
- GraphAligner (alternative - minigraph)

SCRIPTS:

pipe.sh - main executable
gen_bash_jobs.sh - generates a file with jobs to run pipe.sh for each species
map_contigs.sh - a script to map the sample to the graph and filter the alignments; uses
    filter_mappings.r - filters out contigs with less than n (default 0.9) of length mapped -> contig ids
    add_init.py - adds missing contigs from initial bin
    fasta_from_ids.py - gets sequences for the list of contig ids and writes to fasta file
gen_jobs.sh - in the pipeline generates a file with jobs to run map_contigs.sh script on every sample
check_job_status.py - is used in the pipeline to check if jobs which are run with submitjobs are completed
0_get_data.py - creates links to input MAGs and assembly files from global data given a species name; 
                gets GUNC, CheckM, and dRep  scores from global data table 
                (/g/scb2/bork/tschmidt/global_data/data/data.genome.v02.tsv)
add_gunc_results.py - adds a "pass GUNC" column from gunc output file to checkm2 output (after step 1.4)
select_covered_nodes.py - selects nodes covered by nmin contigs and writes them (node id per line) to file
merge_gunc_checkm2_3.py - merges GUNC and checkm2 scores from initial MAGs, GUNC refined bins, and final bins into one tsv table
plots.r - creates from the results.tsv and results3.tsv files:
            alluvial plot for MIMAG quality categories
            completeness-contamination arrow plots for initial-final and initial-gunc_refined-final bins
            completeness and contamination boxplots
+ (located at /g/scb2/bork/groudko/scripts)
refine_gunc.py - perform gunc refinement
step2_refine_gunc.py - process gunc refinement output: create species direcories and write fasta files



REQUIRED PACKAGES:
(in pipeline are used as separate conda enviroments and bork group modules)

GUNC
checkm2
BBMap
pggb (contains odgi)
GraphAligner
python: (pymongo), pandas
R: ggplot2, ggaluvial, dplyr, stringr, scales, munsell
