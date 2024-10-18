# Genome metadata
genome="danRer11"
version="GRCz11.104"
 
 
# Set up source and build directories
build="${genome}-${version}-build"
mkdir -p "$build"
 
 
# Download source files if they do not exist in reference-sources/ folder
 
source="${genome}-${version}-reference-sources"
mkdir -p "$source"
 
 
fasta_url="http://ftp.ensembl.org/pub/release-104/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz"
fasta_in="${source}/Danio_rerio.GRCz11.dna.primary_assembly.fa"
gtf_url="http://ftp.ensembl.org/pub/release-104/gtf/danio_rerio/Danio_rerio.GRCz11.104.gtf.gz"
gtf_in="${source}/Danio_rerio.GRCz11.104.annotation.gtf"
motifs_url="http://jaspar.genereg.net/download/CORE/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
motifs_in="${source}/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
 
 
if [ ! -f "$fasta_in" ]; then
    curl -sS "$fasta_url" | zcat > "$fasta_in"
fi
if [ ! -f "$gtf_in" ]; then
    curl -sS "$gtf_url" | zcat > "$gtf_in"
fi
if [ ! -f "$motifs_in" ]; then
    curl -sS "$motifs_url" > "$motifs_in"
fi
 
 
# Modify sequence headers in the Ensembl FASTA to match the file
# "GRCh38.primary_assembly.genome.fa" from GENCODE. Unplaced and unlocalized
# sequences such as "KI270728.1" have the same names in both versions.
#
# Input FASTA:
#   >1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF
#
# Output FASTA:
#   >chr1 1
fasta_modified="$build/$(basename "$fasta_in").modified"
# sed commands:
# 1. Replace metadata after space with original contig name, as in GENCODE
# 2. Add "chr" to names of autosomes and sex chromosomes
# 3. Handle the mitochrondrial chromosome #Note: the chr prefix is not added to unplaced and unlocalized
cat "$fasta_in" | sed -E 's/^>(\S+).*/>\1 \1/' | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' | sed -E 's/^>MT />chrM /' > "$fasta_modified"
 
 
############### Do not need this step for zebrafish Ensembl annotation #####################
# Remove version suffix from transcript, gene, and exon IDs in order to match
# previous Cell Ranger reference packages
#
# Input GTF:
#     ... gene_id "ENSG00000223972.5"; ...
# Output GTF:
#     ... gene_id "ENSG00000223972"; gene_version "5"; ...
## gtf_modified="$build/$(basename "$gtf_in").modified"
# Pattern matches Ensembl gene, transcript, and exon IDs for human or mouse:
## ID="(ENS(DAR)?[GTEP][0-9]+)\.([0-9]+)"
##cat "$gtf_in" \
##    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
##    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
##    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
##    | sed -E 's/protein_id "'"$ID"'";/protein_id "\1"; protein_version "\3";/' \
 
##    > "$gtf_modified"
############################################################################################
############## Instead, add Chr prefix to zebrafish ensemble annotation ####################
## INPUT
#gtf_in="${source}/Danio_rerio.GRCz11.104.annotation.gtf"
 
## OUTPUT
gtf_modified="$build/$(basename "$gtf_in").modified"
 
#Note: the chr prefix is not added to unplaced and unlocalized
cat $gtf_in |
    awk -F"\t" -vOFS="\t" '$1~/^[1-9]/ {$1="chr"$1}
                           $1=="MT" {$1="chrM"}
                           {print $0}' > $gtf_modified
 
 
# Define string patterns for GTF tags
# NOTES:
# - Since GENCODE release 31/M22 (Ensembl 97), the "lincRNA" and "antisense"
#   biotypes are part of a more generic "lncRNA" biotype.
# - These filters are relevant only to GTF files from GENCODE. The GTFs from
#   Ensembl release 98 have the following differences:
#   - The names "gene_biotype" and "transcript_biotype" are used instead of
#     "gene_type" and "transcript_type".
#   - Readthrough transcripts are present but are not marked with the
#     "readthrough_transcript" tag.
#   - Only the X chromosome versions of genes in the pseudoautosomal regions
#     are present, so there is no "PAR" tag.
BIOTYPE_PATTERN=\
"(protein_coding|lincRNA|antisense|\
IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
TR_V_pseudogene|TR_J_pseudogene)"
GENE_PATTERN="gene_biotype \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_biotype \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""
PAR_PATTERN="tag \"PAR\""
 
 
# Construct the gene ID allowlist. We filter the list of all transcripts
# based on these criteria:
#   - allowable gene_type (biotype)
#   - allowable transcript_type (biotype)
#   - no "PAR" tag (only present for Y chromosome PAR)
#   - no "readthrough_transcript" tag
# We then collect the list of gene IDs that have at least one associated
# transcript passing the filters.
cat "$gtf_modified" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | grep -Ev "$PAR_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "${build}/gene_allowlist"
 
 
# Filter the GTF file based on the gene allowlist
gtf_filtered="${build}/$(basename "$gtf_in").filtered"
# Copy header lines beginning with "#"
grep -E "^#" "$gtf_modified" > "$gtf_filtered"
# Filter to the gene allowlist
grep -Ff "${build}/gene_allowlist" "$gtf_modified" \
    >> "$gtf_filtered"
 
 
 
# Change motif headers so the human-readable motif name precedes the motif
# identifier. So ">MA0004.1    Arnt" -> ">Arnt_MA0004.1".
motifs_modified="$build/$(basename "$motifs_in").modified"
awk '{
    if ( substr($1, 1, 1) == ">" ) {
        print ">" $2 "_" substr($1,2)
    } else {
        print
    }
}' "$motifs_in" > "$motifs_modified"

############# Add EGFP information to FASTA and GTF file ###################

# Create EGFP fasta file
EGFP_fa="${build}/genome.EGFP.fa"
echo ">EGFP
ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAA
GTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGC
TGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAG
CAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTA
CAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGG
ACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAAC
GGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACAC
CCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACG
AGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA" > "$EGFP_fa"

# Check sequence length:
cat "$EGFP_fa" | grep -v "^>" | tr -d "\n" | wc -c
720

# Create EGFP gtf file
EGFP_gtf="${build}/genes.EGFP.gtf"
echo "EGFP    unknown gene    1       720     .       +       .       gene_id "EGFP"; gene_name "EGFP"; gene_biotype "protein_coding";
EGFP    unknown transcript      1       720     .       +       .       gene_id "EGFP"; transcript_id "EGFP"; gene_name "EGFP"; gene_biotype "protein_coding";
EGFP    unknown exon    1       720     .       +       .       gene_id "EGFP"; transcript_id "EGFP"; gene_name "EGFP"; gene_biotype "protein_coding";" > "$EGFP_gtf"

# Add EGFP to fasta and gtf:
gtf_filtered2="${build}/$(basename "$gtf_filtered")2"
fasta_modified2="${build}/$(basename "$fasta_modified")2"

cat "$gtf_filtered" "$EGFP_gtf" > "$gtf_filtered2" 
cat "$fasta_modified" "$EGFP_fa" > "$fasta_modified2" 


#######################################
 
# Create a config file
config_in="${build}/config"
echo """{
    organism: \"Danio_rerio\"
    genome: [\""$genome"\"]
    input_fasta: [\""$fasta_modified2"\"]
    input_gtf: [\""$gtf_filtered2"\"]
    input_motifs: \""$motifs_modified"\"
    non_nuclear_contigs: [\"chrM\"]
}""" > "$config_in"
 

# Create reference package
export PATH=/bar/ichen/cellranger-arc-2.0.0:$PATH
cellranger-arc mkref --config="$config_in" --nthreads=16 --memgb=80 > cellranger-arc-2.0.0_mkref.log 2>&1 &


# Make gtf only with standard chromosomes and EGFP 
cat Danio_rerio.GRCz11.104.annotation.gtf.filtered2 | grep "#" > head.gtf
cat Danio_rerio.GRCz11.104.annotation.gtf.filtered2 | grep -v "#" |awk -F"\t" -vOFS="\t" '{if(length($1) < 6) print $0}' > Danio_rerio.GRCz11.104.annotation.onlyChr_EGFP.noHead.gtf
cat head.gtf Danio_rerio.GRCz11.104.annotation.onlyChr_EGFP.noHead.gtf > Danio_rerio.GRCz11.104.annotation.onlyChr_EGFP.gtf
cp config ./config_onlyChr # then change gtf file inside config_onlyChr

# Use only Chr and EGFP fasta

ichen@paleAle:/scratch/ichen/cellranger-arc-2.0.0/danRer11-GRCz11.104-build$ less config_onlyChr
{
    organism: "Danio_rerio"
    genome: ["danRer11"]
    input_fasta: ["danRer11-GRCz11.104-build/danRer11_onlyChr_wEGFP.fa"]
    input_gtf: ["danRer11-GRCz11.104-build/Danio_rerio.GRCz11.104.annotation.onlyChr_EGFP.gtf"]
    input_motifs: "danRer11-GRCz11.104-build/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt.modified"
    non_nuclear_contigs: ["chrM"]
}
export PATH=/bar/ichen/cellranger-arc-2.0.0:$PATH
cellranger-arc mkref --config=/scratch/ichen/cellranger-arc-2.0.0/danRer11-GRCz11.104-build/config_onlyChr --nthreads=24 --memgb=90 > cellranger-arc-2.0.0_mkref_onlyChr_EGFP.log 2>&1 &
