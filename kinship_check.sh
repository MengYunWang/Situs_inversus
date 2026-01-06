#!/usr/bin/env bash
#$ -N kinship_check
#$ -j y
#$ -cwd
#$ -q cudabig.q
#$ -S /bin/bash

module load samtools/1.11 bcftools/1.18 plink/2.00a java/jre1.8.0_201 gatk/4.3.0.0
set -euo pipefail

########## USER SETTINGS ##########
VCF37="SIT_all_HC.all.snps.tranch90.indels.tranch99.recalibrated.filtered.norm.vcf.gz"  # GRCh37
VCF38="raw_variants.snp.vcf.gz"                                                        # GRCh38

OUTDIR="kinship_crossbuild_out"
PREFIX="${OUTDIR}/combined_76subj"

# Kinship QC (be gentle here; intersection already reduces missingness)
MAF=0.05
GENO=0.02
MIND=0.1

# LD pruning
LD_WIN=200
LD_STEP=50
LD_R2=0.2
###################################

mkdir -p "${OUTDIR}/refs" "${OUTDIR}/tmp"

log(){ echo "[$(date -Is)] $*"; }
need_cmd(){ command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not in PATH"; exit 1; }; }

log "Checking tools..."
need_cmd curl; need_cmd samtools; need_cmd bcftools; need_cmd plink2; need_cmd gatk; need_cmd bgzip; need_cmd tabix

###############################################################################
# 1) References + chain
###############################################################################
log "Downloading references + chain (if missing)..."

REF38="${OUTDIR}/refs/Homo_sapiens_assembly38.fasta"
if [[ ! -s "${REF38}" ]]; then
  curl -fL \
    "https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta" \
    -o "${REF38}"
fi

REF37_GZ="${OUTDIR}/refs/human_g1k_v37.fasta.gz"
REF37="${OUTDIR}/refs/human_g1k_v37.fasta"
if [[ ! -s "${REF37}" ]]; then
  if [[ ! -s "${REF37_GZ}" ]]; then
    curl -fL \
      "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz" \
      -o "${REF37_GZ}"
  fi
  gunzip -c "${REF37_GZ}" 2>/dev/null > "${REF37}"
fi

CHAIN="${OUTDIR}/refs/hg19ToHg38.over.chain.gz"
if [[ ! -s "${CHAIN}" ]]; then
  curl -fL \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz" \
    -o "${CHAIN}"
fi

index_fasta(){
  local fa="$1"
  if [[ ! -s "${fa}.fai" ]]; then
    log "samtools faidx ${fa}"
    samtools faidx "${fa}"
  fi
  local dict="${fa%.*}.dict"
  if [[ ! -s "${dict}" ]]; then
    log "gatk CreateSequenceDictionary ${fa}"
    gatk CreateSequenceDictionary -R "${fa}" -O "${dict}"
  fi
}
index_fasta "${REF37}"
index_fasta "${REF38}"

###############################################################################
# 2) Index inputs + sample counts
###############################################################################
log "Indexing inputs (if needed)..."
# bcftools index -ft "${VCF37}" >/dev/null 2>&1 || true
# bcftools index -ft "${VCF38}" >/dev/null 2>&1 || true

log "Input sample counts:"
log "  GRCh37 = $(bcftools query -l "${VCF37}" | wc -l)"
log "  GRCh38 = $(bcftools query -l "${VCF38}" | wc -l)"

###############################################################################
# 3) Filter to PRIMARY autosomal biallelic SNPs only
###############################################################################
# Make region files (one contig per line)
log "Filtering GRCh37 to 1..22 biallelic SNPs..."

REG37="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22"
REG38="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"

V37_SNP="${OUTDIR}/tmp/v37.1-22.snps.biallelic.vcf.gz"
# bcftools view \
  # -r "${REG37}" \
  # -m2 -M2 -v snps \
  # -i 'ALT!="<NON_REF>"' \
  # -Oz -o "${V37_SNP}" "${VCF37}"
# tabix -f -p vcf "${V37_SNP}"

log "Filtering GRCh38 to chr1..chr22 biallelic SNPs..."
V38_SNP="${OUTDIR}/tmp/v38.chr1-22.snps.biallelic.vcf.gz"
# bcftools view \
  # -r "${REG38}" \
  # -m2 -M2 -v snps \
  # -i 'ALT!="<NON_REF>"' \
  # -Oz -o "${V38_SNP}" "${VCF38}"
# tabix -f -p vcf "${V38_SNP}"

###############################################################################
# 4) Rename GRCh37 contigs: 1->chr1 and strip INFO (liftover-safe)
###############################################################################
MAP_37_TO_CHR="${OUTDIR}/tmp/map37_to_chr.txt"
cat > "${MAP_37_TO_CHR}" <<'EOF'
1	chr1
2	chr2
3	chr3
4	chr4
5	chr5
6	chr6
7	chr7
8	chr8
9	chr9
10	chr10
11	chr11
12	chr12
13	chr13
14	chr14
15	chr15
16	chr16
17	chr17
18	chr18
19	chr19
20	chr20
21	chr21
22	chr22
X	chrX
Y	chrY
MT	chrM
M	chrM
EOF

log "Renaming GRCh37 1->chr1 and stripping INFO..."
V37_CHR_MININFO="${OUTDIR}/tmp/v37.chr1-22.mininfo.vcf.gz"
# bcftools annotate --rename-chrs "${MAP_37_TO_CHR}" "${V37_SNP}" -Ou \
  # | bcftools annotate -x INFO -Oz -o "${V37_CHR_MININFO}"
# tabix -f -p vcf "${V37_CHR_MININFO}"

###############################################################################
# 5) Liftover GRCh37 -> GRCh38
###############################################################################
log "LiftoverVcf GRCh37->GRCh38..."
V37_LIFT="${OUTDIR}/tmp/v37.lift38.vcf.gz"
V37_REJ="${OUTDIR}/tmp/v37.lift38.rejected.vcf.gz"

# gatk --java-options "-Xmx16g" LiftoverVcf \
  # -I "${V37_CHR_MININFO}" \
  # -O "${V37_LIFT}" \
  # -CHAIN "${CHAIN}" \
  # -REJECT "${V37_REJ}" \
  # -R "${REF38}" \
  # --RECOVER_SWAPPED_REF_ALT false \
  # --WARN_ON_MISSING_CONTIG true

# tabix -f -p vcf "${V37_LIFT}" || true
# tabix -f -p vcf "${V37_REJ}" 2>/dev/null || true

###############################################################################
# 6) Normalize both to same GRCh38 reference
###############################################################################
log "Normalizing lifted GRCh37->GRCh38..."
V37_NORM="${OUTDIR}/tmp/v37.lift38.norm.vcf.gz"
# bcftools norm -f "${REF38}" -m-any -Oz -o "${V37_NORM}" "${V37_LIFT}"
# tabix -f -p vcf "${V37_NORM}"

log "Normalizing native GRCh38..."
V38_NORM="${OUTDIR}/tmp/v38.norm.vcf.gz"
# bcftools norm -f "${REF38}" -m-any -Oz -o "${V38_NORM}" "${V38_SNP}"
# tabix -f -p vcf "${V38_NORM}"

###############################################################################
# 7) INTERSECT SITES (key step to avoid extreme missingness)
###############################################################################
log "Creating intersection of sites shared by BOTH cohorts..."

# Create site-only VCFs
V37_SITES="${OUTDIR}/tmp/v37.sites.vcf.gz"
V38_SITES="${OUTDIR}/tmp/v38.sites.vcf.gz"
# bcftools view -G -Oz -o "${V37_SITES}" "${V37_NORM}"
# bcftools view -G -Oz -o "${V38_SITES}" "${V38_NORM}"
# tabix -f -p vcf "${V37_SITES}"
# tabix -f -p vcf "${V38_SITES}"

# Compute intersection positions (bcftools isec outputs 0000/0001/0002/0003)

log "bcftools isec (shared variants)..."

SHARED_SITES="${OUTDIR}/tmp/shared_sites.vcf.gz"

# bcftools isec \
  # -n=2 \
  # -w1 \
  # -Oz \
  # -o "${SHARED_SITES}" \
  # "${V37_SITES}" "${V38_SITES}"

# tabix -f -p vcf "${SHARED_SITES}"


log "Shared site count: $(bcftools view -H "${SHARED_SITES}" | wc -l)"

# Keep only those shared sites in each cohort (preserve genotypes)
V37_SHARED="${OUTDIR}/tmp/v37.shared.vcf.gz"
V38_SHARED="${OUTDIR}/tmp/v38.shared.vcf.gz"
# bcftools view -T "${SHARED_SITES}" -Oz -o "${V37_SHARED}" "${V37_NORM}"
# bcftools view -T "${SHARED_SITES}" -Oz -o "${V38_SHARED}" "${V38_NORM}"
# tabix -f -p vcf "${V37_SHARED}"
# tabix -f -p vcf "${V38_SHARED}"

###############################################################################
# 8) Merge cohorts on shared sites
###############################################################################
log "Merging cohorts on shared sites..."
MERGED_CHR="${OUTDIR}/tmp/all76.GRCh38.shared.chr.vcf.gz"
# bcftools merge -m none -Oz -o "${MERGED_CHR}" "${V38_SHARED}" "${V37_SHARED}"
# tabix -f -p vcf "${MERGED_CHR}"

log "Merged sample count: $(bcftools query -l "${MERGED_CHR}" | wc -l)"

###############################################################################
# 9) Rename chr* -> numeric, keep ONLY 1..22 for PLINK
###############################################################################
MAP_CHR_TO_NUM="${OUTDIR}/tmp/mapchr_to_num.txt"
cat > "${MAP_CHR_TO_NUM}" <<'EOF'
chr1	1
chr2	2
chr3	3
chr4	4
chr5	5
chr6	6
chr7	7
chr8	8
chr9	9
chr10	10
chr11	11
chr12	12
chr13	13
chr14	14
chr15	15
chr16	16
chr17	17
chr18	18
chr19	19
chr20	20
chr21	21
chr22	22
chrX	X
chrY	Y
chrM	MT
EOF

log "Renaming chr* -> numeric..."
MERGED_NUM="${OUTDIR}/tmp/all76.shared.numeric.vcf.gz"
# bcftools annotate --rename-chrs "${MAP_CHR_TO_NUM}" -Oz -o "${MERGED_NUM}" "${MERGED_CHR}"
# tabix -f -p vcf "${MERGED_NUM}"

log "Keeping only 1..22..."
MERGED_NUM_MAIN="${OUTDIR}/tmp/all76.shared.numeric.main22.vcf.gz"
# bcftools view -r "${REG37}" -Oz -o "${MERGED_NUM_MAIN}" "${MERGED_NUM}"
# tabix -f -p vcf "${MERGED_NUM_MAIN}"

log "Final VCF sample count: $(bcftools query -l "${MERGED_NUM_MAIN}" | wc -l)"
log "Final VCF variant count: $(bcftools view -H "${MERGED_NUM_MAIN}" | wc -l)"

###############################################################################
# 10) PLINK2 -> LD prune -> KING
###############################################################################
log "PLINK2 conversion + QC filters..."
plink2 \
  --vcf "${MERGED_NUM_MAIN}" \
  --vcf-half-call missing \
  --chr 1-22 \
  --snps-only just-acgt \
  --max-alleles 2 \
  --maf "${MAF}" \
  --geno "${GENO}" \
  --mind "${MIND}" \
  --make-bed \
  --out "${PREFIX}"

log "LD pruning..."
plink2 \
  --bfile "${PREFIX}" \
  --indep-pairwise "${LD_WIN}" "${LD_STEP}" "${LD_R2}" \
  --out "${PREFIX}.prune"

log "Making pruned dataset..."
plink2 \
  --bfile "${PREFIX}" \
  --extract "${PREFIX}.prune.prune.in" \
  --make-bed \
  --out "${PREFIX}.pruned"

log "KING kinship..."
plink2 \
  --bfile "${PREFIX}.pruned" \
  --make-king-table \
  --out "${PREFIX}.king"

log "DONE. Key output: ${PREFIX}.king.kin0"
log "Liftover rejects: ${V37_REJ}"

log "All done."
