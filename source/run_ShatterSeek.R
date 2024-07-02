library(ShatterSeek)

# SAMPLE chr start end event copy.number bins median
# copy.number is what is needed
CNA <- read.table("data/Copy_Number_Estimates/MMRF_CoMMpass_IA21_genome_ichorCNA.seg", header = T)

# SAMPLE SVTYPE CHROM POS ID REF CHR2 ENDPOSSV CIPOS MAPQ PE SR SRQ CE CT IMPRECISE PRECISE HOMLEN RDRATIO LENSV RGENUPS_ENSG RGENUPS HUGO RGENUPS_Distance	RGENUPS_Strand	RGENISEC_ENSG	RGENISEC_HUGO	RGENISEC_Distance	RGENISEC_Strand	RGENDNS_ENSG	RGENDNS_HUGO	RGENDNS_Distance	RGENDNS_Strand	LGENUPS_ENSG	LGENUPS_HUGO	LGENUPS_Distance	LGENUPS_Strand	LGENISEC_ENSG	LGENISEC_HUGO	LGENISEC_Distance	LGENISEC_Strand	LGENDNS_ENSG	LGENDNS_HUGO	LGENDNS_Distance	LGENDNS_Strand	NORMAL_GT	TUMOR_GT	NORMAL_GL_RR	TUMOR_GL_RR	NORMAL_GL_RA	TUMOR_GL_RA	NORMAL_GL_AA	TUMOR_GL_AA	NORMAL_GQ	TUMOR_GQ	NORMAL_FT	TUMOR_FT	NORMAL_RC	TUMOR_RC	NORMAL_RCL	TUMOR_RCL	NORMAL_RCR	TUMOR_RCR	NORMAL_CN	TUMOR_CN	NORMAL_DR	TUMOR_DR	NORMAL_DV	TUMOR_DV	NORMAL_RR	TUMOR_RR	NORMAL_RV	TUMOR_RV	NORMAL_RCALT	TUMOR_RCALT	NORMAL_RDISTDISC1	TUMOR_RDISTDISC1	NORMAL_RDISTDISC2	TUMOR_RDISTDISC2	NORMAL_RCDIS1	TUMOR_RCDIS1	NORMAL_RCDIS2	TUMOR_RCDIS2	CONSENSUS	isBalTra	dupBalTraToRemove
SV_d <- read.table("data/Structural_Event_Files/MMRF_CoMMpass_IA21_genome_delly.tsv", header = T)

# SAMPLE	CHROM	POS	ID	REF	ALT	QUAL	FILTER	IMPRECISE	SVTYPE	SVLEN	END	CIPOS	CIEND	CIGAR	MATEID	EVENT	HOMLEN	HOMSEQ	SVINSLEN	SVINSSEQ	LEFT_SVINSSEQ	RIGHT_SVINSSEQ	BND_DEPTH	MATE_BND_DEPTH	SOMATIC	SOMATICSCORE	JUNCTION_SOMATICSCORE	TRA	INV	INVDUP	DUPTANDEM	CHR2	ENDPOSSV	RO	UNTESTED	RGENUPS	RGENISEC	RGENDNS	LGENUPS	LGENISEC	LGENDNS	GENELISTbtwBP	LENSV	NORMAL_PR	NORMAL_SR	NORMAL_RUT1	NORMAL_RUT2	NORMAL_RDISTDISC1	NORMAL_RDISTDISC2	NORMAL_RCDIS1	NORMAL_RCDIS2	NORMAL_DRNOISE1	NORMAL_DRNOISE2	NORMAL_EVALRCNOISE	NORMAL_CLMQ1	NORMAL_CLMQ2	TUMOR_PR	TUMOR_SR	TUMOR_RUT1	TUMOR_RUT2	TUMOR_RDISTDISC1	TUMOR_RDISTDISC2	TUMOR_RCDIS1	TUMOR_RCDIS2	TUMOR_DRNOISE1	TUMOR_DRNOISE2	TUMOR_EVALRCNOISE	TUMOR_CLMQ1	TUMOR_CLMQ2
SV_m <- read.table("data/Structural_Event_Files/MMRF_CoMMpass_IA21_genome_manta.tsv", header = T)

SV_m <- SV_m |> dplyr::mutate(strand1 = substr(RO, 1, 1), strand2 = substr(RO, 2, 2))

SV_m <- SV_m |> dplyr::mutate(CHROM = gsub("^chr", "", CHROM))
SV_m <- SV_m |> dplyr::mutate(CHR2 = gsub("^chr", "", CHR2))

SV_d <- SV_d |> dplyr::mutate(CHROM = gsub("^chr", "", CHROM))
SV_d <- SV_d |> dplyr::mutate(CHR2 = gsub("^chr", "", CHR2))

CNA <- CNA |> dplyr::mutate(chr = gsub("^chr", "", chr))

# todo - see how to combine the mantra + delly results

all_visits <- unique(CNA$SAMPLE)
chromothripsis_mantra <- list()
chromothripsis_delly <- list()


for (i in all_visits) {
    CNA_tmp <- CNA |> dplyr::filter(SAMPLE == i)
    SV_d_tmp <- SV_d |> dplyr::filter(SAMPLE == i)
    SV_m_tmp <- SV_m |> dplyr::filter(SAMPLE == i)

    CN_data <- CNVsegs(
        chrom = as.character(CNA_tmp$chr),
        start = CNA_tmp$start,
        end = CNA_tmp$end,
        total_cn = CNA_tmp$copy.number
    )


    # remove Y chromosome?
    SV_d_tmp <- SV_d_tmp |> dplyr::filter(CHROM != "Y", CHR2 != "Y")

    SV_data_delly <- SVs(
        chrom1 = as.character(SV_d_tmp$CHROM),
        pos1 = as.numeric(SV_d_tmp$POS),
        chrom2 = as.character(SV_d_tmp$CHR2),
        pos2 = as.numeric(SV_d_tmp$ENDPOSSV),
        SVtype = as.character(SV_d_tmp$SVTYPE),
        strand1 = as.character(SV_d_tmp$RGENUPS_Strand), # Note - not standard outputs? RGENUPS_Strand	Strand_of_the_gene_upstream_of_the_left_(5')_breakpoint
        strand2 = as.character(SV_d_tmp$RGENDNS_Strand) # RGENDNS_Strand	Gene_Strand_of_the_gene_downstream_of_the_left_(5')_breakpoint
    )



    # SV_data_mantra <- SVs(
    #     chrom1 = as.character(SV_m_tmp$CHROM),
    #     pos1 = as.numeric(SV_m_tmp$POS),
    #     chrom2 = as.character(SV_m_tmp$CHR2),
    #     pos2 = as.numeric(SV_m_tmp$ENDPOSSV),
    #     SVtype = as.character(SV_m_tmp$SVTYPE),
    #     strand1 = as.character(SV_m_tmp$strand1), # Note - not standard outputs? RGENUPS_Strand	Strand_of_the_gene_upstream_of_the_left_(5')_breakpoint
    #     strand2 = as.character(SV_m_tmp$strand2) # RGENDNS_Strand	Gene_Strand_of_the_gene_downstream_of_the_left_(5')_breakpoint
    # )


    # chromothripsis_mantra[[i]] <- shatterseek( # Has invalid calls?
    #     SV.sample = SV_data_mantra,
    #     seg.sample = CN_data,
    #     genome = "hg19"
    # )

    chromothripsis_delly[[i]] <- shatterseek(
        SV.sample = SV_data_delly,
        seg.sample = CN_data,
        genome = "hg19"
    )
}
