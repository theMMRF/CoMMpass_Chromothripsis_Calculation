# As mentioned here: https://ashpublications.org/blood/article/136/Supplement%201/52/471067/Copy-Number-Signatures-Predict-Chromothripsis-and
# Shatterseek seems to not predict chromothripsis well on the CoMMpass dataset

library(ShatterSeek)

# +: 3', -: 5'

# SAMPLE chr start end event copy.number bins median
# copy.number is what is needed
# CNA <- read.table("data/Copy_Number_Estimates/MMRF_CoMMpass_IA21_genome_ichorCNA.seg", header = T)
# CNAb <- read.table("data/Copy_Number_Estimates/MMRF_CoMMpass_IA21_genome_gatk_cna_PerGene_LargestOverlap.tsv", header = T)
CNAc <- read.table("data/Copy_Number_Estimates/MMRF_CoMMpass_IA21_genome_gatk_cna.seg", header = T)
CNAc <- CNAc |> dplyr::mutate(copy.number = 2 * 2^Segment_Mean)

# SAMPLE SVTYPE CHROM POS ID REF CHR2 ENDPOSSV CIPOS MAPQ PE SR SRQ CE CT IMPRECISE PRECISE HOMLEN RDRATIO LENSV RGENUPS_ENSG RGENUPS HUGO RGENUPS_Distance	RGENUPS_Strand	RGENISEC_ENSG	RGENISEC_HUGO	RGENISEC_Distance	RGENISEC_Strand	RGENDNS_ENSG	RGENDNS_HUGO	RGENDNS_Distance	RGENDNS_Strand	LGENUPS_ENSG	LGENUPS_HUGO	LGENUPS_Distance	LGENUPS_Strand	LGENISEC_ENSG	LGENISEC_HUGO	LGENISEC_Distance	LGENISEC_Strand	LGENDNS_ENSG	LGENDNS_HUGO	LGENDNS_Distance	LGENDNS_Strand	NORMAL_GT	TUMOR_GT	NORMAL_GL_RR	TUMOR_GL_RR	NORMAL_GL_RA	TUMOR_GL_RA	NORMAL_GL_AA	TUMOR_GL_AA	NORMAL_GQ	TUMOR_GQ	NORMAL_FT	TUMOR_FT	NORMAL_RC	TUMOR_RC	NORMAL_RCL	TUMOR_RCL	NORMAL_RCR	TUMOR_RCR	NORMAL_CN	TUMOR_CN	NORMAL_DR	TUMOR_DR	NORMAL_DV	TUMOR_DV	NORMAL_RR	TUMOR_RR	NORMAL_RV	TUMOR_RV	NORMAL_RCALT	TUMOR_RCALT	NORMAL_RDISTDISC1	TUMOR_RDISTDISC1	NORMAL_RDISTDISC2	TUMOR_RDISTDISC2	NORMAL_RCDIS1	TUMOR_RCDIS1	NORMAL_RCDIS2	TUMOR_RCDIS2	CONSENSUS	isBalTra	dupBalTraToRemove
SV_d <- read.table("data/Structural_Event_Files/MMRF_CoMMpass_IA21_genome_delly.tsv", header = T)

# SAMPLE	CHROM	POS	ID	REF	ALT	QUAL	FILTER	IMPRECISE	SVTYPE	SVLEN	END	CIPOS	CIEND	CIGAR	MATEID	EVENT	HOMLEN	HOMSEQ	SVINSLEN	SVINSSEQ	LEFT_SVINSSEQ	RIGHT_SVINSSEQ	BND_DEPTH	MATE_BND_DEPTH	SOMATIC	SOMATICSCORE	JUNCTION_SOMATICSCORE	TRA	INV	INVDUP	DUPTANDEM	CHR2	ENDPOSSV	RO	UNTESTED	RGENUPS	RGENISEC	RGENDNS	LGENUPS	LGENISEC	LGENDNS	GENELISTbtwBP	LENSV	NORMAL_PR	NORMAL_SR	NORMAL_RUT1	NORMAL_RUT2	NORMAL_RDISTDISC1	NORMAL_RDISTDISC2	NORMAL_RCDIS1	NORMAL_RCDIS2	NORMAL_DRNOISE1	NORMAL_DRNOISE2	NORMAL_EVALRCNOISE	NORMAL_CLMQ1	NORMAL_CLMQ2	TUMOR_PR	TUMOR_SR	TUMOR_RUT1	TUMOR_RUT2	TUMOR_RDISTDISC1	TUMOR_RDISTDISC2	TUMOR_RCDIS1	TUMOR_RCDIS2	TUMOR_DRNOISE1	TUMOR_DRNOISE2	TUMOR_EVALRCNOISE	TUMOR_CLMQ1	TUMOR_CLMQ2
SV_m <- read.table("data/Structural_Event_Files/MMRF_CoMMpass_IA21_genome_manta.tsv", header = T)

SV_m <- SV_m |> dplyr::mutate(strand1 = substr(RO, 1, 1), strand2 = substr(RO, 2, 2))

SV_m <- SV_m |> dplyr::mutate(CHROM = gsub("^chr", "", CHROM))
SV_m <- SV_m |> dplyr::mutate(CHR2 = gsub("^chr", "", CHR2))

SV_d <- SV_d |> dplyr::mutate(CHROM = gsub("^chr", "", CHROM))
SV_d <- SV_d |> dplyr::mutate(CHR2 = gsub("^chr", "", CHR2))

CNAc <- CNAc |> dplyr::mutate(Chromosome = gsub("^chr", "", Chromosome))
CNAc <- CNAc |> dplyr::filter(Chromosome != "Y") # shatterseek can't process Y chromosome
CNAc <- CNAc |> dplyr::mutate(copy.number = round(2 * 2^(as.numeric(Segment_Mean))))

CNAc <- CNAc |> dplyr::mutate(copy.number = dplyr::case_when(
    copy.number > 2 ~ 3,
    copy.number < 2 ~ 1,
    copy.number == 2 ~ 2,
    .default = -1
))

# Convert CT to +/- format
SV_d <- SV_d |> dplyr::mutate(strand1 = dplyr::case_when(
    CT == "3to3" ~ "+",
    CT == "3to5" ~ "+",
    CT == "5to3" ~ "-",
    CT == "5to5" ~ "-",
    .default = "ERROR"
))

SV_d <- SV_d |> dplyr::mutate(strand2 = dplyr::case_when(
    CT == "3to3" ~ "+",
    CT == "3to5" ~ "-",
    CT == "5to3" ~ "+",
    CT == "5to5" ~ "-",
    .default = "ERROR"
))

# Split INV to h2h and t2t
SV_d <- SV_d |> dplyr::mutate(SVTYPE = dplyr::case_when(
    CT == "3to3" & SVTYPE == "INV" ~ "h2hINV",
    CT == "5to5" & SVTYPE == "INV" ~ "t2tINV",
    .default = SVTYPE
))

# todo - see how to combine the mantra + delly results

all_visits <- unique(CNAc$SAMPLE)
chromothripsis_mantra <- list()
chromothripsis_delly <- list()



for (i in all_visits) {
    CNA_tmp <- CNAc |> dplyr::filter(SAMPLE == i)
    SV_d_tmp <- SV_d |> dplyr::filter(SAMPLE == i)
    SV_m_tmp <- SV_m |> dplyr::filter(SAMPLE == i)

    # CN_data <- CNVsegs(
    #     chrom = as.character(CNA_tmp$chr),
    #     start = CNA_tmp$start,
    #     end = CNA_tmp$end,
    #     total_cn = CNA_tmp$copy.number
    # )


    CN_data <- CNVsegs(
        chrom = as.character(CNA_tmp$Chromosome),
        start = CNA_tmp$Start,
        end = CNA_tmp$End,
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
        strand1 = as.character(SV_d_tmp$strand1), # Note - not standard outputs? RGENUPS_Strand	Strand_of_the_gene_upstream_of_the_left_(5')_breakpoint
        strand2 = as.character(SV_d_tmp$strand2) # RGENDNS_Strand	Gene_Strand_of_the_gene_downstream_of_the_left_(5')_breakpoint
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

# Criteria: https://github.com/parklab/ShatterSeek/issues/7
# https://github.com/parklab/ShatterSeek/issues/22
# pval_fragment_joins > 0.05
# chr_breakpoint_enrichment < 0.05
# pval_exp_cluster < 0.05
# number_DEL + number_DUP + number_h2hINV + number_t2tINV >=3
# number_TRA>=4
# max_number_oscillating_CN_segments_2_states>=7

HC1_calls <- c()
HC2_calls <- c()
HC3_calls <- c()
LC_calls <- c()

chromSum_ind <- list()
for (sample in names(chromothripsis_delly)) {
    chromSum <- chromothripsis_delly[[sample]]@chromSummary
    chromSum <- chromSum |> dplyr::mutate(
        number_intrachromosomal_SVs = number_DEL + number_DUP + number_h2hINV + number_t2tINV
    )
    # 1st High Confidence filter
    # filt1 = all_data$number_intrachromosomal_SVs >= 6
    # filt2 = all_data$max_number_oscillating_CN_segments_2_states >= 7
    # filt3 = all_data$pval_fragment_joins >= 0.05
    # filt4 = (all_data$chr_breakpoint_enrichment <= 0.05) | (all_data$pval_exp_chr <= 0.05)

    # 2nd High Confidence filter
    # filt1 = all_data$number_intrachromosomal_SVs >= 3
    # filt2 = all_data$number_TRA >= 4
    # filt3 = all_data$max_number_oscillating_CN_segments_2_states >= 7
    # filt4 = all_data$pval_fragment_joins >= 0.05
    # HC2 = (filt1) & (filt2) & (filt3) & (filt4)
    # HC2[is.na(HC2)] <- FALSE
    # all_data$HC2 <- HC2

    chromSum <- chromSum |> dplyr::mutate(
        HC1 = (
            pval_fragment_joins > 0.05 &
                number_intrachromosomal_SVs >= 6 &
                max_number_oscillating_CN_segments_2_states >= 7 &
                (chr_breakpoint_enrichment < 0.05 | pval_exp_chr < 0.05)
        ), HC2 = (
            number_intrachromosomal_SVs >= 3 &
                number_TRA >= 4 &
                max_number_oscillating_CN_segments_2_states > 7 &
                pval_fragment_joins > 0.05
        ),
        HC3 = (
            clusterSize_including_TRA >= 40 &
                pval_fragment_joins >= 0.05
        ),
        LC = (
            number_intrachromosomal_SVs >= 6 &
                (max_number_oscillating_CN_segments_2_states >= 4 | max_number_oscillating_CN_segments_2_states <= 6) &
                pval_fragment_joins >= 0.05 &
                (chr_breakpoint_enrichment <= 0.05 | pval_exp_chr <= 0.05)
        )
    )

    chromSum_ind[[sample]] <- chromSum
    HC1_calls <- c(HC1_calls, any(chromSum$HC1, na.rm = T))
    HC2_calls <- c(HC2_calls, any(chromSum$HC2, na.rm = T))
    HC3_calls <- c(HC3_calls, any(chromSum$HC3, na.rm = T))
    LC_calls <- c(LC_calls, any(chromSum$LC, na.rm = T))
}
