these phyloseq objects are copied the https://github.com/ferreira-scm/Eimeria_AmpSeq repository, commit 091f1eb https://github.com/ferreira-scm/Eimeria_AmpSeq/commit/091f1ebdee2f2d3f3f8d53a3a1cb46deb8e830ce

preprocessing steps are https://github.com/ferreira-scm/Eimeria_AmpSeq/R


phyloseq object description:

PhyloSeqCombi_HMHZ_All.Rds: Multiamplicon dataset of faecal samples from wild-caught mice. All amplicon sequencing variants (ASVs) from all amplicons have been pooled per sample.

PhyloSeqList_HMHZ_All.Rds: Same but amplicons are not pooled per sample. resulting in a list of phyloseq objects for each amplicon.

PhyloSeqData18S_SILVA.Rds: Single amplicon datasets targeting 18S RNA gene from faecal samples of laboratory mice infected with Eimeria ferrisi. Each mice has been sampled each day during infection progression.

PhyloSeqList18S_SILVA.Rds: Same, but includes also 2 single amplicon datasets against 2 different regions of the 16S gene.

PhyloSeqData_All_Tax_New.Rds: Same samples but multi-amplicon has been used. ASVs from all amplicons are pooled per sample

PhyloSeqList_All_Tax_New.Rds: same as above but ASVs are not pooled, resulting in a list of phyloseq objects for each amplicon
