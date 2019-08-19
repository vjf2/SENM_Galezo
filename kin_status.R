#merge relatedness and kindat_pos

kin_peddata<-expected_kinship("dolphin_id", 
                 "mother_id",
                 "father_id", "sex", life_history_lookup)

#relatedness

comb_rel<-merge_pairs(kin_peddata, relatedness, 
                      "ID1", "ID2", "xID1", "ID2")

#If relatedness = 0 then 0, if relatednes > threshold then kin 

threshold<-0.0362

#remove animals not in analysis
comb_rel<-comb_rel[comb_rel$ID1 %in% dolphins &
                     comb_rel$ID2 %in% dolphins,]

comb_rel$kin_status<-ifelse(comb_rel$relatedness > threshold |
                              comb_rel$biparental > threshold, "kin", "non_kin")

comb_rel$kin_status<-ifelse(is.na(comb_rel$relatedness) &
                              comb_rel$biparental  < threshold, "unknown", comb_rel$kin_status)

#20% of pairs have a known relationship status

write.csv(comb_rel, "comb_rel.csv")
