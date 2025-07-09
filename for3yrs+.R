######################################################################################
### DATA PREPARATION
# Rename variables: height → Ht, sbp → msbp, dbp → mdbp
# Truncate age and name as age_trun, categorize sex: 1=male, 2=female
# Download 'ht.lms.csv' and 'bp.refcn.csv' to your working directory
######################################################################################

# Load required library
library(data.table)

# 1. Read reference data and convert to data.table
ht.ref <- fread('ht.lms.csv')
bp.ref <- fread('bp.refcn.csv')
setDT(data)  # Convert original data to data.table

# Merge with height reference data
d1 <- merge(data, ht.ref, by = c("age_trun", "sex"), all.x = TRUE)
setDT(d1)  # Ensure d1 is data.table

# 2. Calculate extreme height boundaries (in-place modification)
d1[, `:=`(
    sd2ps = ht_m * ((1 + 2 * ht_l * ht_s) ^ (1 / ht_l)),
    sd3ps = ht_m * ((1 + 3 * ht_l * ht_s) ^ (1 / ht_l)),
    sd3ng = ht_m * ((1 - 3 * ht_l * ht_s) ^ (1 / ht_l)),
    sd2ng = ht_m * ((1 - 2 * ht_l * ht_s) ^ (1 / ht_l))
)]

# 3. Calculate height z-scores (split into two steps to avoid dependency issues)
d1[, zht_raw := ((Ht / ht_m) ^ ht_l - 1) / (ht_s * ht_l)]  # First compute raw z-score
d1[, zht := fcase(                                         # Then adjust extremes
    zht_raw > 3, 3 + (Ht - sd3ps) / (sd3ps - sd2ps),
    zht_raw < -3, -3 - (sd3ng - Ht) / (sd2ng - sd3ng),
    default = zht_raw
)]

# 4. Convert z-score to height percentiles
d1[, htpct := fcase(
    zht < -1.4632, 5,
    zht < -0.97802, 10,
    zht < -0.33724, 25,
    zht < 0.33724, 50,
    zht < 0.97802, 75,
    zht < 1.4632, 90,
    zht >= 1.4632, 95,
    default = NA
)]

# 5. Merge with blood pressure reference data
d2 <- merge(d1, bp.ref, by = c("age_trun", "sex", "htpct"), all.x = TRUE)
setDT(d2)

# 6. Blood pressure grading
# 6.1 Childhood BP grading
d2[, hbp_ch := fcase(
    msbp < sbp90 & mdbp < dbp90, 1,
    (sbp90 <= msbp & msbp < sbp95) | (dbp90 <= mdbp & mdbp < dbp95) | 
        (sbp90 >= 120 & msbp < sbp95) | (mdbp >= 80 & mdbp < dbp95), 2,
    (sbp95 <= msbp & msbp < sbp99 + 5) | (dbp95 <= mdbp & mdbp < dbp99 + 5), 3,
    msbp >= sbp99 + 5 | mdbp >= dbp99 + 5, 4,
    default = NA
)]

# 6.2 Adolescent BP grading
d2[, hbp_ad := fcase(
    msbp < 120 & mdbp < 80, 1,
    (120 <= msbp & msbp < 140) | (80 <= mdbp & mdbp < 90), 2,
    (140 <= msbp & msbp < 160) | (90 <= mdbp & mdbp < 100), 3,
    (160 <= msbp & msbp < 180) | (100 <= mdbp & mdbp < 110), 4,
    msbp >= 180 | mdbp >= 110, 5,
    default = NA  
)]

# 7. Determine final BP grade (age-specific)
d2[, bp_grade := fifelse(age_trun < 16, hbp_ch, hbp_ad)]
