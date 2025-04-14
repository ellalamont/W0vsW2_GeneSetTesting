# Making new gene set lists from the literature
# E. Lamont
# 4/1/25

# Once I save the rda files they should just load when I load Import_data.R and I shouldn't have to run this code again



source("Import_data.R")

###########################################################
######################## WALTER 2015 ######################
# Use the Walter 2015 sets: https://pmc.ncbi.nlm.nih.gov/articles/PMC4548467/

# Needs to be called allGeneSets so it is easier to load with all the others
allGeneSets <- list("Aerobic respiration" = c("Rv1542c", "Rv1854c", "Rv2193", "Rv2194", "Rv2195", "Rv2196", "Rv2249c", "Rv2470", "Rv3145", "Rv3146", "Rv3147", "Rv3148", "Rv3149", "Rv3150", "Rv3151", "Rv3152", "Rv3153", "Rv3154", "Rv3155", "Rv3156", "Rv3157", "Rv3158"),
                           "TCA cycle" = c("Rv0066c", "Rv0889c", "Rv0896", "Rv0951", "Rv0952", "Rv1098c", "Rv1131", "Rv1240", "Rv1248c", "Rv1475c", "Rv2215", "Rv2967c", "Rv3316", "Rv3318", "Rv3339c"),
                           "ATP synthesis" = c("Rv1304", "Rv1308"),
                           "Anaerobic respiration" = c("Rv0252", "Rv1161"),
                           "Glyoxylate bypass" = c("Rv0467", "Rv1837c", "Rv1915", "Rv1916"),
                           "Triacylglycerol synthases" = c("Rv3234c", "Rv3734c"),
                           "Ribosomal protein genes" = c("Rv0105c", "Rv0705", "Rv0706", "Rv0707", "Rv0708", "Rv0709", "Rv0710", "Rv0714", "Rv0716", "Rv0995", "Rv1015c", "Rv1298", "Rv1630", "Rv1642", "Rv2055c", "Rv2056c", "Rv2057c", "Rv2058c", "Rv2442c", "Rv2904c", "Rv2909c", "Rv3241c", "Rv3420c"),
                           "Polyketide synthases" = c("Rv1013", "Rv1180", "Rv1181", "Rv1182", "Rv1527c", "Rv1528c", "Rv1660", "Rv1661", "Rv1662", "Rv1663", "Rv1664", "Rv1665", "Rv2931", "Rv2932", "Rv2934", "Rv2935", "Rv2939", "Rv2940c", "Rv2946c", "Rv3800c", "Rv3820c"),
                           "ESX genes" = c("Rv0287", "Rv0288", "Rv1037c", "Rv1038c", "Rv1198", "Rv1792", "Rv1793", "Rv3019c", "Rv3444c", "Rv3445c", "Rv3904c", "Rv3905c"),
                           "Oxidative stress response" = c("Rv0338c", "Rv0816c", "Rv0848", "Rv1286", "Rv1317c", "Rv2386c", "Rv2737c", "Rv3202c", "Rv3585", "Rv3841", "Rv3914"),
                           "PE/PPE" = c("Rv0109", "Rv0152c", "Rv0159c", "Rv0160c", "Rv0256c", "Rv0285", "Rv0335c", "Rv0354c", "Rv0355c", "Rv0453", "Rv0578c", "Rv0755c", "Rv0916c", "Rv1039c", "Rv1087", "Rv1091", "Rv1196", "Rv1386", "Rv1430", "Rv1441c", "Rv1452c", "Rv1651c", "Rv1706c", "Rv1759c", "Rv1791", "Rv1800", "Rv1802", "Rv1803c", "Rv1809", "Rv1918c", "Rv1983", "Rv2108", "Rv2328", "Rv2353c", "Rv2356c", "Rv2396", "Rv2408", "Rv2431c", "Rv2519", "Rv2591", "Rv2615c", "Rv2634c", "Rv2741", "Rv2768c", "Rv2769c", "Rv3097c", "Rv3125c", "Rv3135", "Rv3350c", "Rv3426", "Rv3429", "Rv3477", "Rv3558", "Rv3622c", "Rv3652", "Rv3739c", "Rv3746c", "Rv3872", "Rv3873", "Rv3892c"),
                           "Phage proteins and insertion sequences" = c("Rv0094c", "Rv0336", "Rv0397", "Rv0605", "Rv0829", "Rv0920c", "Rv1034c", "Rv1036c", "Rv1047", "Rv1054", "Rv1128c", "Rv1148c", "Rv1313c", "Rv1573", "Rv1581c", "Rv1584c", "Rv1701", "Rv1702c", "Rv1765A", "Rv2013", "Rv2085", "Rv2087", "Rv2309c", "Rv2355", "Rv2647", "Rv2650c", "Rv2652c", "Rv2655c", "Rv2656c", "Rv2659c", "Rv2666", "Rv2792c", "Rv2810c", "Rv2943", "Rv2961", "Rv2979c", "Rv3349c", "Rv3386", "Rv3427c", "Rv3428c", "Rv3431c", "Rv3640c", "Rv3751", "Rv3770B"),
                           "Enduring hypoxic response" = c("Rv0048c", "Rv0084", "Rv0116c", "Rv0140", "Rv0171", "Rv0186", "Rv0233", "Rv0244c", "Rv0251c", "Rv0268c", "Rv0327c", "Rv0347", "Rv0384c", "Rv0520", "Rv0521", "Rv0551c", "Rv0563", "Rv0575c", "Rv0614", "Rv0654", "Rv0687", "Rv0726c", "Rv0766c", "Rv0767c", "Rv0793", "Rv0826", "Rv0847", "Rv0986", "Rv0989c", "Rv0991c", "Rv1048c", "Rv1082", "Rv1256c", "Rv1284", "Rv1375", "Rv1403c", "Rv1405c", "Rv1472", "Rv1623c", "Rv1628c", "Rv1652", "Rv1669", "Rv1755c", "Rv1760", "Rv1766", "Rv1874", "Rv1875", "Rv1888c", "Rv1891", "Rv1954c", "Rv1974", "Rv2012", "Rv2050", "Rv2122c", "Rv2299c", "Rv2389c", "Rv2390c", "Rv2465c", "Rv2497c", "Rv2499c", "Rv2500c", "Rv2504c", "Rv2557", "Rv2583c", "Rv2662", "Rv2663", "Rv2744c", "Rv2876", "Rv2913c", "Rv2964", "Rv3007c", "Rv3016", "Rv3161c", "Rv3176c", "Rv3206c", "Rv3221A", "Rv3406", "Rv3418c", "Rv3480c", "Rv3503c", "Rv3515c", "Rv3531c", "Rv3537", "Rv3571", "Rv3848", "Rv3864"),
                           "DosR regulon" = c("Rv0079", "Rv0080", "Rv0082", "Rv0083", "Rv0569", "Rv0570", "Rv0571c", "Rv0572c", "Rv0574c", "Rv1734c", "Rv1737c", "Rv1812c", "Rv1813c", "Rv1996", "Rv1997", "Rv2003c", "Rv2004c", "Rv2005c", "Rv2006", "Rv2007c", "Rv2028c", "Rv2029c", "Rv2030c", "Rv2032", "Rv2623", "Rv2624c", "Rv2625c", "Rv2629", "Rv2630", "Rv2631", "Rv3126c", "Rv3128c", "Rv3129", "Rv3131", "Rv3134c"),
                           "Sigma factors" = c("Rv0182c", "Rv0445c", "Rv0735", "Rv1189", "Rv2069", "Rv2703", "Rv3328c", "Rv3414c", "Rv3911"),
                           "Transcription factors" = c("Rv0014c", "Rv0015c", "Rv0018c", "Rv0019c", "Rv0020c", "Rv0022c", "Rv0042c", "Rv0043c", "Rv0078", "Rv0117", "Rv0135c", "Rv0158", "Rv0165c", "Rv0195", "Rv0212c", "Rv0238", "Rv0260c", "Rv0273c", "Rv0275c", "Rv0302", "Rv0339c", "Rv0386", "Rv0410c", "Rv0452", "Rv0465c", "Rv0490", "Rv0491", "Rv0586", "Rv0600c", "Rv0601c", "Rv0602c", "Rv0681", "Rv0691c", "Rv0737", "Rv0758", "Rv0818", "Rv0844c", "Rv0845", "Rv0880", "Rv0894", "Rv0902c", "Rv0903c", "Rv0981", "Rv0982", "Rv1019", "Rv1028c", "Rv1032c", "Rv1033c", "Rv1151c", "Rv1167c", "Rv1219c", "Rv1266c", "Rv1358", "Rv1359", "Rv1404", "Rv1423", "Rv1453", "Rv1473A", "Rv1474c", "Rv1479", "Rv1556", "Rv1626", "Rv1675c", "Rv1743", "Rv1776c", "Rv1816", "Rv1846c", "Rv2011c", "Rv2027c", "Rv2088", "Rv2160c", "Rv2175c", "Rv2232", "Rv2250c", "Rv2258c", "Rv2282c", "Rv2488c", "Rv2506", "Rv2711", "Rv2718c", "Rv2720", "Rv2779c", "Rv2884", "Rv2887", "Rv2914c", "Rv3050c", "Rv3058c", "Rv3060c", "Rv3095", "Rv3124", "Rv3143", "Rv3160c", "Rv3167c", "Rv3208", "Rv3219", "Rv3220c", "Rv3245c", "Rv3260c", "Rv3291c", "Rv3295", "Rv3405c", "Rv3416", "Rv3557c", "Rv3575c", "Rv3676", "Rv3736", "Rv3744", "Rv3764c", "Rv3765c", "Rv3830c", "Rv3833", "Rv3849", "Rv3855"),
                           "Toxins" = c("Rv0240", "Rv0299", "Rv0301", "Rv0549c", "Rv0582", "Rv0595c", "Rv0609", "Rv0617", "Rv0624", "Rv0627", "Rv0656c", "Rv0659c", "Rv0661c", "Rv0665", "Rv0749", "Rv0919", "Rv1102c", "Rv1114", "Rv1246c", "Rv1397c", "Rv1495", "Rv1561", "Rv1720c", "Rv1741", "Rv1838c", "Rv1942c", "Rv1953", "Rv1959c", "Rv1989c", "Rv2010", "Rv2035", "Rv2103c", "Rv2527", "Rv2530c", "Rv2546", "Rv2548", "Rv2549c", "Rv2596", "Rv2602", "Rv2757c", "Rv2759c", "Rv2866", "Rv2872", "Rv3189", "Rv3320c", "Rv3358", "Rv3384c", "Rv3408", "Rv3749c"),
                           "Antitoxins" = c("Rv0239", "Rv0298", "Rv0300", "Rv0550c", "Rv0581", "Rv0596c", "Rv0608", "Rv0623", "Rv0626", "Rv0657c", "Rv0660c", "Rv0662c", "Rv0664", "Rv0748", "Rv0918", "Rv1103c", "Rv1241", "Rv1247c", "Rv1398c", "Rv1494", "Rv1560", "Rv1721c", "Rv1839c", "Rv1943c", "Rv1952", "Rv2009", "Rv2018", "Rv2104c", "Rv2493", "Rv2530A", "Rv2547", "Rv2550c", "Rv2601A", "Rv2758c", "Rv2760c", "Rv2865", "Rv2871", "Rv3181c", "Rv3321c", "Rv3357", "Rv3385c", "Rv3407"),
                           "Stringent response" = c("Rv0009", "Rv0046c", "Rv0120c", "Rv0315", "Rv0462", "Rv0685", "Rv0803", "Rv0815c", "Rv0867c", "Rv1076", "Rv1269c", "Rv1274", "Rv1303", "Rv1368", "Rv1435c", "Rv1477", "Rv1483", "Rv1639c", "Rv1829", "Rv2163c", "Rv2169c", "Rv2213", "Rv2272", "Rv2406c", "Rv2468c", "Rv2633c", "Rv2638", "Rv2641", "Rv2903c", "Rv2945c", "Rv3029c", "Rv3048c", "Rv3083", "Rv3084", "Rv3400", "Rv3457c", "Rv3462c", "Rv3487c", "Rv3510c", "Rv3569c", "Rv3661", "Rv3762c", "Rv3763", "Rv3804c", "Rv3918c"),
                           "Isoniazid stress resonse" = c("Rv1592c", "Rv1772", "Rv2244", "Rv2247", "Rv3139"),
                           "Drug activating enzymes and drug targets" = c("Rv0005", "Rv0006", "Rv1484", "Rv2043c", "Rv3854c"),
                           "Putative drug efflux pumps and transporters" = c("Rv0783c", "Rv0933", "Rv1002c", "Rv1348", "Rv1456c", "Rv1634", "Rv1747", "Rv1877", "Rv2686c", "Rv2936", "Rv2937", "Rv2938", "Rv2994", "Rv3065", "Rv3679"),
                           "Cluster A" = c("Rv0362", "Rv0811c", "Rv1224", "Rv1268c", "Rv1717", "Rv1834", "Rv2832c", "Rv2863"),
                           "Cluster B" = c("Rv0221", "Rv0817c", "Rv1236", "Rv1273c", "Rv1371", "Rv1745c", "Rv1814", "Rv1817", "Rv1819c", "Rv1857", "Rv2609c", "Rv3896c"),
                           "Cluster D" = c("Rv1097c", "Rv1794", "Rv3197"),
                           "Cluster E" = c("Rv0470c", "Rv1117", "Rv1898", "Rv1925", "Rv2401A", "Rv2455c", "Rv2459", "Rv2562", "Rv2588c", "Rv2941", "Rv3528c", "Rv3614c", "Rv3615c", "Rv3616c", "Rv3687c", "Rv3724A", "Rv3865"),
                           "Non-specific stress responses" = c("Rv0191", "Rv0375c", "Rv0674", "Rv0838", "Rv0841", "Rv1044", "Rv1101c", "Rv1424c", "Rv1714", "Rv1718", "Rv1726", "Rv1735c", "Rv1774", "Rv1811", "Rv1961", "Rv1964", "Rv1969", "Rv1995", "Rv2038c", "Rv2156c", "Rv2501c", "Rv2643", "Rv2661c", "Rv2667", "Rv2742c", "Rv2836c", "Rv3288c", "Rv3289c", "Rv3450c"),
                    "ESX Genes 2" = c("Rv3017c", "Rv3019c", "Rv3904c", "Rv3905c", "Rv3445c", "Rv3444c", "Rv3020c", "Rv3891c", "Rv1198", "Rv1792", "Rv2347c", "Rv1197c", "Rv0288", "Rv0287", "Rv3890c", "Rv1037c", "Rv3619c", "Rv1793", "Rv3875", "Rv3874", "Rv3620c", "Rv1038c", "Rv2346c")
                           )

# #25 (cluster C) missing from the supplemental excel file! 

# SAVE AS RDA FOR LATER 
save(allGeneSets, file = "GeneSet_Data/Walter2015GeneSets.rda")

# ADDING TO GENE SET LIST OF LISTS 
# Want to add it to the allGeneSetList
allGeneSetList[["Walter2015GeneSets"]] <- Walter2015GeneSets


###########################################################
########################## ELLA ###########################
# Pulling from Honeyborne2016 and Cole1998
# https://pmc.ncbi.nlm.nih.gov/articles/PMC6851703/ for tca cycle genes
# Mycobrowsered "zinc" 

# Needs to be called allGeneSets so it is easier to load with all the others
allGeneSets <- list("glyoxylate bypass and methylcitrate cycle" = c("Rv0467", "Rv1915", "Rv1916", "Rv1837c", "Rv3323c", "Rv1131", "Rv1129c", "Rv1130"),
                    "ATP synthase" = c("Rv1308", "Rv1304", "Rv1311", "Rv1310", "Rv1305", "Rv1306", "Rv1309", "Rv1307"),
                    "Electron transport" = c("Rv0409", "Rv1623c", "Rv1622c", "Rv1620c", "Rv1621c", "Rv2007c", "Rv3554", "Rv1177", "Rv3503c", "Rv3029c", "Rv3028c", "Rv3106", "Rv0886", "Rv3251c", "Rv3250c"),
                    "TCA cycle" = c("Rv1475c", "Rv0889c", "Rv2498c", "Rv1098c", "Rv1131", "Rv0896", "Rv3339c", "Rv0066c", "Rv0794c", "Rv1240", "Rv2967c", "Rv3318", "Rv3319", "Rv3316", "Rv3317", "Rv1248c", "Rv2215", "Rv0951", "Rv0952"),
                    "Lipid biosynthesis" = c("Rv3285", "Rv0904c", "Rv3799c", "Rv3280", "Rv2247", "Rv2244", "Rv2523c", "Rv2243", "Rv0639", "Rv1483", "Rv1350", "Rv2002", "Rv0242c", "Rv2524c", "Rv1484", "Rv2245", "Rv2246", "Rv1618", "Rv2605c", "Rv0033", "Rv1344", "Rv1722", "Rv3221c", "Rv3472"),
                    "Zinc related" = c("Rv0198c", "Rv2359", "Rv2782c", "Rv3610c", "Rv0162c", "Rv0761c", "Rv3086"))

# SAVE AS RDA FOR LATER 
save(allGeneSets, file = "GeneSet_Data/EllaGeneSets.rda")









