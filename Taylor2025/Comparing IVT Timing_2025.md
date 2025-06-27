## Comparing IVT Timing Across Past Sequencing Experiments 


R1/C1 - August Overnight IVT:
summary(colSums(A))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 102886  178941  293916  336359  488810  712679
 summary(colSums(A>0))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  12120   14170   15079   15294   16618   18020

R1/C1 - November 4 HR IVT:
summary(colSums(A))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
   1561   22184   32576   46408   64799  144264
summary(colSums(A>0))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    969    6974    7958    7682    9799   12359

Single Tube NTS1 Overnight IVT (Same sequencing batch as November 4hr IVT RC1 Exp)
summary(colSums(A))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 171365  228084  290721  296514  341803  477746
summary(colSums(A > 0))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  14142   15873   16506   16458   17358   18271

  NTS1 4 Hr IVT (Same Sequencing Batch as November 4 hr IVT RC1 Experiment and Single Tube NTS1 Experiment)
summary(colSums(A))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
      0     325    3770    4373    6729   20093
summary(colSums(A>0))
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
      0     305    2706    2608    4119    7557