import numpy as np
import matplotlib.pyplot as plt

import vireoSNP
print("vireoSNP version: %s" %vireoSNP.__version__)

res = vireoSNP.vcf.match_VCF_samples('malaria_scRNA_kai/EXP1_l1/GT_donors.vireo.vcf.gz',
                                     'malaria_scRNA_kai/EXP1_l2/GT_donors.vireo.vcf.gz',
                                     GT_tag1 = 'PL',
                                     GT_tag2 = 'PL')

#Shape for Geno Prob in VCF1: (117813, 2, 3)
#Shape for Geno Prob in VCF2: (117813, 2, 3)
#n_variants in VCF1, VCF2 and matched: 117813, 27360, 17664
#aligned donors:
#['donor0' 'donor1']
#['donor1' 'donor0']

res = vireoSNP.vcf.match_VCF_samples('malaria_scRNA_kai/EXP2_l1/GT_donors.vireo.vcf.gz',
                                     'malaria_scRNA_kai/EXP2_l2_reseq/GT_donors.vireo.vcf.gz',
                                     GT_tag1 = 'PL',
                                     GT_tag2 = 'PL')
#Shape for Geno Prob in VCF1: (59229, 3, 3)
#Shape for Geno Prob in VCF2: (59229, 3, 3)
#n_variants in VCF1, VCF2 and matched: 59229, 27360, 24896
#aligned donors:
#['donor0' 'donor1' 'donor2']
#['donor1' 'donor2' 'donor0']

res = vireoSNP.vcf.match_VCF_samples('malaria_scRNA_kai/EXP2_l1/GT_donors.vireo.vcf.gz',
                                     'malaria_scRNA_kai/EXP2_l3/GT_donors.vireo.vcf.gz',
                                     GT_tag1 = 'PL',
                                     GT_tag2 = 'PL')

#Shape for Geno Prob in VCF1: (59229, 3, 3)
#Shape for Geno Prob in VCF2: (59229, 3, 3)
#n_variants in VCF1, VCF2 and matched: 59229, 59618, 47587
#aligned donors:
#['donor0' 'donor1' 'donor2']
#['donor2' 'donor1' 'donor0']
fig = plt.figure()
vireoSNP.plot.heat_matrix(res['matched_GPb_diff'],
                          res['matched_donors1'] ,
                          res['matched_donors2'] )
plt.title("Geno Prob Delta: %d SNPs" %(res['matched_n_var']))
plt.tight_layout()
plt.show()
