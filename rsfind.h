#ifndef RSFIND_H
#define RSFIND_H

unsigned int init_rsid_search(char * genoFilePath, char * dbFilePath, char * rsid);
unsigned int get_bgen_genotype_probs(float * res);

#endif /* RSFIND_H */
