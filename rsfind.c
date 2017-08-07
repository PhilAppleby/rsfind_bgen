#include <sqlite3.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <zlib.h>
#include <math.h>
#include "rsfind.h"
/*
 * What this needs to do (only):
 * - Accept cmd line args genofile, dbfile, and rsid (both strings)
 * - Access sqlitedb using the dbfile path: obtain connection 
 *
 *
 */

static struct option long_options[] = {
  {"genofile",  required_argument, 0, 'g'},
  {"dbfile",    required_argument, 0, 'd'},
  {"rsid",      required_argument, 0, 'r'},
  {0, 0, 0, 0}
};
// Local functions and data
static char *parse_bgen_header_data(FILE * fPtr);
static unsigned int parse_bgen_variant_buf(void  *vBuf, long int bSize);
static int inflate_data(const void *src, int srcLen, void *dst, int dstLen);
static int ipow(int base, int exp);

static unsigned int probLen;
static void * bufPtr;
static void *ucompDest;
static unsigned int vremLen;
static unsigned char B;
static unsigned short int numAlleles;
static unsigned int N;

int main(int argc, char **argv) {

  char *genoFilePath;
  char *dbFilePath;
  char *rsid;

  int c;

  while(1) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "d:",
                     long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1) {
      fprintf(stderr, "End of cmd line opts\n");
      break; 
    }

    printf ("switch stmt\n");
    switch (c) {
      case 0:
        break;
      case 'g':
        genoFilePath = optarg;
        break;
      case 'd':
        dbFilePath = optarg;
        break;
      case 'r':
        rsid = optarg;
        break;
      default:
        abort();
    }
  }
  /* Preamble - get some size metrics */

  printf("Short Int %d\n", sizeof(short int));
  printf("Int %d\n", sizeof(int));
  printf("UInt %d\n", sizeof(unsigned int));
  printf("Long Int %d\n", sizeof(long int));

  probLen = init_rsid_search(genoFilePath, dbFilePath, rsid);

  float *probBuff = malloc(sizeof(float) * probLen);

  int rtn = get_bgen_genotype_probs(probBuff);

  free(probBuff);
    
  return 0;
}
/*---------------------------------------------------------------- 
 *
 *
 *
 ----------------------------------------------------------------*/
unsigned int init_rsid_search(char * genoFilePath, char * dbFilePath, char * rsid) {    
  sqlite3 *db;
  char *err_msg = 0;
  sqlite3_stmt *res;
    
  /* 1) sqlite lookup phase */
  int rc = sqlite3_open(dbFilePath, &db);
    
  if (rc != SQLITE_OK) {
    fprintf(stderr, "Cannot open database: %s\n", sqlite3_errmsg(db));
    sqlite3_close(db);
    return (0);
  }
  printf ("Opened %s successfully\n", dbFilePath);
  char *sql = "SELECT file_start_position, size_in_bytes from Variant where rsid=?";

  rc = sqlite3_prepare_v2(db, sql, -1, &res, 0);

  if (rc == SQLITE_OK) {
    sqlite3_bind_text(res, 1, rsid, strlen(rsid),SQLITE_STATIC);
  } else {
    fprintf(stderr, "Failed to execute statement: %s\n", sqlite3_errmsg(db));
    return (0);
  }
      
  int step = sqlite3_step(res);
  long int fOffset = 0L; // file offset
  long int bSize = 0L;   // buffer size
          
  if (step == SQLITE_ROW) {
    fOffset = atol(sqlite3_column_text(res, 0));
    printf("%ld,", fOffset);
    bSize = atol(sqlite3_column_text(res, 1));
    printf("%ld\n", bSize);
  }
  else {
    fprintf(stderr, "rsid not found: %s\n", rsid);
    return (0);
  }
  sqlite3_finalize(res);
  sqlite3_close(db);
  /* 2) geno file handling phase */
  void *fBuf = malloc(bSize);

  FILE *fPtr = fopen(genoFilePath, "rb");
  char * hdr_data = parse_bgen_header_data(fPtr);

  fseek(fPtr, fOffset, SEEK_SET);
  size_t ret_code = fread(fBuf, 1, bSize, fPtr);
  if (ret_code != bSize) {
    fprintf(stderr, "Read Error - %ld vs %ld\n", ret_code, bSize);
    return (0);
  } else {
    fprintf(stderr, "Read OK - %ld vs %ld\n", ret_code, bSize);
  }
  fclose(fPtr);
  /* 3) results parsing phase */
  return (parse_bgen_variant_buf(fBuf, bSize));
}
/*---------------------------------------------------------------- 
 *
 *
 *
 ----------------------------------------------------------------*/
char *parse_bgen_header_data(FILE *fPtr) {
  unsigned int inBuff[256];  
  char magicBytes[5];  

  fread(inBuff, 1, 4, fPtr);
  unsigned int *bufPtr = inBuff;
  printf("Offset of var data start (add 4) for actual file offset %d\n", *bufPtr);

  fread(inBuff, 1, 16, fPtr);
  bufPtr = inBuff;
  unsigned int hLen = *bufPtr;
  printf("Len hdr data %d\n", *bufPtr);
  bufPtr++;
  printf("Num variants in file %d\n", *bufPtr);
  bufPtr++;
  printf("Num samples in samples in file %d\n", *bufPtr);
  bufPtr++;
  strncpy(magicBytes, (char *)bufPtr, 4);
  printf("Magic bytes %s\n", magicBytes);
  /* Skip special area */
  fread(inBuff, 1, hLen-20, fPtr);
  /* flags */
  fread(inBuff, 1, 4, fPtr);
  bufPtr = inBuff;
  unsigned int n = *bufPtr;
  while (n) {
    if (n & 1)
      fprintf(stderr, "1");
    else
      fprintf(stderr, "0");
    n >>= 1;
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "Hexval %x\n", *bufPtr);

  return ((char *)NULL);
}
/*---------------------------------------------------------------- 
 *
 *
 *
 ----------------------------------------------------------------*/
unsigned int parse_bgen_variant_buf(void *vBuf, long int bSize) {
  bufPtr = vBuf;
  // variant id 
  unsigned short int vIdLen = *(unsigned short int *)bufPtr;
  fprintf(stderr, "varId Len %d\n", vIdLen);  
  char *varId = (char *)malloc(vIdLen + 1);
  bufPtr += sizeof(unsigned short int);
  strncpy(varId, bufPtr, vIdLen);
  fprintf(stderr, "VarId %s\n", varId);  
  bufPtr += vIdLen;
  // end variant id
  // rs id 
  unsigned short int rsidLen = *(unsigned short int *)bufPtr;
  fprintf(stderr, "rsid Len %d\n", rsidLen);  
  char *rsid = (char *)malloc(rsidLen + 1);
  rsid[rsidLen] = '\0';
  bufPtr += sizeof(unsigned short int);
  strncpy(rsid, (char *)bufPtr, rsidLen);
  fprintf(stderr, "rsid %s\n", rsid);  
  bufPtr += rsidLen;
  // end variant id
  // chromosome
  unsigned short int chrLen = *(unsigned short int *)bufPtr;
  fprintf(stderr, "chr Len %d\n", chrLen);  
  char *chr = (char *)malloc(chrLen + 1);
  chr[chrLen] = '\0';
  bufPtr += sizeof(unsigned short int);
  strncpy(chr, bufPtr, chrLen);
  fprintf(stderr, "chr %s\n", chr);  
  bufPtr += chrLen;
  // end chromosome
  // posn
  unsigned int posn = *(unsigned int *)bufPtr;
  fprintf(stderr, "var posn %d\n", posn);  
  bufPtr += sizeof(unsigned int);
  // end posn
  // num alleles
  unsigned short int k = *(unsigned short int *)bufPtr;
  fprintf(stderr, "K=%d\n", k);  
  bufPtr += sizeof(unsigned short int);
  // end posn
  // for each allele
  for (int i = 0; i < k; i++) {
    unsigned int aLen = *(unsigned int *)bufPtr;
    char *allele = (char *)malloc(aLen + 1);
    allele[aLen] = '\0';
    bufPtr += sizeof(unsigned int);
    strncpy(allele, bufPtr, aLen);
    fprintf(stderr, "allele=%s, num=%d\n", allele, (i+1));  
    free(allele);
    bufPtr += aLen;
  }
  // end for
  //
  free(varId);
  free(rsid);
  free(chr);
  // The genotype part of the block
  unsigned int vdataLen = *(unsigned int *)bufPtr;
  fprintf(stderr, "vdata Len %d\n", vdataLen);  
  bufPtr += sizeof(unsigned int);
  unsigned int ucompLen = *(unsigned int *)bufPtr;
  fprintf(stderr, "ucomp Len %d\n", ucompLen);  
  bufPtr += sizeof(unsigned int);

  ucompDest = malloc(ucompLen);
  vremLen = ucompLen;

  int ret_count = inflate_data(bufPtr, (vdataLen - 4), ucompDest, ucompLen);
  fprintf(stderr, "return from inflate %d\n", ret_count);  
  bufPtr = ucompDest;
  N = *(unsigned int *)bufPtr;
  fprintf(stderr, "N=%d\n", N);  
  bufPtr += sizeof(unsigned int);
  vremLen -= sizeof(unsigned int);
  numAlleles = *(unsigned short int *)bufPtr;
  fprintf(stderr, "numAlleles=%d\n", numAlleles);  
  bufPtr += sizeof(unsigned short int);
  vremLen -= sizeof(unsigned short int);
  unsigned char minPloidy = *(unsigned char *)bufPtr;
  fprintf(stderr, "minPloidy=%d\n", minPloidy);  
  bufPtr += sizeof(unsigned char);
  vremLen -= sizeof(unsigned char);
  unsigned char maxPloidy = *(unsigned char *)bufPtr;
  fprintf(stderr, "maxPloidy=%d\n", maxPloidy);  
  bufPtr += sizeof(unsigned char);
  vremLen -= sizeof(unsigned char);
  char * sampleFlags = (char *)bufPtr;
  bufPtr += (sizeof(char) * N);
  vremLen -= (sizeof(char) * N);
  unsigned char phased = *(unsigned char *)bufPtr;
  fprintf(stderr, "phased=%d\n", phased);  
  bufPtr += sizeof(unsigned char);
  vremLen -= sizeof(unsigned char);
  B = *(unsigned char *)bufPtr;
  fprintf(stderr, "B=%d\n", B);  
  bufPtr += sizeof(unsigned char);
  vremLen -= sizeof(unsigned char);
  return(N*3); // number of floats required
}

unsigned int get_bgen_genotype_probs(float *probResults) {
  fprintf(stderr, "remaining length=%d, B=%d, numAlleles=%d\n", vremLen, B, numAlleles);  
  // float flDiv = powf((float)2,(float)B);
  float flDiv = ipow(2, B) - 1;
  unsigned int ctr = 0;
  // we're seriously in diploid mode here (should use an array of probs of lenth numAlles - 1)
  float prob1;
  int homrefctr = 0;
  unsigned int l = 0;
  for (int i = vremLen; i > 0; i--) {
    unsigned char iByte = *(unsigned char *)bufPtr;
    float prob = (float)iByte / flDiv;
    if (ctr % numAlleles == 0) {
      prob1 = prob;
    } else {
      //printf("%.3f\n", prb);
      //printf("%g,%g,%g\n", prob1, prob, (1 - (prob1 + prob)));
      probResults[l] = prob1;
      probResults[l + 1] = prob;
      probResults[l + 2] = (1 - (prob1 + prob));
      l += 3;
      if ((prob1 == 1.0) && (prob == 0.0)) {
        homrefctr++;
      }
    }
    ctr++;
    bufPtr += sizeof(unsigned char);
  }

  fprintf(stderr, "homref=%d, N=%d, prop=%.3f\n", homrefctr, N, ((float)homrefctr/(float)N));

  free(ucompDest); // dangerous to put this in a separate call
  return (l);
}
//
// Wrap zlib inflate for a one-off buffer
//
// Lifted from code found on stack overflow:
// https://stackoverflow.com/questions/4901842/in-memory-decompression-with-zlib
//
//
int inflate_data(const void *src, int srcLen, void *dst, int dstLen) {
  z_stream strm  = {0};
  strm.total_in  = strm.avail_in  = srcLen;
  strm.total_out = strm.avail_out = dstLen;
  strm.next_in   = (Bytef *) src;
  strm.next_out  = (Bytef *) dst;

  strm.zalloc = Z_NULL;
  strm.zfree  = Z_NULL;
  strm.opaque = Z_NULL;

  int err = -1;
  int ret = -1;

  err = inflateInit2(&strm, (15 + 32)); //15 window bits, and the +32 tells zlib to to detect if using gzip or zlib
  if (err == Z_OK) {
    err = inflate(&strm, Z_FINISH);
    if (err == Z_STREAM_END) {
      ret = strm.total_out;
    }
    else {
      inflateEnd(&strm);
      return err;
    }
  }
  else {
    inflateEnd(&strm);
    return err;
  }

  inflateEnd(&strm);
  return ret;
}
//
// Raise an integer to a power
//
// Lifted from code found on stack overflow:
// https://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int
//
//
int ipow(int base, int exp)
{
  int result = 1;
  while (exp) {
    if (exp & 1)
      result *= base;
    exp >>= 1;
    base *= base;
  }
  return result;
}
