



#ifndef BAM_BAM_H
#define BAM_BAM_H



#define BAM_VERSION "0.1.18 (r982:295)"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifndef BAM_LITE
#define BAM_VIRTUAL_OFFSET16

#include "bgzf.h"


typedef BGZF *bamFile;
#define bam_open(fn, mode) bgzf_open(fn, mode)
#define bam_dopen(fd, mode) bgzf_fdopen(fd, mode)
#define bam_close(fp) bgzf_close(fp)
#define bam_read(fp, buf, size) bgzf_read(fp, buf, size)
#define bam_write(fp, buf, size) bgzf_write(fp, buf, size)
#define bam_tell(fp) bgzf_tell(fp)
#define bam_seek(fp, pos, dir) bgzf_seek(fp, pos, dir)
#else
#define BAM_TRUE_OFFSET
#include <zlib.h>
typedef gzFile bamFile;
#define bam_open(fn, mode) gzopen(fn, mode)
#define bam_dopen(fd, mode) gzdopen(fd, mode)
#define bam_close(fp) gzclose(fp)
#define bam_read(fp, buf, size) gzread(fp, buf, size)

#endif


typedef struct
{
  int32_t n_targets;
  char **target_name;
  uint32_t *target_len;
  void *dict, *hash, *rg2lib;
  size_t l_text, n_text;
  char *text;
} bam_header_t;


#define BAM_FPAIRED        1

#define BAM_FPROPER_PAIR   2

#define BAM_FUNMAP         4

#define BAM_FMUNMAP        8

#define BAM_FREVERSE      16

#define BAM_FMREVERSE     32

#define BAM_FREAD1        64

#define BAM_FREAD2       128

#define BAM_FSECONDARY   256

#define BAM_FQCFAIL      512

#define BAM_FDUP        1024

#define BAM_OFDEC          0
#define BAM_OFHEX          1
#define BAM_OFSTR          2


#define BAM_DEF_MASK (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)

#define BAM_CORE_SIZE   sizeof(bam1_core_t)


#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  ((1 << BAM_CIGAR_SHIFT) - 1)



#define BAM_CMATCH      0

#define BAM_CINS        1

#define BAM_CDEL        2

#define BAM_CREF_SKIP   3

#define BAM_CSOFT_CLIP  4

#define BAM_CHARD_CLIP  5

#define BAM_CPAD        6

#define BAM_CEQUAL        7

#define BAM_CDIFF        8


typedef struct
{
  int32_t tid;
  int32_t pos;
  uint32_t bin: 16, qual: 8, l_qname: 8;
  uint32_t flag: 16, n_cigar: 16;
  int32_t l_qseq;
  int32_t mtid;
  int32_t mpos;
  int32_t isize;
} bam1_core_t;


typedef struct
{
  bam1_core_t core;
  int l_aux, data_len, m_data;
  uint8_t *data;
} bam1_t;

typedef struct __bam_iter_t *bam_iter_t;

#define bam1_strand(b) (((b)->core.flag&BAM_FREVERSE) != 0)
#define bam1_mstrand(b) (((b)->core.flag&BAM_FMREVERSE) != 0)


#define bam1_cigar(b) ((uint32_t*)((b)->data + (b)->core.l_qname))


#define bam1_qname(b) ((char*)((b)->data))


#define bam1_seq(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname)


#define bam1_qual(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1))


#define bam1_seqi(s, i) ((s)[(i)/2] >> 4*(1-(i)%2) & 0xf)


#define bam1_aux(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname + (b)->core.l_qseq + ((b)->core.l_qseq + 1)/2)

#ifndef kroundup32

#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif


extern int bam_is_be;


extern int bam_verbose;


extern unsigned char bam_nt16_table[256];


extern char *bam_nt16_rev_table;

extern char bam_nt16_nt4_table[];

#ifdef __cplusplus
extern "C" {
#endif




typedef struct __tamFile_t *tamFile;


tamFile sam_open(const char *fn);


void sam_close(tamFile fp);


int sam_read1(tamFile fp, bam_header_t *header, bam1_t *b);


bam_header_t *sam_header_read2(const char *fn_list);


bam_header_t *sam_header_read(tamFile fp);


int sam_header_parse(bam_header_t *h);
int32_t bam_get_tid(const bam_header_t *header, const char *seq_name);


int sam_header_parse_rg(bam_header_t *h);

#define sam_write1(header, b) bam_view1(header, b)




int bam_strmap_put(void *strmap, const char *rg, const char *lib);
const char *bam_strmap_get(const void *strmap, const char *rg);
void *bam_strmap_dup(const void *);
void *bam_strmap_init();
void bam_strmap_destroy(void *strmap);





bam_header_t *bam_header_init();


void bam_header_destroy(bam_header_t *header);


bam_header_t *bam_header_read(bamFile fp);


int bam_header_write(bamFile fp, const bam_header_t *header);


int bam_read1(bamFile fp, bam1_t *b);


int bam_write1_core(bamFile fp, const bam1_core_t *c, int data_len, uint8_t *data);


int bam_write1(bamFile fp, const bam1_t *b);


#define bam_init1() ((bam1_t*)calloc(1, sizeof(bam1_t)))


#define bam_destroy1(b) do {          \
      if (b) { free((b)->data); free(b); }  \
    } while (0)


char *bam_format1(const bam_header_t *header, const bam1_t *b);

char *bam_format1_core(const bam_header_t *header, const bam1_t *b, int of);


int bam_validate1(const bam_header_t *header, const bam1_t *b);

const char *bam_get_library(bam_header_t *header, const bam1_t *b);





typedef struct
{
  bam1_t *b;
  int32_t qpos;
  int indel, level;
  uint32_t is_del: 1, is_head: 1, is_tail: 1, is_refskip: 1, aux: 28;
} bam_pileup1_t;

typedef int (*bam_plp_auto_f)(void *data, bam1_t *b);

struct __bam_plp_t;
typedef struct __bam_plp_t *bam_plp_t;

bam_plp_t bam_plp_init(bam_plp_auto_f func, void *data);
int bam_plp_push(bam_plp_t iter, const bam1_t *b);
const bam_pileup1_t *bam_plp_next(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp);
const bam_pileup1_t *bam_plp_auto(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp);
void bam_plp_set_mask(bam_plp_t iter, int mask);
void bam_plp_set_maxcnt(bam_plp_t iter, int maxcnt);
void bam_plp_reset(bam_plp_t iter);
void bam_plp_destroy(bam_plp_t iter);

struct __bam_mplp_t;
typedef struct __bam_mplp_t *bam_mplp_t;

bam_mplp_t bam_mplp_init(int n, bam_plp_auto_f func, void **data);
void bam_mplp_destroy(bam_mplp_t iter);
void bam_mplp_set_maxcnt(bam_mplp_t iter, int maxcnt);
int bam_mplp_auto(bam_mplp_t iter, int *_tid, int *_pos, int *n_plp, const bam_pileup1_t **plp);


typedef int (*bam_pileup_f)(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data);

typedef struct
{
  bam_plp_t iter;
  bam_pileup_f func;
  void *data;
} bam_plbuf_t;

void bam_plbuf_set_mask(bam_plbuf_t *buf, int mask);
void bam_plbuf_reset(bam_plbuf_t *buf);
bam_plbuf_t *bam_plbuf_init(bam_pileup_f func, void *data);
void bam_plbuf_destroy(bam_plbuf_t *buf);
int bam_plbuf_push(const bam1_t *b, bam_plbuf_t *buf);

int bam_pileup_file(bamFile fp, int mask, bam_pileup_f func, void *func_data);

struct __bam_lplbuf_t;
typedef struct __bam_lplbuf_t bam_lplbuf_t;

void bam_lplbuf_reset(bam_lplbuf_t *buf);


bam_lplbuf_t *bam_lplbuf_init(bam_pileup_f func, void *data);


void bam_lplbuf_destroy(bam_lplbuf_t *tv);


int bam_lplbuf_push(const bam1_t *b, bam_lplbuf_t *buf);




struct __bam_index_t;
typedef struct __bam_index_t bam_index_t;


int bam_index_build(const char *fn);


bam_index_t *bam_index_load(const char *fn);


void bam_index_destroy(bam_index_t *idx);


typedef int (*bam_fetch_f)(const bam1_t *b, void *data);


int bam_fetch(bamFile fp, const bam_index_t *idx, int tid, int beg, int end, void *data, bam_fetch_f func);

bam_iter_t bam_iter_query(const bam_index_t *idx, int tid, int beg, int end);
int bam_iter_read(bamFile fp, bam_iter_t iter, bam1_t *b);
void bam_iter_destroy(bam_iter_t iter);


int bam_parse_region(bam_header_t *header, const char *str, int *ref_id, int *begin, int *end);





uint8_t *bam_aux_get(const bam1_t *b, const char tag[2]);

int32_t bam_aux2i(const uint8_t *s);
float bam_aux2f(const uint8_t *s);
double bam_aux2d(const uint8_t *s);
char bam_aux2A(const uint8_t *s);
char *bam_aux2Z(const uint8_t *s);

int bam_aux_del(bam1_t *b, uint8_t *s);
void bam_aux_append(bam1_t *b, const char tag[2], char type, int len, uint8_t *data);
uint8_t *bam_aux_get_core(bam1_t *b, const char tag[2]);





uint32_t bam_calend(const bam1_core_t *c, const uint32_t *cigar);


int32_t bam_cigar2qlen(const bam1_core_t *c, const uint32_t *cigar);

#ifdef __cplusplus
}
#endif


static inline int bam_reg2bin(uint32_t beg, uint32_t end)
{
  --end;

  if (beg >> 14 == end >> 14) return 4681 + (beg >> 14);

  if (beg >> 17 == end >> 17) return 585 + (beg >> 17);

  if (beg >> 20 == end >> 20) return 73 + (beg >> 20);

  if (beg >> 23 == end >> 23) return 9 + (beg >> 23);

  if (beg >> 26 == end >> 26) return 1 + (beg >> 26);

  return 0;
}


static inline bam1_t *bam_copy1(bam1_t *bdst, const bam1_t *bsrc)
{
  uint8_t *data = bdst->data;
  int m_data = bdst->m_data;

  if (m_data < bsrc->data_len)
    {
      m_data = bsrc->data_len;
      kroundup32(m_data);
      data = (uint8_t *) realloc(data, m_data);
    }

  memcpy(data, bsrc->data, bsrc->data_len);
  *bdst = *bsrc;

  bdst->m_data = m_data;
  bdst->data = data;
  return bdst;
}


static inline bam1_t *bam_dup1(const bam1_t *src)
{
  bam1_t *b;
  b = bam_init1();
  *b = *src;
  b->m_data = b->data_len;
  b->data = (uint8_t *) calloc(b->data_len, 1);
  memcpy(b->data, src->data, b->data_len);
  return b;
}

static inline int bam_aux_type2size(int x)
{
  if (x == 'C' || x == 'c' || x == 'A') return 1;
  else if (x == 'S' || x == 's') return 2;
  else if (x == 'I' || x == 'i' || x == 'f') return 4;
  else return 0;
}


#endif
