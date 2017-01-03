

#ifndef __BGZF_H
#define __BGZF_H

#include <stdint.h>
#include <stdio.h>
#include <zlib.h>

#ifdef _USE_KNETFILE
#include "knetfile.h"
#endif



typedef struct
{
  int file_descriptor;
  char open_mode;
  int16_t owned_file, compress_level;
#ifdef _USE_KNETFILE
  union
  {
    knetFile *fpr;
    FILE *fpw;
  } x;
#else
  FILE *file;
#endif
  int uncompressed_block_size;
  int compressed_block_size;
  void *uncompressed_block;
  void *compressed_block;
  int64_t block_address;
  int block_length;
  int block_offset;
  int cache_size;
  const char *error;
  void *cache;
} BGZF;

#ifdef __cplusplus
extern "C" {
#endif


BGZF *bgzf_fdopen(int fd, const char *__restrict mode);


BGZF *bgzf_open(const char *path, const char *__restrict mode);


int bgzf_close(BGZF *fp);


int bgzf_read(BGZF *fp, void *data, int length);


int bgzf_write(BGZF *fp, const void *data, int length);


#define bgzf_tell(fp) ((fp->block_address << 16) | (fp->block_offset & 0xFFFF))


int64_t bgzf_seek(BGZF *fp, int64_t pos, int where);


void bgzf_set_cache_size(BGZF *fp, int cache_size);

int bgzf_check_EOF(BGZF *fp);
int bgzf_read_block(BGZF *fp);
int bgzf_flush(BGZF *fp);
int bgzf_flush_try(BGZF *fp, int size);
int bgzf_check_bgzf(const char *fn);

#ifdef __cplusplus
}
#endif

static inline int bgzf_getc(BGZF *fp)
{
  int c;

  if (fp->block_offset >= fp->block_length)
    {
      if (bgzf_read_block(fp) != 0) return -2;

      if (fp->block_length == 0) return -1;
    }

  c = ((unsigned char *) fp->uncompressed_block)[fp->block_offset++];

  if (fp->block_offset == fp->block_length)
    {
#ifdef _USE_KNETFILE
      fp->block_address = knet_tell(fp->x.fpr);
#else
      fp->block_address = ftello(fp->file);
#endif
      fp->block_offset = 0;
      fp->block_length = 0;
    }

  return c;
}

#endif
