#ifndef _SIMPLE_MEMORY_POOL_H_
#define _SIMPLE_MEMORY_POOL_H_

#include <stdlib.h>

typedef struct MemoryPool
{
  unsigned int curPtr;
  size_t capacity;
  void *address;
} MemoryPool;

MemoryPool *createPool(size_t size);

void destroyPool(MemoryPool *pool);

unsigned int pmallocIndex(MemoryPool *pool, size_t size);

unsigned int premallocIndex(MemoryPool *pool, void *address, size_t size);

void *getAddress(MemoryPool *pool, unsigned int index);

#endif
