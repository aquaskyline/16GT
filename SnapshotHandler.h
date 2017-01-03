#ifndef _SNAPSHOT_HANDLER_H_
#define _SNAPSHOT_HANDLER_H_

#include "struct.h"
#include "coreStruct.h"

void writeSnpInfoSnapshot(SnpBundle *snpStatBundle,
                          ExomeRegion *region, unsigned int regionSize,
                          unsigned char isExome,
                          const char *snapshotFilename);

void readSnpInfoSnapshot(SnpBundle *snpStatBundle, unsigned int &textLength,
                         const char *snapshotFilename);


#endif
