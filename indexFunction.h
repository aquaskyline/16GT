#ifndef _INDEX_FUNCTION_H_
#define _INDEX_FUNCTION_H_

#define GRID_SAMPLING_FACTOR_2_POWER 18
#define LOOP_BOUND 1000

#include "readIndex.h"


unsigned int getAmbPos(unsigned short chr_id, unsigned int offset,
                       unsigned int *ambiguityMap, Translate *translate,
                       unsigned int dnaLength);


unsigned int getChrIDFromName(Annotation *annotation, unsigned int numOfSeq, const char *chrName);

#endif
