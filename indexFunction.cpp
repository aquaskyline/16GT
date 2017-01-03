#include <stdio.h>
#include <string.h>

#include "indexFunction.h"

unsigned int getAmbPos(unsigned short chr_id, unsigned int offset,
                       unsigned int *ambiguityMap, Translate *translate,
                       unsigned int dnaLength) {
  long long low = 0;
  long long high = dnaLength - 1;
  int loop_count = 0;
  while (low <= high) {
    loop_count++;
    if (loop_count == LOOP_BOUND) {
      break;
    }
    unsigned int ambPos = (low + high) / 2;
    unsigned int approxIndex, approxValue;
    unsigned int correctPosition = ambPos;
    approxIndex = ambPos >> GRID_SAMPLING_FACTOR_2_POWER;
    approxValue = ambiguityMap[approxIndex];
    while (translate[approxValue].startPos > ambPos) {
      approxValue--;
    }
    correctPosition -= translate[approxValue].correction;
    unsigned short chr = translate[approxValue].chrID;

    if (chr < chr_id || (chr == chr_id && correctPosition < offset)) {
      low = ambPos + 1;
    } else if (chr > chr_id || (chr == chr_id && correctPosition > offset)) {
      high = ambPos - 1;
    } else {
      return ambPos;
    }
  }

  return dnaLength + 2;
}


unsigned int getChrIDFromName(Annotation *annotation, unsigned int numOfSeq, const char *chrName) {
  unsigned int i;
  for (i = 0; i < numOfSeq; ++i) {
    if (strcmp(annotation[i].text, chrName) == 0) {
      return i + 1;
    }
  }
  if (strlen(chrName) >= 3 && chrName[0] == 'c' && chrName[1] == 'h' && chrName[2] == 'r') {
    for (i = 0; i < numOfSeq; ++i) {
      if (strcmp(annotation[i].text, chrName + 3) == 0) {
        return i + 1;
      }
    }
  } else {
    char tmp[1000];
    tmp[0] = 0;
    strcat(tmp, "chr");
    strcat(tmp, chrName);
    for (i = 0; i < numOfSeq; ++i) {
      if (strcmp(annotation[i].text, chrName) == 0) {
        return i + 1;
      }
    }
  }
  return i + 1;
}
