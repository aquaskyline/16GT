#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "readIndex.h"
#include "indexFunction.h"
#include "struct.h"

#define BUFFER_SIZE 1024
#define INIT_REGION_SIZE (1048576)

#define FT_GFF 1
#define FT_BED 2

char isOverlap(unsigned int a0, unsigned int a1, unsigned int b0, unsigned int b1,
               unsigned int &c0, unsigned int &c1) {
  char flag = 0;

  unsigned int l0, l1, r0, r1;
  if (a0 < b0) {
    l0 = a0;
    l1 = a1;
    r0 = b0;
    r1 = b1;
  } else {
    l0 = b0;
    l1 = b1;
    r0 = a0;
    r1 = a1;
  }

  if (r0 < l1) {
    flag = 1;
    c0 = l0;
    if (r1 < l1) {
      c1 = l1;
    } else {
      c1 = r1;
    }
  }

  return flag;
}

void parseInput(const char *dataFileName, char fileType,
                unsigned int *ambMap, Translate *translate, unsigned int textLength,
                Annotation *annotation, unsigned int numOfSeq,
                const char *geneListFileName) {
#define checkBuffer() \
if ( bufferPtr >= readNum ) { readNum = fread ( buffer, sizeof ( char ), BUFFER_SIZE, dataFile ); bufferPtr = 0; }

#define EXTEND_ARRAY(_size, _cap, _ext, _array, _type) \
do { \
    if ( ( _size ) == ( _cap ) ) { ( _cap ) += ( _ext ); _array = ( _type * ) realloc ( _array, ( _cap ) * sizeof ( _type ) ); } \
} while ( 0 )

  FILE *dataFile = fopen(dataFileName, "r");
  if (!dataFile) {
    fprintf(stderr, "Cannot open file : %s\n", dataFileName);
    exit(1);
  }
  FILE *geneListFile = fopen(geneListFileName, "wb");
  if (!geneListFile) {
    fprintf(stderr, "Cannot open file : %s\n", geneListFileName);
    exit(1);
  }

  unsigned int readNum = 0;
  char buffer[BUFFER_SIZE] = {'\0'};
  unsigned int bufferPtr = 0;
  readNum = fread(buffer, sizeof(char), BUFFER_SIZE, dataFile);

  char chrName[MAX_SEQ_NAME_LENGTH];
  unsigned int chrPos0, chrPos1;
  unsigned int chrID;
  unsigned int ambPos0, ambPos1;

  unsigned int maxRegions = INIT_REGION_SIZE;
  unsigned int numRegions = 0;
  IndexRegion *indexRegion = (IndexRegion *) malloc(maxRegions * sizeof(IndexRegion));
  printf("\e[?25l");
  while (1) {
    
    checkBuffer();
    if (readNum <= 0) {
      break;
    }
    while (buffer[bufferPtr] == '#') {
      do {
        checkBuffer();
      } while (buffer[bufferPtr++] != '\n');
    }

    checkBuffer();

    {
      if (numRegions % 10000 == 0) printf("Handled %u entries    \n", numRegions);
      unsigned int namePtr = 0;
      do {
        if (namePtr >= MAX_SEQ_NAME_LENGTH - 1) {
          bufferPtr++;
          checkBuffer();
          continue;
        }
        chrName[namePtr++] = buffer[bufferPtr++];
        checkBuffer();
      } while (buffer[bufferPtr] != ' ' && buffer[bufferPtr] != '\t');
      chrName[namePtr] = '\0'; 

      if (fileType == FT_GFF) {
        
        do {
          bufferPtr++;
          checkBuffer();
        } while (buffer[bufferPtr] != ' ' && buffer[bufferPtr] != '\t');

        
        do {
          bufferPtr++;
          checkBuffer();
        } while (buffer[bufferPtr] != ' ' && buffer[bufferPtr] != '\t');
      }
      bufferPtr++;

      chrPos0 = 0;
      checkBuffer();
      while (buffer[bufferPtr] != ' ' && buffer[bufferPtr] != '\t') {
        if (buffer[bufferPtr] < '0' || buffer[bufferPtr] > '9') {
          fprintf(stderr, "Error in parsing gff file : position is not a number\n");
          printf("\e[?25h");
          exit(1);
        }
        chrPos0 = chrPos0 * 10 + buffer[bufferPtr++] - '0';
        checkBuffer();
      }
      bufferPtr++;

      chrPos1 = 0;
      checkBuffer();
      while (buffer[bufferPtr] != ' ' && buffer[bufferPtr] != '\t' && buffer[bufferPtr] != '\n') {
        if (buffer[bufferPtr] < '0' || buffer[bufferPtr] > '9') {
          fprintf(stderr, "Error in parsing gff file : position is not a number\n");
          printf("\e[?25h");
          exit(1);
        }
        chrPos1 = chrPos1 * 10 + buffer[bufferPtr++] - '0';
        checkBuffer();
      }

      if (readNum <= 0) {
        
        printf("File Format Error\n");
        break;
      }

      if (buffer[bufferPtr] != '\n') {
        bufferPtr++;
        
        do {
          checkBuffer();
        } while (buffer[bufferPtr++] != '\n' && readNum);
        checkBuffer();
      } else {
        bufferPtr++;
      }

      chrID = getChrIDFromName(annotation, numOfSeq, chrName);
      if (chrID <= numOfSeq) {
        ambPos0 = getAmbPos(chrID, chrPos0, ambMap, translate, textLength);
        ambPos1 = getAmbPos(chrID, chrPos1, ambMap, translate, textLength);
        if (ambPos0 < textLength && ambPos1 < textLength) {
          EXTEND_ARRAY (numRegions, maxRegions, maxRegions, indexRegion, IndexRegion);
          indexRegion[numRegions].startPos = ambPos0;
          indexRegion[numRegions].endPos = ambPos1;
          numRegions++;
        }
      }
    }
  }
  if (numRegions % 10000 == 0) printf("Handled %u entries    \n", numRegions);


#define NUM_BUCKETS 65536
#define SORT_DIGIT(array, var, shift_bits, num_elements, buffer, record_type) \
    do { unsigned int buckets[NUM_BUCKETS]; memset(buckets, 0, NUM_BUCKETS * sizeof(unsigned int)); \
    for (unsigned int i=0; i<num_elements; ++i) buckets[((array)[i].var>>shift_bits)&0xFFFF]++; \
    for (unsigned int b=0, acc=0; b<NUM_BUCKETS; ++b) {unsigned int t=acc; acc+=buckets[b]; buckets[b]=t;} \
    for (unsigned int i=0; i<num_elements; ++i) buffer[buckets[((array)[i].var>>shift_bits)&0xFFFF]++] = array[i]; \
    { record_type *t = buffer; buffer = array; array = t; } } while (0)

  IndexRegion *sortBuffer = (IndexRegion *) malloc(numRegions * sizeof(IndexRegion));

  SORT_DIGIT (indexRegion, startPos, 0, numRegions, sortBuffer, IndexRegion);
  SORT_DIGIT (indexRegion, startPos, 16, numRegions, sortBuffer, IndexRegion);

  free(sortBuffer);

  unsigned int newSize = 0;
  unsigned int i = 0;
  unsigned int nextI = i + 1;
  unsigned int _startPos, _endPos;
  while (1) {
    if (i == numRegions - 1 || nextI == numRegions) {
      newSize++;
      break;
    }
    if (isOverlap(indexRegion[i].startPos, indexRegion[i].endPos,
                  indexRegion[nextI].startPos, indexRegion[nextI].endPos,
                  _startPos, _endPos)) {
      indexRegion[newSize].startPos = _startPos;
      indexRegion[newSize].endPos = _endPos;
      nextI++;
    } else {
      newSize++;
      i = nextI++;
      indexRegion[newSize].startPos = indexRegion[i].startPos;
      indexRegion[newSize].endPos = indexRegion[i].endPos;
    }
  }

  fwrite(&newSize, sizeof(unsigned int), 1, geneListFile);
  fwrite(indexRegion, sizeof(IndexRegion), newSize, geneListFile);

  printf("\n");
  printf("Removed Overlapped entries : %u remained\n", newSize);
  printf("\e[?25h");

  fclose(dataFile);
  fclose(geneListFile);

  free(indexRegion);
}

int main(int argc, char **argv) {
  if (argc != 5) {
    fprintf(stderr, "Usage:\n\t%s [Index] [Input Region List] [Output List Name] [format]\n", argv[0]);
    fprintf(stderr, "\tformat: -bed OR -gff\n");
    return 1;
  }

  char *translateFileName;
  char *annotationFileName;

  unsigned int indexLength = strlen(argv[1]);
  translateFileName = (char *) malloc(indexLength + 5);
  annotationFileName = (char *) malloc(indexLength + 5);
  sprintf(translateFileName, "%s.tra", argv[1]);
  sprintf(annotationFileName, "%s.ann", argv[1]);

  printf("Loading Index\n");
  unsigned int dnaLength = 0;
  unsigned int *ambiguityMap;
  Translate *translate;
  Annotation *annotation;
  unsigned int numOfSeq = 0;

  loadTranslate(translateFileName, dnaLength, &ambiguityMap, &translate);
  printf("Loaded Translate : %u\n", dnaLength);
  loadSeqInfo(annotationFileName, dnaLength, &annotation, NULL, numOfSeq);
  printf("Loaded Annotation : %u, #seqeunce : %u\n", dnaLength, numOfSeq);

  char *geneListFileName = argv[3];

  char fileType = 0;
  if (strcmp(argv[4], "-gff") == 0) {
    fileType = FT_GFF;
  } else if (strcmp(argv[4], "-bed") == 0) {
    fileType = FT_BED;
  }

  if (fileType == 0) {
    printf("Invalid input format\n");
    exit(1);
  }

  parseInput(argv[2], fileType, ambiguityMap, translate, dnaLength, annotation, numOfSeq, geneListFileName);

  freeTranslate(ambiguityMap, translate);
  freeSeqInfo(annotation, NULL);

  free(translateFileName);
  free(annotationFileName);

  return 0;
}
