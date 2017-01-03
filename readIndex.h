#ifndef _READ_INDEX_H_
#define _READ_INDEX_H_

#define MAX_SEQ_NAME_LENGTH 256

#pragma pack(2)
typedef struct Translate
{
  unsigned int startPos;
  unsigned int correction;
  unsigned short chrID;
} Translate;
#pragma pack()

#pragma pack(2)
typedef struct Annotation
{
  int gi;
  char text[MAX_SEQ_NAME_LENGTH + 1];
  char decoratedText[MAX_SEQ_NAME_LENGTH + 1];
} Annotation;
#pragma pack()

#pragma pack(4)
typedef struct SeqOffset
{
  unsigned int startPos;
  unsigned int endPos;
  int firstAmbiguityIndex;
  int lastAmbiguityIndex;
} SeqOffset;
#pragma pack()

size_t loadTranslateWithTranslateSize(const char *inputFileName, unsigned int &dnaLength,
                                      unsigned int **ambiguityMap, Translate **translate);

void loadTranslate(const char *inputFileName, unsigned int &dnaLength,
                   unsigned int **ambiguityMap, Translate **translate);

void freeTranslate(unsigned int *ambiguityMap, Translate *translate);

void loadSeqInfo(const char *inputFileName, unsigned int &dnaLength,
                 Annotation **annotation, SeqOffset **seqOffset, unsigned int &numOfSeq);

void freeSeqInfo(Annotation *annotation, SeqOffset *seqOffset);

#endif
