#ifndef _INDEX_STRUCT_H_
#define _INDEX_STRUCT_H_

#pragma pack(4)
typedef struct IndexRegion
{
  unsigned int startPos;
  unsigned int endPos;
} IndexRegion;
#pragma pack()

typedef IndexRegion GeneRegion;
typedef IndexRegion ExomeRegion;

#pragma pack(2)
typedef struct IndelDataBase
{
  unsigned int startPos;
  unsigned int pattArrayIndex;
  unsigned short indelTypeLength;
} IndelDB;
#pragma pack()

#endif
