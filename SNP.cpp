#include <unistd.h>
#include <pthread.h>
#include <zlib.h>
#include <sys/mman.h>

#include "SNP.h"
#include "CPUfunctions.h"

#include <assert.h>




#define INIT_OVERFLOW_ARRAY_SIZE 1024*1024
#define EXTRA_OVERFLOW_ARRAY_SIZE 125829120
#define OVERFLOW_BATCH_SIZE 1048576u

#define OVERFLOW_BUFFER_SIZE 1048576u

#define COUNTER_NOT_OVERFLOW 0
#define COUNTER_OVERFLOW 1

#define UPDATE_INVALID_POS

#define CIGAR_MATCHMISMATCH 'M'
#define CIGAR_INSERT 'I'
#define CIGAR_DELETE 'D'
#define CIGAR_SOFT_CLIP 'S'

SnpOverflowBufferArray ** snpOverflowBufferArray;





void destroySnpCounter ( SnpCounter * snpCounter, unsigned int dnaLength )
{
  
  munmap ( snpCounter, dnaLength * sizeof ( SnpCounter ) );
}

void constructSnpOverflowCounter ( SnpOverflowCounterArray ** snpOverflowCounterArrayPtr, int numOfCPUThreads )
{
  
  *snpOverflowCounterArrayPtr = ( SnpOverflowCounterArray * ) xmalloc ( numOfCPUThreads * sizeof(SnpOverflowCounterArray) );

  unsigned int initSize = INIT_OVERFLOW_ARRAY_SIZE / numOfCPUThreads;

  int i;
  for ( i = 0; i < numOfCPUThreads; ++i )
  {
    
    ( * snpOverflowCounterArrayPtr )[i].counters = ( SnpOverflowCounter * ) xmalloc ( initSize * sizeof(SnpOverflowCounter) );
    ( * snpOverflowCounterArrayPtr )[i].limit = initSize; 
    ( * snpOverflowCounterArrayPtr )[i].size = 0;
    memset ( ( * snpOverflowCounterArrayPtr )[i].counters, 0, initSize * sizeof(SnpOverflowCounter) );
  }
}


unsigned int addSnpOverflowCounterToArray ( SnpOverflowCounterArray & snpOverflowCounterArray, int numOfCPUThreads )
{
  if ( snpOverflowCounterArray.size == snpOverflowCounterArray.limit )
  {
    SnpOverflowCounter * oldCounters = snpOverflowCounterArray.counters;

    snpOverflowCounterArray.limit = snpOverflowCounterArray.limit + EXTRA_OVERFLOW_ARRAY_SIZE / numOfCPUThreads; 
    snpOverflowCounterArray.counters = ( SnpOverflowCounter * ) xmalloc ( snpOverflowCounterArray.limit * sizeof(SnpOverflowCounter) );
    memset ( snpOverflowCounterArray.counters, 0, snpOverflowCounterArray.limit * sizeof(SnpOverflowCounter) );

    memcpy ( snpOverflowCounterArray.counters, oldCounters, snpOverflowCounterArray.size * sizeof(SnpOverflowCounter) );

    xfree ( oldCounters );
  }

  return snpOverflowCounterArray.size++;
}

void destroySnpOverflowCounter ( SnpOverflowCounterArray * snpOverflowCounterArray, int numOfCPUThreads )
{
  int i;
  for ( i = 0; i < numOfCPUThreads; ++i )
  {
    xfree ( snpOverflowCounterArray[i].counters );
    
  }

  xfree ( snpOverflowCounterArray );
}





void readFileForDirectionalSoftclip ( FILE * scFile, SnpDirectionalSoftclip * softclips, unsigned int &size )
{
  fseek( scFile, 0L, SEEK_END );
  int numberOfElement = ftell( scFile );
  fseek( scFile, 0L, SEEK_SET);
  fread( (softclips + size), sizeof ( SnpDirectionalSoftclip ), numberOfElement / sizeof ( SnpDirectionalSoftclip ), scFile );
#ifdef DEBUG_DIRECTIONAL_SOFTCLIP
  fprintf ( stderr, "Size of softclip entry file: %u\n", numberOfElement );
  fprintf ( stderr, "Size of SnpDirectionalSoftclip: %u\n", sizeof ( SnpDirectionalSoftclip ) );
  fprintf ( stderr, "Number of softclip entries: %u\n", numberOfElement/ sizeof ( SnpDirectionalSoftclip ) );
#endif
  size += numberOfElement / sizeof ( SnpDirectionalSoftclip );    
}

void constructDirectionalSoftclipCounterFromFiles ( int numOfCPUThreads, 
    FILE ** snpDirectionalSoftclipFilePtr,
    FILE * &snpDirectionalSoftclipFileDpPtr,
    FILE * &snpDirectionalSoftclipFileUnpairPtr,
    FILE ** snpDirectionalSoftclipDedupFilePtr,
    FILE * &snpDirectionalSoftclipDedupFileDpPtr,
    FILE * &snpDirectionalSoftclipDedupFileUnpairPtr,
    SnpDirectionalSoftclipCounter * &snpDirectionalSoftclipCounter,
    unsigned int & snpDirectionalSoftclipCounterSize,
    unsigned int & snpDirectionalSoftclipCounterCapacity )
{
  unsigned int arraySize = 0, dedupArraySize = 0;
  unsigned int arrayCapacity = 0, dedupArrayCapacity = 0;
  for ( int i=0; i<numOfCPUThreads;++i )
  {
    fseek(snpDirectionalSoftclipFilePtr[i], 0L, SEEK_END);
    arrayCapacity += ftell(snpDirectionalSoftclipFilePtr[i]);
    fseek(snpDirectionalSoftclipFilePtr[i], 0L, SEEK_SET);
  }
  fseek(snpDirectionalSoftclipFileDpPtr, 0L, SEEK_END);
  arrayCapacity += ftell(snpDirectionalSoftclipFileDpPtr);
  fseek(snpDirectionalSoftclipFileDpPtr, 0L, SEEK_SET);
  fseek(snpDirectionalSoftclipFileUnpairPtr, 0L, SEEK_END);
  arrayCapacity += ftell(snpDirectionalSoftclipFileUnpairPtr);
  fseek(snpDirectionalSoftclipFileUnpairPtr, 0L, SEEK_SET);
  arrayCapacity = arrayCapacity / sizeof ( SnpDirectionalSoftclip );

  for ( int i=0; i<numOfCPUThreads;++i )
  {
    fseek(snpDirectionalSoftclipDedupFilePtr[i], 0L, SEEK_END);
    dedupArrayCapacity += ftell(snpDirectionalSoftclipDedupFilePtr[i]);
    fseek(snpDirectionalSoftclipDedupFilePtr[i], 0L, SEEK_SET);
  }
  fseek(snpDirectionalSoftclipDedupFileDpPtr, 0L, SEEK_END);
  dedupArrayCapacity += ftell(snpDirectionalSoftclipDedupFileDpPtr);
  fseek(snpDirectionalSoftclipDedupFileDpPtr, 0L, SEEK_SET);
  fseek(snpDirectionalSoftclipDedupFileUnpairPtr, 0L, SEEK_END);
  dedupArrayCapacity += ftell(snpDirectionalSoftclipDedupFileUnpairPtr);
  fseek(snpDirectionalSoftclipDedupFileUnpairPtr, 0L, SEEK_SET);
  dedupArrayCapacity = dedupArrayCapacity / sizeof ( SnpDirectionalSoftclip );

  SnpDirectionalSoftclip * softclips = ( SnpDirectionalSoftclip * ) xmalloc ( arrayCapacity * sizeof ( SnpDirectionalSoftclip ) );
  SnpDirectionalSoftclip * softclipsBuffer = ( SnpDirectionalSoftclip * ) xmalloc ( arrayCapacity * sizeof ( SnpDirectionalSoftclip ) );
  SnpDirectionalSoftclip * dedupSoftclips = ( SnpDirectionalSoftclip * ) xmalloc ( dedupArrayCapacity * sizeof ( SnpDirectionalSoftclip ) );
  SnpDirectionalSoftclip * dedupSoftclipsBuffer = ( SnpDirectionalSoftclip * ) xmalloc ( dedupArrayCapacity * sizeof ( SnpDirectionalSoftclip ) );

  for ( int i=0; i<numOfCPUThreads;++i )
  {
    readFileForDirectionalSoftclip( snpDirectionalSoftclipFilePtr[i], softclips, arraySize );
  }
  readFileForDirectionalSoftclip( snpDirectionalSoftclipFileDpPtr, softclips, arraySize );
  readFileForDirectionalSoftclip( snpDirectionalSoftclipFileUnpairPtr, softclips, arraySize );
  assert ( arraySize == arrayCapacity );
  unsigned int sortingBuckets[NUM_BUCKETS];
  SORT_DIGIT_BY_VAR ( softclips, ambPosition, 0, arraySize, sortingBuckets, softclipsBuffer, SnpDirectionalSoftclip );
  SORT_DIGIT_BY_VAR ( softclips, ambPosition, 16, arraySize, sortingBuckets, softclipsBuffer, SnpDirectionalSoftclip );

  for ( int i=0; i<numOfCPUThreads;++i )
  {
    readFileForDirectionalSoftclip( snpDirectionalSoftclipDedupFilePtr[i], dedupSoftclips, dedupArraySize );
  }
  readFileForDirectionalSoftclip ( snpDirectionalSoftclipDedupFileDpPtr, dedupSoftclips, dedupArraySize );
  readFileForDirectionalSoftclip( snpDirectionalSoftclipDedupFileUnpairPtr, dedupSoftclips, dedupArraySize );
  assert ( dedupArraySize == dedupArrayCapacity );   
  SORT_DIGIT_BY_VAR ( dedupSoftclips, ambPosition, 0, dedupArraySize, sortingBuckets, dedupSoftclipsBuffer, SnpDirectionalSoftclip );
  SORT_DIGIT_BY_VAR ( dedupSoftclips, ambPosition, 16, dedupArraySize, sortingBuckets, dedupSoftclipsBuffer, SnpDirectionalSoftclip );

  
  int posCount = ( arraySize != 0 ? 1 : 0 );
  for ( arraySize = 1; arraySize < arrayCapacity; ++arraySize )
  {
    if ( softclips[arraySize].ambPosition != softclips[arraySize-1].ambPosition )
    { ++posCount; }
  }

  snpDirectionalSoftclipCounter = ( SnpDirectionalSoftclipCounter * ) xmalloc ( posCount * sizeof ( SnpDirectionalSoftclipCounter ) );
  memset ( snpDirectionalSoftclipCounter, 0, posCount * sizeof ( SnpDirectionalSoftclipCounter ) );

  unsigned int idx = 0;
  arraySize = 0, dedupArraySize = 0;
  while ( arraySize < arrayCapacity )
  {
    SnpDirectionalSoftclipCounter tmp;
    tmp.ambPosition = softclips[arraySize].ambPosition;
    tmp.leftSoftclipCount = 0;
    tmp.rightSoftclipCount = 0;
    while ( arraySize < arrayCapacity && softclips[arraySize].ambPosition == tmp.ambPosition )
    {
      switch ( softclips[arraySize].leftOrRight )
      {
        case 0: ++tmp.leftSoftclipCount; break;
        case 1: ++tmp.rightSoftclipCount; break;
      }
      ++arraySize;
    }
    while ( dedupArraySize < dedupArrayCapacity && dedupSoftclips[dedupArraySize].ambPosition == tmp.ambPosition )
    {
      switch ( dedupSoftclips[dedupArraySize].leftOrRight )
      {
        case 0: --tmp.leftSoftclipCount; break;
        case 1: --tmp.rightSoftclipCount; break;
      }
      ++dedupArraySize;            
    }
    if ( tmp.leftSoftclipCount > 0 || tmp.rightSoftclipCount > 0 )
    { 
      snpDirectionalSoftclipCounter[idx] = tmp;
      ++idx; 
    }
  }
  assert ( dedupArraySize == dedupArrayCapacity );
  snpDirectionalSoftclipCounterSize = idx;
  snpDirectionalSoftclipCounterCapacity = posCount;
  xfree ( softclipsBuffer );
  xfree ( dedupSoftclipsBuffer );
  xfree ( softclips );
  xfree ( dedupSoftclips );
}

void constructDirectionalSoftclipFileNames( char * outputPrefix, int numOfCPUThreads,
    char ** snpDirectionalSoftclipFileName, char * & snpDirectionalSoftclipDpFileName, char * &snpDirectionalSoftclipUnpairFileName,
    char ** snpDirectionalSoftclipDedupFileName, char * &snpDirectionalSoftclipDedupDpFileName , char * &snpDirectionalSoftclipDedupUnpairFileName )
{
  
  for ( int i = 0; i < numOfCPUThreads; i++ )
  {
    snpDirectionalSoftclipFileName[i] = ( char * ) xmalloc ( strlen ( outputPrefix ) + 7 );
  }
  snpDirectionalSoftclipDpFileName = ( char * ) xmalloc ( strlen ( outputPrefix ) + 7 );
  snpDirectionalSoftclipUnpairFileName = ( char * ) xmalloc ( strlen ( outputPrefix ) + 11 );

  for ( int i = 0; i < numOfCPUThreads; i++ )
  {
    sprintf ( snpDirectionalSoftclipFileName[i], "%s.sc.%d", outputPrefix, i );
  }
  sprintf( snpDirectionalSoftclipDpFileName, "%s.sc.dp", outputPrefix);
  sprintf ( snpDirectionalSoftclipUnpairFileName, "%s.sc.unpair", outputPrefix );

  
  for ( int i = 0; i < numOfCPUThreads; i++ )
  {
    snpDirectionalSoftclipDedupFileName[i] = ( char * ) xmalloc ( strlen ( outputPrefix ) + 10 );
  }
  snpDirectionalSoftclipDedupDpFileName = ( char * ) xmalloc ( strlen ( outputPrefix ) + 10 );
  snpDirectionalSoftclipDedupUnpairFileName = ( char * ) xmalloc ( strlen ( outputPrefix ) + 14 );

  for ( int i = 0; i < numOfCPUThreads; i++ )
  {
    sprintf ( snpDirectionalSoftclipDedupFileName[i], "%s.rm_sc.%d", outputPrefix, i );
  }
  sprintf ( snpDirectionalSoftclipDedupDpFileName, "%s.rm_sc.dp", outputPrefix );
  sprintf ( snpDirectionalSoftclipDedupUnpairFileName, "%s.rm_sc.unpair", outputPrefix );

}

unsigned int getDirectionalSoftclipSumByPos ( SnpDirectionalSoftclipCounter * dscCounter,
    unsigned int dscSize,
    unsigned int pos )
{
  unsigned int low = 0;
  unsigned int high = dscSize;
  unsigned int mid;
  do
  {
    mid = ( low + high ) / 2;
    if ( dscCounter[mid].ambPosition > pos )
    { high = mid; }
    else if ( dscCounter[mid].ambPosition < pos )
    { low = mid + 1; }
    else
    { return dscCounter[mid].leftSoftclipCount + dscCounter[mid].rightSoftclipCount; }
  } while ( low < high );
  return 0;
}

void destroyMemoryPool ( MemoryPool * pool )
{
  destroyPool ( pool );
}


void constructSnpOverflowBufferArray ( int numOfCPUThreads )
{
  snpOverflowBufferArray = ( SnpOverflowBufferArray ** ) xmalloc ( numOfCPUThreads * sizeof(SnpOverflowBufferArray *) );

  int i;
  for ( i = 0; i < numOfCPUThreads; ++i )
  {
    snpOverflowBufferArray[i] = ( SnpOverflowBufferArray * ) malloc ( sizeof(SnpOverflowBufferArray) );
    snpOverflowBufferArray[i]->buffer = ( SnpOverflowBuffer * ) malloc ( OVERFLOW_BATCH_SIZE * sizeof(SnpOverflowBuffer) );
    snpOverflowBufferArray[i]->limit = OVERFLOW_BATCH_SIZE; 
    snpOverflowBufferArray[i]->size = 0;
    memset ( snpOverflowBufferArray[i]->buffer, 0, OVERFLOW_BATCH_SIZE * sizeof(SnpOverflowBuffer) );
  }
}


unsigned int addSnpOverflowBufferToArray ( SnpOverflowBufferArray * bufferArray,
    unsigned int position, char strand, unsigned char info )
{
  if ( bufferArray->size == bufferArray->limit )
  {
    SnpOverflowBuffer * oldBuffer = bufferArray->buffer;

    bufferArray->limit = bufferArray->limit * 2;
    bufferArray->buffer = ( SnpOverflowBuffer * ) malloc ( bufferArray->limit * sizeof(SnpOverflowBuffer) );

    memcpy ( bufferArray->buffer, oldBuffer, bufferArray->size * sizeof(SnpOverflowBuffer) );

    free ( oldBuffer );
  }

  bufferArray->buffer[bufferArray->size].position = position;
  bufferArray->buffer[bufferArray->size].strand = strand;
  bufferArray->buffer[bufferArray->size].info = info;
  return bufferArray->size++;
}

unsigned int addSnpInsertionBufferToArray ( SnpOverflowBufferArray * bufferArray,
    unsigned int position, char strand, char avgWeight,
    unsigned char insertLength, unsigned char * insertSeq )
{
  if ( bufferArray->size == bufferArray->limit )
  {
    SnpOverflowBuffer * oldBuffer = bufferArray->buffer;

    bufferArray->limit = bufferArray->limit * 2;
    bufferArray->buffer = ( SnpOverflowBuffer * ) malloc ( bufferArray->limit * sizeof(SnpOverflowBuffer) );

    memcpy ( bufferArray->buffer, oldBuffer, bufferArray->size * sizeof(SnpOverflowBuffer) );

    free ( oldBuffer );
  }

  bufferArray->buffer[bufferArray->size].position = position;
  bufferArray->buffer[bufferArray->size].strand = strand;
  bufferArray->buffer[bufferArray->size].info = 0x80;
  bufferArray->buffer[bufferArray->size].weight = avgWeight;
  bufferArray->buffer[bufferArray->size].length = insertLength;
  memcpy ( bufferArray->buffer[bufferArray->size].insertSeq, insertSeq, ( insertLength + 3 ) / 4 );
  return bufferArray->size++;
}

unsigned int addSnpDeletionBufferToArray ( SnpOverflowBufferArray * bufferArray,
    unsigned int position, char strand, char weight, unsigned char deleteLength )
{
  if ( bufferArray->size == bufferArray->limit )
  {
    SnpOverflowBuffer * oldBuffer = bufferArray->buffer;

    bufferArray->limit = bufferArray->limit * 2;
    bufferArray->buffer = ( SnpOverflowBuffer * ) malloc ( bufferArray->limit * sizeof(SnpOverflowBuffer) );

    memcpy ( bufferArray->buffer, oldBuffer, bufferArray->size * sizeof(SnpOverflowBuffer) );

    free ( oldBuffer );
  }

  bufferArray->buffer[bufferArray->size].position = position;
  bufferArray->buffer[bufferArray->size].strand = strand;
  bufferArray->buffer[bufferArray->size].info = 0xC0;
  bufferArray->buffer[bufferArray->size].weight = weight;
  bufferArray->buffer[bufferArray->size].length = deleteLength;
  return bufferArray->size++;
}


void emptySnpOverflowBuffer ( int numOfCPUThreads )
{
  int i;
  for ( i = 0; i < numOfCPUThreads; ++i )
  {
    free ( snpOverflowBufferArray[i]->buffer );
    snpOverflowBufferArray[i]->buffer = ( SnpOverflowBuffer * ) malloc ( OVERFLOW_BATCH_SIZE * sizeof(SnpOverflowBuffer) );
    snpOverflowBufferArray[i]->limit = OVERFLOW_BATCH_SIZE; 
    snpOverflowBufferArray[i]->size = 0;
    memset ( snpOverflowBufferArray[i]->buffer, 0, OVERFLOW_BATCH_SIZE * sizeof(SnpOverflowBuffer) );
  }
}

void addSnpOverflowBufferFromInfo ( unsigned int dnaLength, int numOfCPUThreads, unsigned int position, char strand, unsigned char info )
{
  int i = findRegionByPosition ( position, dnaLength, numOfCPUThreads );
  addSnpOverflowBufferToArray ( snpOverflowBufferArray[i], position, strand, info );
}

void addSnpOverflowBufferForInsertionFromInfo ( unsigned int dnaLength, int numOfCPUThreads, unsigned int position, char strand,
    char avgWeight, unsigned char insertLength, unsigned char * insertedSeq )
{
  int i = findRegionByPosition ( position, dnaLength, numOfCPUThreads );

  addSnpInsertionBufferToArray ( snpOverflowBufferArray[i], position, strand, avgWeight, insertLength, insertedSeq );
}

void addSnpOverflowBufferForDeletionFromInfo ( unsigned int dnaLength, int numOfCPUThreads, unsigned int position, char strand,
    char weight, unsigned char deleteLength )
{
  int i = findRegionByPosition ( position, dnaLength, numOfCPUThreads );
  addSnpDeletionBufferToArray ( snpOverflowBufferArray[i], position, strand, weight, deleteLength );
}

void destroySnpOverflowBufferArray ( int numOfCPUThreads )
{
  int i;
  for ( i = 0; i < numOfCPUThreads; ++i )
  {
    free ( snpOverflowBufferArray[i]->buffer );
    free ( snpOverflowBufferArray[i] );
  }

  xfree ( snpOverflowBufferArray );
}

















void updateSnpCounterForReads ( SnpBundle snpBundle, VcSetting * vcSetting, int trimHeadSize, int trimTailSize,
    unsigned char * recalScores,
    LongSoftClipPositionArray * longSoftClipArray, unsigned int longSoftClipThreshold,
    unsigned int position, const char strand,
    unsigned char * query, char * qualities, unsigned int readLength, char * cigar,
    FILE * overflowFilePtr, FILE * softclipFilePtr )
{
  

  unsigned int readPos = 0; 
  unsigned int genomePosition = position; 
  char overflowFlag = 0;
  
  char writeFileFlag = 0;
  unsigned char buffer[MAX_READ_LENGTH * 2] = { 0 };
  unsigned int bufferPos = 0;
  if ( cigar == NULL )  
  {
    while ( readPos < readLength )
    {
      if ( readPos < trimHeadSize )
      {
        ++genomePosition;
        ++readPos;
        continue;
      }

      if ( readPos >= readLength - trimTailSize )
      {
        break;
      }

      char baseBits = strand ? query[readPos] : 3 - query[readLength - readPos - 1];
      char qScore = strand ? qualities[readPos] : qualities[readLength - readPos - 1];
      if ( qScore < 0 || qScore > 40 )
      {
        if ( vcSetting->enableQualityCorrection )
        {
          qScore = qScore < 0 ? 0 : 40;
        }
        else
        {
          fprintf ( stderr, "Error in Base Quality. Please check your input read file.\n" );
          exit ( 1 );
        }
      }
      char recalScore = recalScores[qScore];
      
      char weight = vcSetting->weightMap[recalScore];
      if ( overflowFlag )
      {
        
        buffer[bufferPos++] = ( weight + ( baseBits << 3 ) + 0x40 );
      }
      else if ( updateSnpCounter ( snpBundle, vcSetting,
            genomePosition, baseBits, weight, strand, 1 ) == COUNTER_OVERFLOW )
      {
        
        overflowFlag = 1;
        buffer[bufferPos++] |= genomePosition & 255;
        buffer[bufferPos++] |= ( genomePosition >> 8 ) & 255;
        buffer[bufferPos++] |= ( genomePosition >> 16 ) & 255;
        buffer[bufferPos++] |= ( genomePosition >> 24 ) & 255;
        buffer[bufferPos++] = strand;

        
        buffer[bufferPos++] = ( weight + ( baseBits << 3 ) + 0x40 );  
      }

      ++genomePosition;
      ++readPos;
    }

    if ( overflowFlag )
    {
      
      buffer[bufferPos - 1] |= 0x20;
      fwrite ( buffer, 1, bufferPos, overflowFilePtr );
    }
  }
  else
  {
    unsigned int i = 0;
    int trimmedHead = 0;
    
    
    

    unsigned int lastMPos = 0;

    while ( readPos < readLength )
    {
      int range = 0;
      while ( cigar[i] <= 57 )
      {
        range = range * 10 + ( cigar[i] - 48 );
        ++i;
      }

      if ( cigar[i] == CIGAR_MATCHMISMATCH )
      {
        if ( readPos < trimHeadSize )
        {
          int trimRange = ( range < trimHeadSize - trimmedHead ? range : trimHeadSize - trimmedHead );
          genomePosition += trimRange;
          readPos += trimRange;
          range -= trimRange;
          trimmedHead += trimRange;
        }

        while ( range > 0 )
        {
          if ( readPos >= readLength - trimTailSize )
          {
            readPos = readLength;
            break;
          }

          char baseBits = strand ? query[readPos] : 3 - query[readLength - readPos - 1];
          char qScore = strand ? qualities[readPos] : qualities[readLength - readPos - 1];
          if ( qScore < 0 || qScore > 40 )
          {
            if ( vcSetting->enableQualityCorrection )
            {
              qScore = qScore < 0 ? 0 : 40;
            }
            else
            {
              fprintf ( stderr, "Error in Base Quality. Please check your input read file.\n" );
              exit ( 1 );
            }
          }
          char recalScore = recalScores[qScore];
          
          char weight = vcSetting->weightMap[recalScore];
          if ( overflowFlag )
          {
            lastMPos = bufferPos;
            buffer[bufferPos++] = ( weight + ( baseBits << 3 ) + 0x40 );
          }
          
          else if ( updateSnpCounter ( snpBundle, vcSetting,
                genomePosition, baseBits, weight, strand, 1 ) == COUNTER_OVERFLOW )
          {
            overflowFlag = 1;
            buffer[bufferPos++] |= genomePosition & 255;
            buffer[bufferPos++] |= ( genomePosition >> 8 ) & 255;
            buffer[bufferPos++] |= ( genomePosition >> 16 ) & 255;
            buffer[bufferPos++] |= ( genomePosition >> 24 ) & 255;
            buffer[bufferPos++] = strand;
            lastMPos = bufferPos;
            buffer[bufferPos++] = ( weight + ( baseBits << 3 ) + 0x40 );  
          }

          
          
          
          ++genomePosition;
          ++readPos;
          --range;
        }
      }
      else if ( cigar[i] == CIGAR_INSERT || cigar[i] == CIGAR_SOFT_CLIP )
      {
        if ( cigar[i] == CIGAR_SOFT_CLIP && range >= vcSetting->softClipThreshold )
        {
          unsigned int clipPos = readPos ? genomePosition - 1 : genomePosition;
          assert ( clipPos < vcSetting->dnaLength );
#ifdef ENABLE_SOFTCLIP_COUNTER
          unsigned char oldCount = 0;
          unsigned char newCount = 0;
          do
          {
            oldCount = snpCounter[clipPos].softClipCount;
            newCount = oldCount + 1;
            if ( newCount == 255 )
            {
#ifdef UPDATE_INVALID_POS
              unsigned int oldBits = 0;
              unsigned int newBits = 0;
              do
              {
                oldBits = invalidPos[clipPos >> 5];
                newBits = oldBits | ( 1 << ( 31 - ( clipPos & 31 ) ) );
              } while ( ! __sync_bool_compare_and_swap ( & ( invalidPos[clipPos >> 5] ), oldBits, newBits ) );
#endif
              break;
            }
          } while ( ! __sync_bool_compare_and_swap ( &( snpCounter[clipPos].softClipCount ), oldCount, newCount ) );
#endif
#ifdef SNP_DIRECTIONAL_SOFTCLIP
          SnpDirectionalSoftclip bufferSoftclip = { clipPos, ( readPos ? 0 : 1) };
          fwrite ( &bufferSoftclip, sizeof ( SnpDirectionalSoftclip ), 1, softclipFilePtr );
#endif
        }
        

        if ( readPos < trimHeadSize )
        {
          readPos += range;
          ++i;
          trimmedHead += range;
          continue;
        }

        if ( readPos + range >= readLength )
        {
          break;
        }

        if ( cigar[i] == CIGAR_INSERT )
        {
          if ( !overflowFlag )
          {
            overflowFlag = 1;
            
            writeFileFlag = 1;
            buffer[bufferPos++] |= genomePosition & 255;
            buffer[bufferPos++] |= ( genomePosition >> 8 ) & 255;
            buffer[bufferPos++] |= ( genomePosition >> 16 ) & 255;
            buffer[bufferPos++] |= ( genomePosition >> 24 ) & 255;
            buffer[bufferPos++] = strand;
          }
          while ( range > 0 )
          {
            char baseBits = strand ? query[readPos] : 3 - query[readLength - readPos - 1];
            char qScore = strand ? qualities[readPos] : qualities[readLength - readPos - 1];
            if ( qScore < 0 || qScore > 40 )
            {
              if ( vcSetting->enableQualityCorrection )
              {
                qScore = qScore < 0 ? 0 : 40;
              }
              else
              {
                fprintf ( stderr, "Error in Base Quality. Please check your input read file.\n" );
                exit ( 1 );
              }
            }
            char recalScore = recalScores[qScore];
            
            char weight = vcSetting->weightMap[recalScore];

            buffer[bufferPos++] = ( weight + ( baseBits << 3 ) + 0x80 );
            ++readPos;
            --range;
          }
        }
        else
        {
          readPos += range;
        }

      }
      else if ( cigar[i] == CIGAR_DELETE )
      {
        if ( readPos < trimHeadSize )
        {
          genomePosition += range;
          ++i;
          continue;
        }

        
        if ( readPos + 1 >= readLength )
        {
          break;
        }

        if ( !overflowFlag )
        {
          overflowFlag = 1;
          
          writeFileFlag = 1;
          buffer[bufferPos++] |= genomePosition & 255;
          buffer[bufferPos++] |= ( genomePosition >> 8 ) & 255;
          buffer[bufferPos++] |= ( genomePosition >> 16 ) & 255;
          buffer[bufferPos++] |= ( genomePosition >> 24 ) & 255;
          buffer[bufferPos++] = strand;
        }
        int range1 = range;
        while ( range1 > 63 )
        {
          buffer[bufferPos++] = ( 0xFF );
          range1 -= 63;
        }
        buffer[bufferPos++] = ( range1 + 0xC0 );

        genomePosition += range;
      }
      ++i;
    }

    if ( overflowFlag || writeFileFlag )
    {
      
      
      if ( lastMPos > 0 )
      {
        buffer[lastMPos] |= 0x20;  
        fwrite ( buffer, 1, lastMPos + 1, overflowFilePtr );
      }
    }
  }
}


char updateSnpCounter ( SnpBundle snpBundle, VcSetting * vcSetting, 
    unsigned int position, char base, char weight, char strand, char strandCount )
{
  SnpCounter * snpCounter = snpBundle.snpCounter;
  SnpOverflowCounterArray * snpOverflowCounterArray = snpBundle.snpOverflowCounterArray;
  unsigned int numOfCPUThreads = snpBundle.numOfCPUThreads;
  unsigned int * invalidPos = snpBundle.invalidSnpCounterPos;

  assert ( position < vcSetting->dnaLength );
  if ( invalidPos[position >> 5] & ( 1 << ( 31 - ( position & 31 ) ) ) )
  {
    return COUNTER_NOT_OVERFLOW;
  }

  if ( ( ( snpCounter[position].weightedCount1 & snpCounter[position].weightedCount2 ) & 0xC000 ) == 0xC000 )
  {
    int region = findRegionByPosition ( position, vcSetting->dnaLength, numOfCPUThreads );

    unsigned short * weightedCounter;
    unsigned short * strandCounter;
    assert ( snpCounter[position].data.arrayIndex < snpOverflowCounterArray[region].limit );
    weightedCounter = & ( snpOverflowCounterArray[region].counters[snpCounter[position].data.arrayIndex].weightedCount[base] );
    strandCounter = strand ?
      & ( snpOverflowCounterArray[region].counters[snpCounter[position].data.arrayIndex].posStrandCount[base] )
      : & ( snpOverflowCounterArray[region].counters[snpCounter[position].data.arrayIndex].negStrandCount[base] );

    unsigned short oldCount = 0;
    unsigned short newCount = 0;
    do
    {
      oldCount = * strandCounter;
      newCount = oldCount + strandCount;
      if ( newCount >= 0xFFFC ) {
#ifdef UPDATE_INVALID_POS 
        unsigned int oldBits = 0;
        unsigned int newBits = 0;
        do
        {
          oldBits = invalidPos[position >> 5];
          newBits = oldBits | ( 1 << ( 31 - ( position & 31 ) ) );
        } while ( ! __sync_bool_compare_and_swap ( & ( invalidPos[position >> 5] ), oldBits, newBits ) );
#endif
        return COUNTER_NOT_OVERFLOW;
      }
    } while ( ! __sync_bool_compare_and_swap ( strandCounter, oldCount, newCount ) );

    do
    {
      oldCount = * weightedCounter;
      newCount = oldCount + weight;
      if ( oldCount >= 0xFFFF - weight ) {
#ifdef UPDATE_INVALID_POS
        
        unsigned int oldBits = 0;
        unsigned int newBits = 0;
        do
        {
          oldBits = invalidPos[position >> 5];
          newBits = oldBits | ( 1 << ( 31 - ( position & 31 ) ) );
        } while ( ! __sync_bool_compare_and_swap ( & ( invalidPos[position >> 5] ), oldBits, newBits ) );
#endif
        return COUNTER_NOT_OVERFLOW;
      }
    } while ( ! __sync_bool_compare_and_swap ( weightedCounter, oldCount, newCount ) );

    return COUNTER_NOT_OVERFLOW;
  }

  unsigned short * weightedCounter;
  unsigned char * strandCounter;

  unsigned char strandCount1 = snpCounter[position].data.strandCounts.posStrandCount1 | snpCounter[position].data.strandCounts.negStrandCount1;
  unsigned char strandCount2 = snpCounter[position].data.strandCounts.posStrandCount2 | snpCounter[position].data.strandCounts.negStrandCount2;
  if ( strandCount1 | strandCount2 )
  {
    if ( strandCount1 && ( ( ( snpCounter[position].weightedCount1 >> 14 ) & 3 ) == base ) )
    {
      weightedCounter = &( snpCounter[position].weightedCount1 );
      strandCounter = strand ? &( snpCounter[position].data.strandCounts.posStrandCount1 ) : &( snpCounter[position].data.strandCounts.negStrandCount1 );        
    }
    else if ( strandCount2 && ( ( ( snpCounter[position].weightedCount2 >> 14 ) & 3 ) == base ) )
    {
      weightedCounter = &( snpCounter[position].weightedCount2 );
      strandCounter = strand ? &( snpCounter[position].data.strandCounts.posStrandCount2 ) : &( snpCounter[position].data.strandCounts.negStrandCount2 );
    }
    else if ( strandCount1 && !strandCount2 )
    {
      snpCounter[position].weightedCount2 |= ( ( ( unsigned short ) base ) << 14 );
      weightedCounter = &( snpCounter[position].weightedCount2 );
      strandCounter = strand ? &( snpCounter[position].data.strandCounts.posStrandCount2 ) : &( snpCounter[position].data.strandCounts.negStrandCount2 );
    }
    else if ( !strandCount1 && strandCount2 )
    {
      snpCounter[position].weightedCount1 |= ( ( ( unsigned short ) base ) << 14 );
      weightedCounter = &( snpCounter[position].weightedCount1 );
      strandCounter = strand ? &( snpCounter[position].data.strandCounts.posStrandCount1 ) : &( snpCounter[position].data.strandCounts.negStrandCount1 );
    }
    else
    {
      return COUNTER_OVERFLOW;        
    }
  }
  else
  {
    
    snpCounter[position].weightedCount1 |= ( ( ( unsigned short ) base ) << 14 );
    weightedCounter = &( snpCounter[position].weightedCount1 );
    strandCounter = strand ? &( snpCounter[position].data.strandCounts.posStrandCount1 ) : &( snpCounter[position].data.strandCounts.negStrandCount1 );
    snpCounter[position].weightedCount2 = 0; 
  }

  unsigned char oldStrandCount = 0;
  unsigned char newStrandCount = 0;
  do
  {
    oldStrandCount = * strandCounter;
    newStrandCount = oldStrandCount + strandCount;
    if ( newStrandCount == 255 ) {
      return COUNTER_OVERFLOW;
    }
  } while ( ! __sync_bool_compare_and_swap ( strandCounter, oldStrandCount, newStrandCount ) );

  unsigned short oldCount = 0;
  unsigned short newCount = 0;
  do
  {
    oldCount = * weightedCounter;
    newCount = ( oldCount & 0x3FFF ) + weight;
    if ( newCount >= 0x3FFF ) {
#ifdef UPDATE_INVALID_POS
      
      unsigned int oldBits = 0;
      unsigned int newBits = 0;
      do
      {
        oldBits = invalidPos[position >> 5];
        newBits = oldBits | ( 1 << ( 31 - ( position & 31 ) ) );
      } while ( ! __sync_bool_compare_and_swap ( & ( invalidPos[position >> 5] ), oldBits, newBits ) );
#endif
      break;
    }
  } while ( ! __sync_bool_compare_and_swap ( weightedCounter, oldCount , newCount | ( ( ( unsigned short ) base ) << 14 ) ) );

  

  return COUNTER_NOT_OVERFLOW;
}

void readFilesForOverflow ( SnpBundle snpBundle, VcSetting * vcSetting,
    const char * filename )
{
  unsigned int numOfCPUThreads = snpBundle.numOfCPUThreads;

  FILE * file = ( FILE * ) fopen ( filename, "rb" );
  unsigned int totalSnpOverflowBufferSize = 0;
  unsigned int bufferSize = OVERFLOW_BUFFER_SIZE;

  unsigned char buffer[OVERFLOW_BUFFER_SIZE];
  int bufferPos = 0;
  unsigned int bRead = fread ( buffer, 1, bufferSize, file );

  char isStartBatch = 1;
  unsigned int startPos = 0;
  char strand;
  unsigned char baseInfo;
  
  
  
  char nextBasePos = 0;

  char isDelete = 0;
  unsigned int deletePos = 0;
  unsigned char deleteLength = 0;

  float totalWeight = 0.0f;
  unsigned char insertLength = 0;
  unsigned char insertSeq[MAX_READ_LENGTH];
  memset ( insertSeq, 0, MAX_READ_LENGTH );

  while ( bRead > 0 )
  {
    if ( isStartBatch )
    {
      unsigned char posBuffer[4];
      int i;
      for ( i = 0; i < 4; ++i )
      {
        posBuffer[i] = buffer[bufferPos++];
        if ( bufferPos >= bRead )
        {
          bRead = fread ( buffer, 1, bufferSize, file );
          bufferPos = 0;
        }
      }
      startPos = ( posBuffer[0]
          + ( ( ( unsigned int ) posBuffer[1] ) << 8 )
          + ( ( ( unsigned int ) posBuffer[2] ) << 16 )
          + ( ( ( unsigned int ) posBuffer[3] ) << 24 ) );

      strand = buffer[bufferPos++];
      if ( bufferPos >= bRead )
      {
        bRead = fread ( buffer, 1, bufferSize, file );
        bufferPos = 0;
      }

      isStartBatch = 0;
    }

    baseInfo = buffer[bufferPos++];

    nextBasePos = 0;

    unsigned char infoType = ( ( baseInfo & 0xC0 ) >> 6 );

    switch ( infoType )
    {
      case 1:
        {
          if ( insertLength > 0 )
          {
            addSnpOverflowBufferForInsertionFromInfo ( vcSetting->dnaLength, numOfCPUThreads, startPos, strand, 
                (unsigned char) ( totalWeight / insertLength ),
                insertLength, insertSeq );

            insertLength = 0;
            totalWeight = 0.0f;
            memset ( insertSeq, 0, MAX_READ_LENGTH );
          }

          if ( isDelete )
          {
            addSnpOverflowBufferForDeletionFromInfo ( vcSetting->dnaLength, numOfCPUThreads, deletePos, strand, baseInfo & 7, deleteLength );
            isDelete = 0;
            deleteLength = 0;
          }

          addSnpOverflowBufferFromInfo ( vcSetting->dnaLength, numOfCPUThreads, startPos, strand, baseInfo );
          nextBasePos = 1;
          break;
        }
      case 2:
        {
          if ( isDelete )
          {
            fprintf ( stderr, "Warning in parsing alignment information : deletion before insertion\n" );
            addSnpOverflowBufferForDeletionFromInfo ( vcSetting->dnaLength, numOfCPUThreads, deletePos, strand, baseInfo & 7, deleteLength );
            isDelete = 0;
            deleteLength = 0;
          }

          char base = ( baseInfo >> 3 ) & 3;
          totalWeight += ( baseInfo & 7 );
          insertSeq[insertLength >> 2] |= ( base << ( ( insertLength & 3 ) << 1 ) );
          insertLength++;
          break;
        }
      case 3:
        {
          if ( isDelete )
          {
            nextBasePos = baseInfo & 63;
            deleteLength += nextBasePos;
            break;
          }
          else if ( insertLength > 0 )
          {
            addSnpOverflowBufferForInsertionFromInfo ( vcSetting->dnaLength, numOfCPUThreads, startPos, strand, 
                (unsigned char) ( totalWeight / insertLength ),
                insertLength, insertSeq );

            insertLength = 0;
            totalWeight = 0.0f;
            memset ( insertSeq, 0, MAX_READ_LENGTH );
          }

          
          deletePos = startPos;
          isDelete = 1;
          nextBasePos = baseInfo & 63;
          deleteLength = nextBasePos;
          break;
        }
      default:
        {
          fprintf ( stderr, "Error in handling alignment information in Overflow : invalid base type [%d]\n", infoType );
          exit ( 1 );
        }
    }

    if ( ++totalSnpOverflowBufferSize >= OVERFLOW_BATCH_SIZE )
    {
      startUpdateOverflowCounter ( snpBundle, vcSetting );
      emptySnpOverflowBuffer ( numOfCPUThreads );
      totalSnpOverflowBufferSize = 0;
    }

    
    
    if ( infoType == 1 && ( baseInfo & 0x20 ) )  
    {
      isStartBatch = 1;
    }
    else
    {
      
      startPos += nextBasePos;
    }

    if ( startPos >= vcSetting->dnaLength )
    {
#ifndef DISABLE_SNP_WARNING
      fprintf ( stderr, "Invalid postiion in readFilesForOverflow(): %u", startPos );
#endif
      break;
    }

    if ( bufferPos >= bRead )
    {
      bRead = fread ( buffer, 1, bufferSize, file );
      bufferPos = 0;
    }
  }
  fclose ( file );

  startUpdateOverflowCounter ( snpBundle, vcSetting );
  emptySnpOverflowBuffer ( numOfCPUThreads );
}

void updateSnpCounterForOverflow ( SnpBundle snpBundle, VcSetting * vcSetting,
    unsigned int position, char base, char weight,
    char strand, char strandCount, unsigned int threadId )
{
  SnpCounter * snpCounter = snpBundle.snpCounter;
  SnpOverflowCounterArray * snpOverflowCounterArray = snpBundle.snpOverflowCounterArray;
  unsigned int numOfCPUThreads = snpBundle.numOfCPUThreads;
  unsigned int * invalidPos = snpBundle.invalidSnpCounterPos;

  assert ( position < vcSetting->dnaLength );
  if ( invalidPos[position >> 5] & ( 1 << ( 31 - ( position & 31 ) ) ) )
  {
    return;
  }

  if ( ( ( snpCounter[position].weightedCount1 & snpCounter[position].weightedCount2 ) & 0xC000 ) != 0xC000 )
  {
    if ( updateSnpCounter ( snpBundle, vcSetting,
          position, base, weight, strand, strandCount ) == COUNTER_NOT_OVERFLOW )
    {
      return;
    }
    unsigned int arrayIndex = addSnpOverflowCounterToArray ( snpOverflowCounterArray[threadId], numOfCPUThreads );
    assert ( arrayIndex < snpOverflowCounterArray[threadId].limit );

    char baseBit = ( snpCounter[position].weightedCount1 >> 14 ) & 3;
    snpOverflowCounterArray[threadId].counters[arrayIndex].weightedCount[baseBit] = ( snpCounter[position].weightedCount1 & 0x00003FFF );
    snpOverflowCounterArray[threadId].counters[arrayIndex].posStrandCount[baseBit] = snpCounter[position].data.strandCounts.posStrandCount1;
    snpOverflowCounterArray[threadId].counters[arrayIndex].negStrandCount[baseBit] = snpCounter[position].data.strandCounts.negStrandCount1;

    if ( snpCounter[position].data.strandCounts.posStrandCount2 | snpCounter[position].data.strandCounts.negStrandCount2 )
    {
      baseBit = ( snpCounter[position].weightedCount2 >> 14 ) & 3;
      snpOverflowCounterArray[threadId].counters[arrayIndex].weightedCount[baseBit] = ( snpCounter[position].weightedCount2 & 0x00003FFF );
      snpOverflowCounterArray[threadId].counters[arrayIndex].posStrandCount[baseBit] = snpCounter[position].data.strandCounts.posStrandCount2;
      snpOverflowCounterArray[threadId].counters[arrayIndex].negStrandCount[baseBit] = snpCounter[position].data.strandCounts.negStrandCount2;
    }

    
    snpCounter[position].weightedCount1 = 0xC000;
    snpCounter[position].weightedCount2 = 0xC000;
    snpCounter[position].data.arrayIndex = arrayIndex;
  }
  unsigned short * weightedCounter = & ( snpOverflowCounterArray[threadId].counters[snpCounter[position].data.arrayIndex].weightedCount[base] );
  unsigned short * strandCounter = strand ? 
    & ( snpOverflowCounterArray[threadId].counters[snpCounter[position].data.arrayIndex].posStrandCount[base] )
    : & ( snpOverflowCounterArray[threadId].counters[snpCounter[position].data.arrayIndex].negStrandCount[base] );

  if ( ( 0xFFFF - ( * weightedCounter ) < weight ) || ( ( * strandCounter ) + strandCount >= 0xFFFC ) )
  {
#ifdef UPDATE_INVALID_POS
    unsigned int oldBits = 0;
    unsigned int newBits = 0;
    do
    {
      oldBits = invalidPos[position >> 5];
      newBits = oldBits | ( 1 << ( 31 - ( position & 31 ) ) );
    } while ( ! __sync_bool_compare_and_swap ( & ( invalidPos[position >> 5] ), oldBits, newBits ) );
#endif
    return;
  }

  * weightedCounter += weight;
  * strandCounter += strandCount;
}

void updateSnpCounterForInsert ( SnpBundle snpBundle, VcSetting * vcSetting,
    unsigned int position, unsigned char insertLength, unsigned char avgWeight,
    unsigned char * insertSeq, unsigned int threadId )
{
  SnpCounter * snpCounter = snpBundle.snpCounter;
  SnpOverflowCounterArray * snpOverflowCounterArray = snpBundle.snpOverflowCounterArray;
  unsigned int numOfCPUThreads = snpBundle.numOfCPUThreads;
  unsigned int * invalidPos = snpBundle.invalidSnpCounterPos;
  MemoryPool * pool = snpBundle.snpMemoryPool;

  if ( invalidPos[position >> 5] & ( 1 << ( 31 - ( position & 31 ) ) ) )
  {
    return;
  }

  unsigned char hqCount = 0;
  unsigned char lqCount = 0;
  if ( avgWeight > vcSetting->indelWeightThreshold )
  {
    hqCount = 1;
  }
  else
  {
    lqCount = 1;
  }

  if ( ( ( snpCounter[position].weightedCount1 & snpCounter[position].weightedCount2 ) & 0xC000 ) != 0xC000 )
  {
    unsigned int arrayIndex = addSnpOverflowCounterToArray ( snpOverflowCounterArray[threadId], numOfCPUThreads );
    assert ( arrayIndex < snpOverflowCounterArray[threadId].limit );

    char baseBit = ( snpCounter[position].weightedCount1 >> 14 ) & 3;
    snpOverflowCounterArray[threadId].counters[arrayIndex].weightedCount[baseBit] = ( snpCounter[position].weightedCount1 & 0x00003FFF );
    snpOverflowCounterArray[threadId].counters[arrayIndex].posStrandCount[baseBit] = snpCounter[position].data.strandCounts.posStrandCount1;
    snpOverflowCounterArray[threadId].counters[arrayIndex].negStrandCount[baseBit] = snpCounter[position].data.strandCounts.negStrandCount1;

    if ( snpCounter[position].data.strandCounts.posStrandCount2 | snpCounter[position].data.strandCounts.negStrandCount2 )
    {
      baseBit = ( snpCounter[position].weightedCount2 >> 14 ) & 3;
      snpOverflowCounterArray[threadId].counters[arrayIndex].weightedCount[baseBit] = ( snpCounter[position].weightedCount2 & 0x00003FFF );
      snpOverflowCounterArray[threadId].counters[arrayIndex].posStrandCount[baseBit] = snpCounter[position].data.strandCounts.posStrandCount2;
      snpOverflowCounterArray[threadId].counters[arrayIndex].negStrandCount[baseBit] = snpCounter[position].data.strandCounts.negStrandCount2;
    }

    
    snpCounter[position].weightedCount1 = 0xC000;
    snpCounter[position].weightedCount2 = 0xC000;
    snpCounter[position].data.arrayIndex = arrayIndex;
  }

  unsigned int arrayIndex = snpCounter[position].data.arrayIndex;
  assert ( arrayIndex < snpOverflowCounterArray[threadId].limit );

  if ( !( snpCounter[position].weightedCount1 & 0x2000 ) )
    
  {
    snpCounter[position].weightedCount1 |= 0x2000;

    if ( insertLength <= 12 )
    {
      snpCounter[position].weightedCount1 |= ( ( ( ( unsigned short ) insertLength ) << 8 ) + hqCount );
      snpOverflowCounterArray[threadId].counters[arrayIndex].insertion.counters.lqCount = lqCount;
      for ( int i = 0; i < ( ( insertLength + 3 ) / 4 ); ++i )
      {
        snpOverflowCounterArray[threadId].counters[arrayIndex].insertion.counters.insertSeq[i] = insertSeq[i];
      }
    }
    else
    {
      snpCounter[position].weightedCount1 |= 0x1000;
      unsigned int ptr;
      unsigned int numBytes = ( insertLength + 3 ) / 4 + 5;
      int lz = __builtin_clz ( numBytes );
      unsigned int numBytes2 = ( 1 << ( 32 - lz - 1 ) );
      if ( numBytes2 < numBytes ) {
        numBytes2 = numBytes2 << 1;
      }
      ptr = pmallocIndex ( pool, numBytes2 * sizeof(char) );

      unsigned char * insertInfo = ( unsigned char * ) getAddress ( pool, ptr );
      int infoIndex = 0;

      insertInfo[infoIndex++] = insertLength;
      for ( int i = 0; i < ( ( insertLength + 3 ) / 4 ); ++i )
      {
        insertInfo[infoIndex++] = insertSeq[i];
      }
      insertInfo[infoIndex++] = hqCount & 0xff;
      insertInfo[infoIndex++] = ( hqCount >> 8 ) & 0xff;

      insertInfo[infoIndex++] = lqCount & 0xff;
      insertInfo[infoIndex++] = ( lqCount >> 8 ) & 0xff;

      snpCounter[position].weightedCount1 |= infoIndex;

      snpOverflowCounterArray[threadId].counters[arrayIndex].insertion.ptr = ptr;
    }
  }
  else
    
  {
    if ( !( snpCounter[position].weightedCount1 & 0x1000 ) )
      
    {

      unsigned char insertLength0 = ( ( snpCounter[position].weightedCount1 >> 8 ) & 0xF );
      char sameSeq = 0;
      if ( insertLength0 == insertLength )
      {    
        sameSeq = 1;
        for ( int i = 0; i < ( ( insertLength + 3 ) / 4 ); ++i )
        {
          if ( snpOverflowCounterArray[threadId].counters[arrayIndex].insertion.counters.insertSeq[i] != insertSeq[i] )
          {
            sameSeq = 0;
            break;
          }
        }

        if ( sameSeq 
            && ( ( unsigned int ) ( snpCounter[position].weightedCount1 & 0xFF ) ) + hqCount < 0xFF 
            && ( ( unsigned int ) snpOverflowCounterArray[threadId].counters[arrayIndex].insertion.counters.lqCount ) + lqCount < 0xFF )
        {
          
          snpCounter[position].weightedCount1 += hqCount;
          snpOverflowCounterArray[threadId].counters[arrayIndex].insertion.counters.lqCount += lqCount;

          return;
        }
      }

      
      unsigned ptr;
      unsigned int numBytes = ( ( insertLength0 + 3 ) / 4 + 5 ) + ( sameSeq ? 0 : ( ( insertLength + 3 ) / 4 + 5 ) );
      int lz = __builtin_clz ( numBytes );
      unsigned int numBytes2 = ( 1 << ( 32 - lz - 1 ) );
      if ( numBytes2 < numBytes ) {
        numBytes2 = numBytes2 << 1;
      }
      ptr = pmallocIndex ( pool, numBytes2 * sizeof(char) );

      unsigned char * insertInfo = ( unsigned char * ) getAddress ( pool, ptr );
      int infoIndex = 0;
      
      insertInfo[infoIndex++] = insertLength0;
      for ( int i = 0; i < ( ( insertLength0 + 3 ) / 4 ); ++i )
      {
        insertInfo[infoIndex++] = snpOverflowCounterArray[threadId].counters[arrayIndex].insertion.counters.insertSeq[i];
      }
      unsigned short mphqCount;
      mphqCount = ( snpCounter[position].weightedCount1 & 0xff ) + ( sameSeq ? hqCount : 0 );
      insertInfo[infoIndex++] = mphqCount & 0xff;
      insertInfo[infoIndex++] = ( mphqCount >> 8 ) & 0xff; 

      unsigned short mplqCount;
      mplqCount = snpOverflowCounterArray[threadId].counters[arrayIndex].insertion.counters.lqCount + ( sameSeq ? lqCount : 0 );
      insertInfo[infoIndex++] = mplqCount & 0xff;
      insertInfo[infoIndex++] = ( mplqCount >> 8 ) & 0xff; 

      if ( !sameSeq )
      {
        
        insertInfo[infoIndex++] = insertLength;
        for ( int i = 0; i < ( ( insertLength + 3 ) / 4 ); ++i )
        {
          insertInfo[infoIndex++] = insertSeq[i];
        }
        insertInfo[infoIndex++] = hqCount & 0xff;
        insertInfo[infoIndex++] = ( hqCount >> 8 ) & 0xff;
        insertInfo[infoIndex++] = lqCount & 0xff;
        insertInfo[infoIndex++] = ( lqCount >> 8 ) & 0xff;
      }
      snpCounter[position].weightedCount1 = ( 0xF000 | infoIndex );

      snpOverflowCounterArray[threadId].counters[arrayIndex].insertion.ptr = ptr;

    }
    else
      
    {
      unsigned short arraySize = ( snpCounter[position].weightedCount1 & 0xFFF );
      unsigned char * insertInfo = ( unsigned char * ) getAddress ( pool, snpOverflowCounterArray[threadId].counters[arrayIndex].insertion.ptr );
      int infoIndex = 0;
      while ( infoIndex < arraySize )
      {
        unsigned char size = insertInfo[infoIndex++];

        if ( insertLength == size )
        {
          char sameSeq = 1;
          for ( int i = 0; i < ( size + 3 ) / 4; ++i )
          {
            if ( insertInfo[infoIndex + i] != insertSeq[i] )
            {
              sameSeq = 0;
              break;
            }
          }
          if ( sameSeq )
          {
            infoIndex += ( size + 3 ) / 4;
            unsigned short mphqCount;
            mphqCount = insertInfo[infoIndex] & 0xff;
            mphqCount |= insertInfo[infoIndex+1] << 8;
            if (  mphqCount >= 0xFFFF - hqCount )
            {
#ifdef UPDATE_INVALID_POS
              unsigned int oldBits = 0;
              unsigned int newBits = 0;
              do
              {
                oldBits = invalidPos[position >> 5];
                newBits = oldBits | ( 1 << ( 31 - ( position & 31 ) ) );
              } while ( ! __sync_bool_compare_and_swap ( & ( invalidPos[position >> 5] ), oldBits, newBits ) );
#endif
              return;
            }
            mphqCount += hqCount;
            insertInfo[infoIndex++] = mphqCount & 0xff;
            insertInfo[infoIndex++] = ( mphqCount >> 8 ) & 0xff;

            unsigned short mplqCount;
            mplqCount = insertInfo[infoIndex] & 0xff;
            mplqCount |= insertInfo[infoIndex+1] << 8; 
            if ( mplqCount >= 0xFFFF - lqCount )
            {
#ifdef UPDATE_INVALID_POS
              unsigned int oldBits = 0;
              unsigned int newBits = 0;
              do
              {
                oldBits = invalidPos[position >> 5];
                newBits = oldBits | ( 1 << ( 31 - ( position & 31 ) ) );
              } while ( ! __sync_bool_compare_and_swap ( & ( invalidPos[position >> 5] ), oldBits, newBits ) );
#endif
              return;
            }
            mplqCount += lqCount;
            insertInfo[infoIndex++] = mplqCount & 0xff;
            insertInfo[infoIndex++] = ( mplqCount >> 8 ) & 0xff;
            return;
          }
        }

        infoIndex += ( size + 3 ) / 4 + 4;
      }

      
      unsigned int ptr;
      int lz0 = __builtin_clz ( arraySize );
      unsigned int curBytes2 = ( 1 << ( 32 - lz0 - 1 ) );
      if ( curBytes2 < arraySize )
      {
        curBytes2 = curBytes2 << 1;
      }
      unsigned int numBytes = ( insertLength + 3 ) / 4 + 5 + arraySize;
      int lz = __builtin_clz ( numBytes );
      unsigned int numBytes2 = ( 1 << ( 32 - lz - 1 ) );
      if ( numBytes2 < numBytes )
      {
        numBytes2 = numBytes2 << 1;
      }

      if ( numBytes2 > 4096 )
      {
#ifdef UPDATE_INVALID_POS
        unsigned int oldBits = 0;
        unsigned int newBits = 0;
        do
        {
          oldBits = invalidPos[position >> 5];
          newBits = oldBits | ( 1 << ( 31 - ( position & 31 ) ) );
        } while ( ! __sync_bool_compare_and_swap ( & ( invalidPos[position >> 5] ), oldBits, newBits ) );
#endif
        return;
      }

      if ( numBytes2 > curBytes2 )
      {
        
        ptr = premallocIndex ( pool, getAddress ( pool, snpOverflowCounterArray[threadId].counters[arrayIndex].insertion.ptr ),
            numBytes2 * sizeof(char) );
        unsigned char * newInsertInfo = ( unsigned char * ) getAddress ( pool, ptr );
        memcpy ( newInsertInfo, insertInfo, infoIndex );
        insertInfo = newInsertInfo;

        snpOverflowCounterArray[threadId].counters[arrayIndex].insertion.ptr = ptr;
      }

      insertInfo[infoIndex++] = insertLength;
      for ( int i = 0; i < ( insertLength + 3 ) / 4; ++i )
      {
        insertInfo[infoIndex++] = insertSeq[i];
      }
      insertInfo[infoIndex++] = hqCount & 0xff;
      insertInfo[infoIndex++] = ( hqCount >> 8 ) & 0xff;

      insertInfo[infoIndex++] = lqCount & 0xff;
      insertInfo[infoIndex++] = ( lqCount >> 8 ) & 0xff; 

      if ( infoIndex >= 4096 )
      {
#ifdef UPDATE_INVALID_POS
        unsigned int oldBits = 0;
        unsigned int newBits = 0;
        do
        {
          oldBits = invalidPos[position >> 5];
          newBits = oldBits | ( 1 << ( 31 - ( position & 31 ) ) );
        } while ( ! __sync_bool_compare_and_swap ( & ( invalidPos[position >> 5] ), oldBits, newBits ) );
#endif
        return;
      }

      snpCounter[position].weightedCount1 = ( 0xF000 | infoIndex );
    }
  }
}

void updateSnpCounterForDelete ( SnpBundle snpBundle, VcSetting * vcSetting,
    unsigned int position, unsigned char deleteLength,
    char weight, unsigned int threadId )
{
  SnpCounter * snpCounter = snpBundle.snpCounter;
  SnpOverflowCounterArray * snpOverflowCounterArray = snpBundle.snpOverflowCounterArray;
  unsigned int numOfCPUThreads = snpBundle.numOfCPUThreads;
  unsigned int * invalidPos = snpBundle.invalidSnpCounterPos;
  MemoryPool * pool = snpBundle.snpMemoryPool;

  if ( invalidPos[position >> 5] & ( 1 << ( 31 - ( position & 31 ) ) ) )
  {
    return;
  }

  unsigned char hqCount = 0;
  unsigned char lqCount = 0;
  if ( weight > vcSetting->indelWeightThreshold )
  {
    hqCount = 1;
  }
  else
  {
    lqCount = 1;
  }

  if ( ( ( snpCounter[position].weightedCount1 & snpCounter[position].weightedCount2 ) & 0xC000 ) != 0xC000 )
  {
    unsigned int arrayIndex = addSnpOverflowCounterToArray ( snpOverflowCounterArray[threadId], numOfCPUThreads );
    assert ( arrayIndex < snpOverflowCounterArray[threadId].limit );

    char baseBit = ( snpCounter[position].weightedCount1 >> 14 ) & 3;
    snpOverflowCounterArray[threadId].counters[arrayIndex].weightedCount[baseBit] = ( snpCounter[position].weightedCount1 & 0x00003FFF );
    snpOverflowCounterArray[threadId].counters[arrayIndex].posStrandCount[baseBit] = snpCounter[position].data.strandCounts.posStrandCount1;
    snpOverflowCounterArray[threadId].counters[arrayIndex].negStrandCount[baseBit] = snpCounter[position].data.strandCounts.negStrandCount1;

    if ( snpCounter[position].data.strandCounts.posStrandCount2 | snpCounter[position].data.strandCounts.negStrandCount2 )
    {
      baseBit = ( snpCounter[position].weightedCount2 >> 14 ) & 3;
      snpOverflowCounterArray[threadId].counters[arrayIndex].weightedCount[baseBit] = ( snpCounter[position].weightedCount2 & 0x00003FFF );
      snpOverflowCounterArray[threadId].counters[arrayIndex].posStrandCount[baseBit] = snpCounter[position].data.strandCounts.posStrandCount2;
      snpOverflowCounterArray[threadId].counters[arrayIndex].negStrandCount[baseBit] = snpCounter[position].data.strandCounts.negStrandCount2;
    }

    
    snpCounter[position].weightedCount1 = 0xC000;
    snpCounter[position].weightedCount2 = 0xC000;
    snpCounter[position].data.arrayIndex = arrayIndex;
  }

  unsigned int arrayIndex = snpCounter[position].data.arrayIndex;
  assert ( arrayIndex < snpOverflowCounterArray[threadId].limit );

  if ( !( snpCounter[position].weightedCount2 & 0x2000 ) )
    
  {
    snpCounter[position].weightedCount2 |= 0x2000;

    snpOverflowCounterArray[threadId].counters[arrayIndex].deletion.counters.deletionLength = deleteLength;
    snpOverflowCounterArray[threadId].counters[arrayIndex].deletion.counters.hqCount = hqCount;
    snpOverflowCounterArray[threadId].counters[arrayIndex].deletion.counters.lqCount = lqCount;
  }
  else
    
  {
    if ( !( snpCounter[position].weightedCount2 & 0x1000 ) )
      
    {
      char sameDel = ( snpOverflowCounterArray[threadId].counters[arrayIndex].deletion.counters.deletionLength == deleteLength );
      if ( sameDel 
          && ( ( unsigned int ) snpOverflowCounterArray[threadId].counters[arrayIndex].deletion.counters.hqCount ) + hqCount < 0xFF
          && ( ( unsigned int ) snpOverflowCounterArray[threadId].counters[arrayIndex].deletion.counters.lqCount ) + lqCount < 0xFF )
      {
        snpOverflowCounterArray[threadId].counters[arrayIndex].deletion.counters.hqCount += hqCount;
        snpOverflowCounterArray[threadId].counters[arrayIndex].deletion.counters.lqCount += lqCount;

        return;
      }

      snpCounter[position].weightedCount2 |= 0x1000;

      unsigned int ptr;
      ptr = pmallocIndex ( pool, 16 * sizeof(char) );
      unsigned char * deleteInfo = ( unsigned char * ) getAddress ( pool, ptr );
      int infoIndex = 0;
      deleteInfo[infoIndex++] = snpOverflowCounterArray[threadId].counters[arrayIndex].deletion.counters.deletionLength;
      unsigned short mphqCount = snpOverflowCounterArray[threadId].counters[arrayIndex].deletion.counters.hqCount + ( sameDel ? hqCount : 0 );
      deleteInfo[infoIndex++] = mphqCount & 0xff;
      deleteInfo[infoIndex++] = ( mphqCount >> 8 ) & 0xff;
      unsigned short mplqCount = snpOverflowCounterArray[threadId].counters[arrayIndex].deletion.counters.lqCount + ( sameDel ? lqCount : 0 );
      deleteInfo[infoIndex++] = mplqCount & 0xff;
      deleteInfo[infoIndex++] = ( mplqCount >> 8 ) & 0xff;

      if ( !sameDel )
      {
        deleteInfo[infoIndex++] = deleteLength;
        deleteInfo[infoIndex++] = hqCount & 0xff;
        deleteInfo[infoIndex++] = ( hqCount >> 8 ) & 0xff;
        deleteInfo[infoIndex++] = lqCount & 0xff;
        deleteInfo[infoIndex++] = ( lqCount >> 8 ) & 0xff;
      }

      snpOverflowCounterArray[threadId].counters[arrayIndex].deletion.ptr = ptr;
      snpCounter[position].weightedCount2 |= infoIndex;

    }
    else
      
    {
      unsigned short arraySize = ( snpCounter[position].weightedCount2 & 0xFFF );
      unsigned char * deleteInfo = ( unsigned char * ) getAddress ( pool, snpOverflowCounterArray[threadId].counters[arrayIndex].deletion.ptr );
      int infoIndex = 0;
      while ( infoIndex < arraySize )
      {
        unsigned char length = deleteInfo[infoIndex++];
        if ( length == deleteLength )
        {
          unsigned short mphqCount;
          mphqCount = deleteInfo[infoIndex] & 0xff;
          mphqCount |= deleteInfo[infoIndex+1] << 8;
          if ( mphqCount >= 0xFFFF - hqCount )
          {
#ifdef UPDATE_INVALID_POS
            unsigned int oldBits = 0;
            unsigned int newBits = 0;
            do
            {
              oldBits = invalidPos[position >> 5];
              newBits = oldBits | ( 1 << ( 31 - ( position & 31 ) ) );
            } while ( ! __sync_bool_compare_and_swap ( & ( invalidPos[position >> 5] ), oldBits, newBits ) );
#endif
            return;
          }
          mphqCount += hqCount;
          deleteInfo[infoIndex++] = mphqCount & 0xff;
          deleteInfo[infoIndex++] = ( mphqCount >> 8 ) & 0xff;

          unsigned short mplqCount;
          mplqCount = deleteInfo[infoIndex] & 0xff;
          mplqCount |= deleteInfo[infoIndex+1] << 8;
          if ( mplqCount >= 0xFFFF - lqCount )
          {
#ifdef UPDATE_INVALID_POS
            unsigned int oldBits = 0;
            unsigned int newBits = 0;
            do
            {
              oldBits = invalidPos[position >> 5];
              newBits = oldBits | ( 1 << ( 31 - ( position & 31 ) ) );
            } while ( ! __sync_bool_compare_and_swap ( & ( invalidPos[position >> 5] ), oldBits, newBits ) );
#endif
            return;
          }
          mplqCount += lqCount;
          deleteInfo[infoIndex++] = mplqCount & 0xff;
          deleteInfo[infoIndex++] = ( mplqCount >> 8 ) & 0xff;
          return;
        }
        else
        {
          infoIndex += 4;
        }
      }

      
      unsigned int ptr;
      int lz0 = __builtin_clz ( arraySize );
      unsigned int curBytes2 = ( 1 << ( 32 - lz0 - 1 ) );
      if ( curBytes2 < arraySize )
      {
        curBytes2 = curBytes2 << 1;
      }
      unsigned int numBytes = arraySize + 5;
      int lz = __builtin_clz ( numBytes );
      unsigned int numBytes2 = ( 1 << ( 32 - lz - 1 ) );
      if ( numBytes2 < numBytes )
      {
        numBytes2 = numBytes2 << 1;
      }

      if ( numBytes2 > 4096 )
      {
#ifdef UPDATE_INVALID_POS
        unsigned int oldBits = 0;
        unsigned int newBits = 0;
        do
        {
          oldBits = invalidPos[position >> 5];
          newBits = oldBits | ( 1 << ( 31 - ( position & 31 ) ) );
        } while ( ! __sync_bool_compare_and_swap ( & ( invalidPos[position >> 5] ), oldBits, newBits ) );
#endif
        return;
      }

      if ( curBytes2 < numBytes2 )
      {
        
        ptr = premallocIndex ( pool, getAddress ( pool, snpOverflowCounterArray[threadId].counters[arrayIndex].deletion.ptr ),
            numBytes2 * sizeof(char) );
        unsigned char * newDeleteInfo = ( unsigned char * ) getAddress ( pool, ptr );
        memcpy ( newDeleteInfo, deleteInfo, infoIndex );
        deleteInfo = newDeleteInfo;

        snpOverflowCounterArray[threadId].counters[arrayIndex].deletion.ptr = ptr;
      }

      deleteInfo[infoIndex++] = deleteLength;
      deleteInfo[infoIndex++] = hqCount & 0xff;
      deleteInfo[infoIndex++] = ( hqCount >> 8 ) & 0xff;
      deleteInfo[infoIndex++] = lqCount & 0xff;
      deleteInfo[infoIndex++] = ( lqCount >> 8 ) & 0xff;

      if ( infoIndex >= 4096 )
      {
#ifdef UPDATE_INVALID_POS
        unsigned int oldBits = 0;
        unsigned int newBits = 0;
        do
        {
          oldBits = invalidPos[position >> 5];
          newBits = oldBits | ( 1 << ( 31 - ( position & 31 ) ) );
        } while ( ! __sync_bool_compare_and_swap ( & ( invalidPos[position >> 5] ), oldBits, newBits ) );
#endif
        return;
      }

      snpCounter[position].weightedCount2 = ( 0xF000 | infoIndex );
    }
  }
}

void updateOverflowCounterFromBuffer ( SnpBundle snpBundle, VcSetting * vcSetting,
    unsigned int threadId )
{
  SnpOverflowBufferArray * bufferArray = snpOverflowBufferArray[threadId];
  if ( bufferArray == NULL )
  {
    fprintf ( stderr, "Error in Overflow Buffer\n" );
  }

  unsigned int i;
  for ( i = 0; i < bufferArray->size; ++i )
  {
    char type;
    unsigned char baseInfo = bufferArray->buffer[i].info;
    type = ( baseInfo >> 6 ) & 3;

    if ( type == 1 )
    {
      char weight, baseBits;
      weight = baseInfo & 7;
      baseBits = ( baseInfo >> 3 ) & 3;
      updateSnpCounterForOverflow ( snpBundle, vcSetting,
          bufferArray->buffer[i].position, baseBits, weight,
          bufferArray->buffer[i].strand, 1, threadId );
    }
    else if ( type == 2 )
    {
      updateSnpCounterForInsert ( snpBundle, vcSetting,
          bufferArray->buffer[i].position, bufferArray->buffer[i].length,
          bufferArray->buffer[i].weight, bufferArray->buffer[i].insertSeq,
          threadId );
    }
    else if ( type == 3 )
    {
      updateSnpCounterForDelete ( snpBundle, vcSetting,
          bufferArray->buffer[i].position, bufferArray->buffer[i].length,
          bufferArray->buffer[i].weight, threadId );
    }
  }
}

void * updateOverflowCounterFromBufferWrapper ( void * ptr )
{
  SnpUpdateOverflowWrapperObj * obj = ( SnpUpdateOverflowWrapperObj * ) ptr;

  updateOverflowCounterFromBuffer ( obj->snpBundle, obj->vcSetting,
      obj->threadId );

  return NULL;
}

void startUpdateOverflowCounter ( SnpBundle snpBundle, VcSetting * vcSetting )
{
  unsigned int numOfCPUThreads = snpBundle.numOfCPUThreads;

  pthread_t thread[numOfCPUThreads];

  SnpUpdateOverflowWrapperObj ** obj = ( SnpUpdateOverflowWrapperObj ** ) malloc ( numOfCPUThreads * sizeof(SnpUpdateOverflowWrapperObj *) );

  unsigned int threadId;
  for ( threadId = 0; threadId < numOfCPUThreads; ++threadId )
  {
    obj[threadId] = ( SnpUpdateOverflowWrapperObj * ) malloc ( sizeof(SnpUpdateOverflowWrapperObj) );
    obj[threadId]->snpBundle = snpBundle;
    obj[threadId]->vcSetting = vcSetting;
    obj[threadId]->threadId = threadId;
    pthread_create ( &thread[threadId], NULL, updateOverflowCounterFromBufferWrapper, ( void * ) obj[threadId] );
  }

  for ( threadId = 0; threadId < numOfCPUThreads; ++threadId )
  {
    pthread_join ( thread[threadId], NULL );
    free ( obj[threadId] );
  }
  free ( obj );
}

void updateDupSnpCounterForReads ( SnpBundle snpBundle, VcSetting * vcSetting, int trimHeadSize, int trimTailSize,
    unsigned char * recalScores,
    LongSoftClipPositionArray * dupLongSoftClipArray, unsigned int longSoftClipThreshold,
    unsigned int position, const char strand,
    unsigned char * query, char * qualities, unsigned int readLength, char * cigar,
    FILE * softclipFilePtr)
{
  

  unsigned int readPos = 0;
  unsigned int genomePosition = position;
  if ( cigar == NULL )
  {
    while ( readPos < readLength )
    {
      if ( readPos < trimHeadSize )
      {
        ++genomePosition;
        ++readPos;
        continue;
      }

      if ( readPos >= readLength - trimTailSize )
      {
        break;
      }

      char baseBits = strand ? query[readPos] : 3 - query[readLength - readPos - 1];
      char qScore = strand ? qualities[readPos] : qualities[readLength - readPos - 1];
      if ( qScore < 0 || qScore > 40 )
      {
        if ( vcSetting->enableQualityCorrection )
        {
          qScore = qScore < 0 ? 0 : 40;
        }
        else
        {
          fprintf ( stderr, "Error in Base Quality. Please check your input read file.\n" );
          exit ( 1 );
        }
      }
      char recalScore = recalScores[qScore];
      
      char weight = vcSetting->weightMap[recalScore];

      updateDupSnpCounter ( snpBundle, vcSetting,
          genomePosition, baseBits, -weight, strand, -1 );

      ++genomePosition;
      ++readPos;
    }
  }
  else
  {
    unsigned int i = 0;

    float totalWeight = 0.0f;
    unsigned char insertSeq[200];
    memset ( insertSeq, 0, 200 );
    unsigned char insertLength = 0;
    unsigned int insertPos = 0;

    unsigned char deleteLength = 0;

    int trimmedHead = 0;
    while ( readPos < readLength )
    {
      int range = 0;
      while ( cigar[i] <= 57 )
      {
        range = range * 10 + ( cigar[i] - 48 );
        ++i;
      }

      if ( cigar[i] == CIGAR_MATCHMISMATCH )
      {
        if ( readPos < trimHeadSize )
        {
          int trimRange = ( range < trimHeadSize - trimmedHead ? range : trimHeadSize - trimmedHead );
          genomePosition += trimRange;
          readPos += trimRange;
          range -= trimRange;
          trimmedHead += trimRange;
        }

        if ( readPos >= readLength - trimTailSize )
        {
          break;
        }

        if ( insertLength > 0 )
        {
          updateDupSnpCounterForInsert ( snpBundle, vcSetting,
              insertPos, insertLength,
              ( char ) ( totalWeight / insertLength ), insertSeq );
          totalWeight = 0.0f;
          memset ( insertSeq, 0, 200 );
          insertLength = 0;
        }
        if ( deleteLength > 0 )
        {
          char qScore = strand ? qualities[readPos] : qualities[readLength - readPos - 1];
          if ( qScore < 0 || qScore > 40 )
          {
            if ( vcSetting->enableQualityCorrection )
            {
              qScore = qScore < 0 ? 0 : 40;
            }
            else
            {
              fprintf ( stderr, "Error in Base Quality. Please check your input read file.\n" );
              exit ( 1 );
            }
          }
          char recalScore = recalScores[qScore];
          
          char weight = vcSetting->weightMap[recalScore];

          updateDupSnpCounterForDelete ( snpBundle, vcSetting,
              genomePosition - deleteLength, deleteLength, weight );

          deleteLength = 0;
        }

        while ( range > 0 )
        {
          if ( readPos >= readLength - trimTailSize )
          {
            readPos = readLength;
            break;
          }

          char baseBits = strand ? query[readPos] : 3 - query[readLength - readPos - 1];
          char qScore = strand ? qualities[readPos] : qualities[readLength - readPos - 1];
          if ( qScore < 0 || qScore > 40 )
          {
            if ( vcSetting->enableQualityCorrection )
            {
              qScore = qScore < 0 ? 0 : 40;
            }
            else
            {
              fprintf ( stderr, "Error in Base Quality. Please check your input read file.\n" );
              exit ( 1 );
            }
          }
          char recalScore = recalScores[qScore];
          
          char weight = vcSetting->weightMap[recalScore];

          updateDupSnpCounter ( snpBundle, vcSetting,
              genomePosition, baseBits, -weight, strand, -1 );

          genomePosition++;
          ++readPos;
          --range;
        }
      }
      else if ( cigar[i] == CIGAR_INSERT || cigar[i] == CIGAR_SOFT_CLIP )
      {
        if ( cigar[i] == CIGAR_SOFT_CLIP && range >= vcSetting->softClipThreshold )
        {
          unsigned int clipPos = readPos ? genomePosition - 1 : genomePosition;
          assert ( clipPos < vcSetting->dnaLength );
#ifdef ENABLE_SOFTCLIP_COUNTER
          unsigned char oldCount = 0;
          unsigned char newCount = 0;
          do
          {
            oldCount = snpCounter[clipPos].softClipCount;
            newCount = oldCount - 1;
            if ( newCount == 254 || oldCount == 0 )
            {
              break;
            }
          } while ( ! __sync_bool_compare_and_swap ( &( snpCounter[clipPos].softClipCount ), oldCount, newCount ) );
#endif
#ifdef SNP_DIRECTIONAL_SOFTCLIP
          SnpDirectionalSoftclip bufferSoftclip = { clipPos, ( readPos ? 0 : 1) };
          fwrite ( &bufferSoftclip, sizeof ( SnpDirectionalSoftclip ), 1, softclipFilePtr );
#endif
        }
        

        if ( readPos < trimHeadSize )
        {
          readPos += range;
          ++i;
          trimmedHead += range;
          continue;
        }

        if ( readPos + range >= readLength )
        {
          break;
        }

        if ( cigar[i] == CIGAR_INSERT )
        {
          while ( range > 0 )
          {
            char baseBits = strand ? query[readPos] : 3 - query[readLength - readPos - 1];
            char qScore = strand ? qualities[readPos] : qualities[readLength - readPos - 1];
            if ( qScore < 0 || qScore > 40 )
            {
              if ( vcSetting->enableQualityCorrection )
              {
                qScore = qScore < 0 ? 0 : 40;
              }
              else
              {
                fprintf ( stderr, "Error in Base Quality. Please check your input read file.\n" );
                exit ( 1 );
              }
            }
            char recalScore = recalScores[qScore];
            
            char weight = vcSetting->weightMap[recalScore];

            totalWeight += weight;
            assert ( ( insertLength >> 2 ) < 200 );
            insertSeq[insertLength >> 2] |= ( baseBits << ( ( insertLength & 3 ) << 1 ) );
            insertLength++;

            ++readPos;
            --range;
          }
          insertPos = genomePosition;
        }
        else
        {
          readPos += range;
        }

      }
      else if ( cigar[i] == CIGAR_DELETE )
      {
        if ( readPos < trimHeadSize )
        {
          genomePosition += range;
          ++i;
          continue;
        }

        
        if ( readPos + 1 >= readLength )
        {
          break;
        }

        deleteLength = range;


        genomePosition += range;
      }
      ++i;
    }
  }
}

void updateDupSnpCounter ( SnpBundle snpBundle, VcSetting * vcSetting,
    unsigned int position, char base, char weight,
    char strand, char strandCount )
{
  SnpCounter * snpCounter = snpBundle.snpCounter;
  SnpOverflowCounterArray * snpOverflowCounterArray = snpBundle.snpOverflowCounterArray;
  unsigned int numOfCPUThreads = snpBundle.numOfCPUThreads;
  unsigned int * invalidPos = snpBundle.invalidSnpCounterPos;

  assert ( position < vcSetting->dnaLength );

  if ( invalidPos[position >> 5] & ( 1 << ( 31 - ( position & 31 ) ) ) )
  {
    return;
  }

  unsigned short * weightedCounter;

  if ( ( ( snpCounter[position].weightedCount1 & snpCounter[position].weightedCount2 ) & 0xC000 ) == 0xC000 )
  {
    unsigned short * strandBigCounter;

    int region = findRegionByPosition ( position, vcSetting->dnaLength, numOfCPUThreads );
    assert ( snpCounter[position].data.arrayIndex < snpOverflowCounterArray[region].limit );
    weightedCounter = & ( snpOverflowCounterArray[region].counters[snpCounter[position].data.arrayIndex].weightedCount[base] );
    if ( strand )
    {
      strandBigCounter = & ( snpOverflowCounterArray[region].counters[snpCounter[position].data.arrayIndex].posStrandCount[base] );
    }
    else
    {
      strandBigCounter = & ( snpOverflowCounterArray[region].counters[snpCounter[position].data.arrayIndex].negStrandCount[base] );
    }

    if ( * strandBigCounter != 0 || strandCount == 0 )
    {
      
    }
#ifndef DISABLE_SNP_WARNING
    else
    {
      fprintf ( stderr, "Warning: Strand Counter (O) at %u for base %c is already 0 (w: %d, s:%d, %c) during duplicate removal.\n", position, dnaChar[base], weight, strandCount, strand ? '+' : '-' );
    }
#endif

    if ( * weightedCounter != 0 || weight == 0 )
    {
      __sync_add_and_fetch ( weightedCounter, weight );
    }
#ifndef DISABLE_SNP_WARNING
    else
    {
      fprintf ( stderr, "1 Warning: Weighted Counter (O) at %u for base %c is already 0 (w: %d, s:%d, %c) during duplicate removal.\n", position, dnaChar[base], weight, strandCount, strand ? '+' : '-'  );
    }
#endif
  }
  else
  {
    unsigned char * strandCounter;

    if ( ( ( snpCounter[position].weightedCount1 >> 14 ) & 3 ) == base )
    {
      weightedCounter = &( snpCounter[position].weightedCount1 );
      strandCounter = strand ? &( snpCounter[position].data.strandCounts.posStrandCount1 ) : &( snpCounter[position].data.strandCounts.negStrandCount1 );
    }
    else if ( ( ( snpCounter[position].weightedCount2 >> 14 ) & 3 ) == base )
    {
      weightedCounter = &( snpCounter[position].weightedCount2 );
      strandCounter = strand ? &( snpCounter[position].data.strandCounts.posStrandCount2 ) : &( snpCounter[position].data.strandCounts.negStrandCount2 );
    }
    else
    {
#ifndef DISABLE_SNP_WARNING
      fprintf ( stderr, "Removing a counter that does not exist at %u for base %c [%hu %hu]\n", position, dnaChar[base], snpCounter[position].weightedCount1, snpCounter[position].weightedCount2 );
#endif
      return;
    }

    if ( * strandCounter != 0 || strandCount == 0 )
    {
      __sync_add_and_fetch ( strandCounter, strandCount );
      
    }
#ifndef DISABLE_SNP_WARNING
    else
    {
      fprintf ( stderr, "Warning: Strand Counter at %u for base %c is already 0 (w: %d, s:%d, %c) during duplicate removal.\n", position, dnaChar[base], weight, strandCount, strand ? '+' : '-'  );
    }
#endif

    if ( ( ( ( * weightedCounter ) & 0x3FFF ) != 0  && ( * weightedCounter ) >= ( -weight ) ) || weight == 0 )
    {
      __sync_add_and_fetch ( weightedCounter, weight );
    }
#ifndef DISABLE_SNP_WARNING
    else
    {
      fprintf ( stderr, "2 Warning: Weighted Counter at %u for base %c is already 0 (w: %d, s:%d, %c) during duplicate removal.\n", position, dnaChar[base], weight, strandCount, strand ? '+' : '-'  );
    }
#endif
  }
}

void updateDupSnpCounterForInsert ( SnpBundle snpBundle, VcSetting * vcSetting,
    unsigned int position, unsigned char insertLength,
    unsigned char avgWeight, unsigned char * insertSeq )
{
  SnpCounter * snpCounter = snpBundle.snpCounter;
  SnpOverflowCounterArray * snpOverflowCounterArray = snpBundle.snpOverflowCounterArray;
  unsigned int numOfCPUThreads = snpBundle.numOfCPUThreads;
  unsigned int * invalidPos = snpBundle.invalidSnpCounterPos;
  MemoryPool * pool = snpBundle.snpMemoryPool;

  assert ( position < vcSetting->dnaLength );

  if ( invalidPos[position >> 5] & ( 1 << ( 31 - ( position & 31 ) ) ) )
  {
    return;
  }

  unsigned char hqCount = 0;
  unsigned char lqCount = 0;
  if ( avgWeight > vcSetting->indelWeightThreshold )
  {
    hqCount = 1;
  }
  else
  {
    lqCount = 1;
  }

  if ( ( ( snpCounter[position].weightedCount1 & snpCounter[position].weightedCount2 ) & 0xC000 ) != 0xC000 )
  {
#ifndef DISABLE_SNP_WARNING
    fprintf ( stderr, "Removing a non-existing Insertion ( 1 ) [%u]\n", position );
#endif
    return;
  }

  int region = findRegionByPosition ( position, vcSetting->dnaLength, numOfCPUThreads );
  unsigned int arrayIndex = snpCounter[position].data.arrayIndex;
  assert ( arrayIndex < snpOverflowCounterArray[region].limit );

  if ( snpCounter[position].weightedCount1 & 0x2000 )
    
  {
    if ( !( snpCounter[position].weightedCount1 & 0x1000 ) )
      
    {
      if ( ( ( snpCounter[position].weightedCount1 >> 8 ) & 15 ) == insertLength )
      {
        char sameSeq = 1;
        for ( int i = 0; i < ( ( insertLength + 3 ) / 4 ); ++i )
        {
          if ( snpOverflowCounterArray[region].counters[arrayIndex].insertion.counters.insertSeq[i] != insertSeq[i] )
          {
            sameSeq = 0;
            break;
          }
        }

        if ( sameSeq )
        {
          if ( hqCount )
          {
            if ( snpCounter[position].weightedCount1 & 0xFFF )
            {
              __sync_sub_and_fetch ( & ( snpCounter[position].weightedCount1 ), hqCount );
            }
#ifndef DISABLE_SNP_WARNING
            else
            {
              fprintf ( stderr, "Removing Insertion with zero count ( 1 ) [%u]\n", position );
            }
#endif
          }
          else if ( lqCount )
          {
            if ( snpOverflowCounterArray[region].counters[arrayIndex].insertion.counters.lqCount > 0 )
            {
              __sync_sub_and_fetch ( & ( snpOverflowCounterArray[region].counters[arrayIndex].insertion.counters.lqCount ), lqCount );
            }
#ifndef DISABLE_SNP_WARNING
            else
            {
              fprintf ( stderr, "Removing Insertion with zero count ( 2 ) [%u]\n", position );
            }
#endif
          }
          return;
        }
      }
#ifndef DISABLE_SNP_WARNING
      fprintf ( stderr, "Removing a non-existing Insertion ( 2 ) [%u]\n", position );
#endif
    }
    else
      
    {
      unsigned short arraySize = ( snpCounter[position].weightedCount1 & 0xFFF );
      unsigned char * insertInfo = ( unsigned char * ) getAddress ( pool, snpOverflowCounterArray[region].counters[arrayIndex].insertion.ptr );
      int infoIndex = 0;
      while ( infoIndex < arraySize )
      {
        unsigned char size = insertInfo[infoIndex++];
        if ( insertLength == size )
        {
          char sameSeq = 1;
          for ( int i = 0; i < ( size + 3 ) / 4; ++i )
          {
            if ( insertInfo[infoIndex + i] != insertSeq[i] )
            {
              sameSeq = 0;
              break;
            }
          }
          if ( sameSeq )
          {
            infoIndex += ( size + 3 ) / 4;
            unsigned short mphqCount;
            mphqCount = insertInfo[infoIndex] & 0xff;
            mphqCount |= insertInfo[infoIndex+1] << 8;

            if ( mphqCount > 0 || hqCount == 0 )
            {
              
              unsigned char _old;
              unsigned char _new;
              do {
                _old = insertInfo[infoIndex];
                _new = _old - hqCount;
              } while ( ! __sync_bool_compare_and_swap ( & ( insertInfo[infoIndex] ), _old, _new ) );
              if ( _old < hqCount )
              {
                --insertInfo[infoIndex+1];
              }
            }
#ifndef DISABLE_SNP_WARNING
            else
            {
              fprintf ( stderr, "Removing Insertion with zero count ( 3 ) [%u]\n", position );
            }
#endif
            infoIndex+=2;

            unsigned short mplqCount;
            mplqCount = insertInfo[infoIndex] & 0xff;
            mplqCount |= insertInfo[infoIndex+1] & 0xff;

            if ( mplqCount > 0 || lqCount == 0 )
            {
              
              unsigned char _old;
              unsigned char _new;
              do {
                _old = insertInfo[infoIndex];
                _new = _old - lqCount;
              } while ( ! __sync_bool_compare_and_swap ( & ( insertInfo[infoIndex] ), _old, _new ) );
              if ( _old < lqCount )
              {
                --insertInfo[infoIndex+1];
              }

            }
#ifndef DISABLE_SNP_WARNING
            else
            {
              fprintf ( stderr, "Removing Insertion with zero count ( 4 ) [%u]\n", position );
            }
#endif
            infoIndex+=2;
            return;
          }
        }
        infoIndex += ( size + 3 ) / 4 + 4;
      }
#ifndef DISABLE_SNP_WARNING
      fprintf ( stderr, "Removing a non-existing Insertion ( 3 ) [%u]\n", position );
      const char dnaMap[] = {'A','C','G','T'};
      for ( int i=0;i<insertLength;++i )
      {
        char baseBit = ( insertSeq[i >> 2] >> ( ( i & 3 ) << 1 ) ) & 3;
        fprintf(stderr, "%c", dnaMap[baseBit]);
      }
      fprintf ( stderr, "\n" );
#endif
    }
  }
#ifndef DISABLE_SNP_WARNING
  else
    
  {
    fprintf ( stderr, "Removing a non-existing Insertion ( 4 ) [%u]\n", position );
  }
#endif
}

void updateDupSnpCounterForDelete ( SnpBundle snpBundle, VcSetting * vcSetting,
    unsigned int position, unsigned char deleteLength, char weight )
{
  SnpCounter * snpCounter = snpBundle.snpCounter;
  SnpOverflowCounterArray * snpOverflowCounterArray = snpBundle.snpOverflowCounterArray;
  unsigned int numOfCPUThreads = snpBundle.numOfCPUThreads;
  unsigned int * invalidPos = snpBundle.invalidSnpCounterPos;
  MemoryPool * pool = snpBundle.snpMemoryPool;

  assert ( position < vcSetting->dnaLength );

  if ( invalidPos[position >> 5] & ( 1 << ( 31 - ( position & 31 ) ) ) )
  {
    return;
  }

  unsigned char hqCount = 0;
  unsigned char lqCount = 0;
  if ( weight > vcSetting->indelWeightThreshold )
  {
    hqCount = 1;
  }
  else
  {
    lqCount = 1;
  }

  if ( ( ( snpCounter[position].weightedCount1 & snpCounter[position].weightedCount2 ) & 0xC000 ) != 0xC000 )
  {
#ifndef DISABLE_SNP_WARNING
    fprintf ( stderr, "Removing a non-existing Deletion ( 1 )\n" );
#endif
    return;
  }

  int region = findRegionByPosition ( position, vcSetting->dnaLength, numOfCPUThreads );
  unsigned int arrayIndex = snpCounter[position].data.arrayIndex;
  assert ( arrayIndex < snpOverflowCounterArray[region].limit );

  if ( snpCounter[position].weightedCount2 & 0x2000 )
  {
    if ( !( snpCounter[position].weightedCount2 & 0x1000 ) )
    {
      if ( snpOverflowCounterArray[region].counters[arrayIndex].deletion.counters.deletionLength == deleteLength )
      {
        if ( hqCount )
        {
          if ( snpOverflowCounterArray[region].counters[arrayIndex].deletion.counters.hqCount > 0 )
          {
            __sync_sub_and_fetch ( & ( snpOverflowCounterArray[region].counters[arrayIndex].deletion.counters.hqCount ), hqCount );
          }
#ifndef DISABLE_SNP_WARNING
          else
          {
            fprintf ( stderr, "Removing Deletion with zero count ( 1 )\n" );
          }
#endif
        }
        else if ( lqCount )
        {
          if ( snpOverflowCounterArray[region].counters[arrayIndex].deletion.counters.lqCount > 0 )
          {
            __sync_sub_and_fetch ( & ( snpOverflowCounterArray[region].counters[arrayIndex].deletion.counters.lqCount ), lqCount );
          }
#ifndef DISABLE_SNP_WARNING
          else
          {
            fprintf ( stderr, "Removing Deletion with zero count ( 2 ) %u\n", position );
          }
#endif
        }
      }
#ifndef DISABLE_SNP_WARNING
      else
      {
        fprintf ( stderr, "Removing a non-existing Deletion ( 3 )\n" );
      }
#endif
    }
    else
    {
      unsigned short arraySize = ( snpCounter[position].weightedCount2 & 0xFFF );
      unsigned char * deleteInfo = ( unsigned char * ) getAddress ( pool, snpOverflowCounterArray[region].counters[arrayIndex].deletion.ptr );
      int infoIndex = 0;
      while ( infoIndex < arraySize )
      {
        unsigned char length = deleteInfo[infoIndex++];
        if ( length == deleteLength )
        {
          unsigned short mphqCount;
          mphqCount = deleteInfo[infoIndex] & 0xff;
          mphqCount |= deleteInfo[infoIndex+1] << 8;
          if ( mphqCount > 0 || hqCount == 0 )
          {
            
            unsigned char _old;
            unsigned char _new;
            do {
              _old = deleteInfo[infoIndex];
              _new = _old - hqCount;
            } while ( ! __sync_bool_compare_and_swap ( & ( deleteInfo[infoIndex] ), _old, _new ) );
            if ( _old < hqCount )
            {
              --deleteInfo[infoIndex+1];
            }
          }
#ifndef DISABLE_SNP_WARNING
          else
          {
            fprintf ( stderr, "Removing Deletion with zero count ( 3 )\n" );
          }
#endif
          infoIndex+=2;
          unsigned short mplqCount;
          mplqCount = deleteInfo[infoIndex] & 0xff;
          mplqCount |= deleteInfo[infoIndex+1] << 8;
          if ( mplqCount > 0 || lqCount == 0 )
          {
            
            unsigned char _old;
            unsigned char _new;
            do {
              _old = deleteInfo[infoIndex];
              _new = _old - lqCount;
            } while ( ! __sync_bool_compare_and_swap ( & ( deleteInfo[infoIndex] ), _old, _new ) );
            if ( _old < lqCount )
            {
              --deleteInfo[infoIndex+1];
            }
          }
#ifndef DISABLE_SNP_WARNING
          else
          {
            fprintf ( stderr, "Removing Deletion with zero count ( 4 ) %u\n", position );
          }
#endif
          infoIndex+=2;
          return;
        }
        else
        {
          infoIndex += 4;
        }
      }
#ifndef DISABLE_SNP_WARNING
      fprintf ( stderr, "Removing a non-existing Deletion ( 4 )\n" );
#endif
    }
  }
#ifndef DISABLE_SNP_WARNING
  else
  {
    fprintf ( stderr, "Removing a non-existing Deletion ( 2 )\n" );
  }
#endif
}

void reupdateSnpCounter ( SnpBundle snpBundle, VcSetting * vcSetting, int trimHeadSize, int trimTailSize,
    unsigned char * query, uint readLength, char * qualities,
    unsigned int position, const char strand,
    unsigned char * cigar, unsigned char * scores,
    unsigned int ra_startPos, unsigned int ra_endPos,
    FILE * softclipFilePtr )
{
  unsigned int numOfCPUThreads = snpBundle.numOfCPUThreads;

  static char stateMap[] = { 'M', 'I', 'D', 'S' };
  unsigned int genomePos = position;
  unsigned int readPos = 0;
  unsigned int i = 0;
  int trimmedHead = 0;

  while ( readPos < readLength )
  {
    if ( i >= RA_CIGAR_LENGTH )
    {
      printf ( "Reupdate SNP counter error: cigar length too long\n" );
      exit ( 1 );
    }

    unsigned int range = ( cigar[i] >> 2 ) & 0x3F;
    char state = stateMap[cigar[i] & 3];

    while ( i + 1 < RA_CIGAR_LENGTH && state == stateMap[cigar[i + 1] & 3] && state != 'M' )
    {
      i++;
      range += ( cigar[i] >> 2 ) & 0x3F;
    }

    if ( range == 0 )
    {
      i++;
      continue;
    }

    if ( state == CIGAR_MATCHMISMATCH )
    {
      if ( readPos < trimHeadSize )
      {
        int trimRange = ( range < trimHeadSize - trimmedHead ? range : trimHeadSize - trimmedHead );
        genomePos += trimRange;
        readPos += trimRange;
        range -= trimRange;
        trimmedHead += trimRange;
      }

      while ( range > 0 )
      {
        if ( readPos >= readLength - trimTailSize )
        {
          readPos = readLength;
          break;
        }

        if ( ra_startPos <= genomePos && genomePos <= ra_endPos )
        {
          char baseBit = strand ? query[readPos] : 3 - query[readLength - readPos - 1];
          char qScore = strand ? qualities[readPos] : qualities[readLength - readPos - 1];
          if ( qScore < 0 || qScore > 40 )
          {
            if ( vcSetting->enableQualityCorrection )
            {
              qScore = qScore < 0 ? 0 : 40;
            }
            else
            {
              fprintf ( stderr, "Error in Base Quality. Please check your input read file.\n" );
              exit ( 1 );
            }
          }
          char recalScore = scores[qScore];
          
          char weight = vcSetting->weightMap[recalScore];
          unsigned int threadId = findRegionByPosition ( genomePos, vcSetting->dnaLength, numOfCPUThreads );

          updateSnpCounterForOverflow ( snpBundle, vcSetting,
              genomePos, baseBit, weight,
              strand, 1, threadId );
        }

        readPos++;
        genomePos++;
        range--;
      }
    }
    else if ( state == CIGAR_INSERT || state == CIGAR_SOFT_CLIP )
    {
      if ( state == CIGAR_SOFT_CLIP && range >= 5 )
      {
        unsigned int clipPos = readPos ? genomePos - 1 : genomePos;
#ifdef ENABLE_SOFTCLIP_COUNTER
        if ( snpCounter[clipPos].softClipCount < 255 )
        {
          snpCounter[clipPos].softClipCount++;
        }
#endif
#ifdef SNP_DIRECTIONAL_SOFTCLIP
        SnpDirectionalSoftclip bufferSoftclip = { clipPos, ( readPos ? 0 : 1) };
        fwrite ( &bufferSoftclip, sizeof ( SnpDirectionalSoftclip ), 1, softclipFilePtr );
#endif
      }

      if ( readPos < trimHeadSize )
      {
        readPos += range;
        i++;
        trimmedHead += range;
        continue;
      }
      if ( readPos + range >= readLength )
      {
        break;
      }

      if ( state == CIGAR_INSERT )
      {
        unsigned char insertLength = range;
        unsigned char insertSeq[MAX_READ_LENGTH] = { 0 };
        unsigned char avgWeight = 0;
        unsigned int threadId = findRegionByPosition ( genomePos, vcSetting->dnaLength, numOfCPUThreads );
        char baseBit;
        char recalScore;

        for ( int k = 0; k < range; k++ )
        {
          baseBit = strand ? query[readPos] : 3 - query[readLength - readPos - 1];
          char qScore = strand ? qualities[readPos] : qualities[readLength - readPos - 1];
          if ( qScore < 0 || qScore > 40 )
          {
            if ( vcSetting->enableQualityCorrection )
            {
              qScore = qScore < 0 ? 0 : 40;
            }
            else
            {
              fprintf ( stderr, "Error in Base Quality. Please check your input read file.\n" );
              exit ( 1 );
            }
          }
          recalScore = scores[qScore];
          
          avgWeight += vcSetting->weightMap[recalScore];
          insertSeq[k >> 2] |= ( baseBit << ( ( k & 3 ) << 1 ) );
          readPos++;
        }
        avgWeight /= range;

        if ( ra_startPos <= genomePos && genomePos <= ra_endPos )
        {
          updateSnpCounterForInsert ( snpBundle, vcSetting, genomePos, insertLength,
              avgWeight, insertSeq, threadId );
        }

        range = 0;
      }
      else
      {
        readPos += range;
      }

    }
    else
    {
      if ( readPos < trimHeadSize )
      {
        genomePos += range;
        i++;
        continue;
      }

      if ( readPos + 1 >= readLength )
      {
        break;
      }
      unsigned char deleteLength = range;
      char qScore = strand ? qualities[readPos] : qualities[readLength - readPos - 1];
      if ( qScore < 0 || qScore > 40 )
      {
        if ( vcSetting->enableQualityCorrection )
        {
          qScore = qScore < 0 ? 0 : 40;
        }
        else
        {
          fprintf ( stderr, "Error in Base Quality. Please check your input read file.\n" );
          exit ( 1 );
        }
      }
      char recalScore = scores[qScore];
      
      char weight = vcSetting->weightMap[recalScore];
      unsigned int threadId = findRegionByPosition ( genomePos, vcSetting->dnaLength, numOfCPUThreads );

      if ( ra_startPos <= genomePos && genomePos <= ra_endPos )
      {
        updateSnpCounterForDelete ( snpBundle, vcSetting, genomePos,
            deleteLength, weight, threadId );
      }

      genomePos += range;

    }
    i++;
  }
}

void reupdateDupSnpCounter ( SnpBundle snpBundle, VcSetting * vcSetting, int trimHeadSize, int trimTailSize,
    unsigned char * query, uint readLength, char * qualities,
    unsigned int position, const char strand,
    unsigned char * cigar, unsigned char * scores,
    FILE * softclipFilePtr )
{

  static char stateMap[] = {'M', 'I', 'D', 'S'};
  unsigned int genomePos = position;
  unsigned int readPos = 0;
  unsigned int i = 0;
  int trimmedHead = 0;

  float totalWeight = 0.0f;
  unsigned char insertSeq[MAX_READ_LENGTH];
  memset ( insertSeq, 0, MAX_READ_LENGTH );
  unsigned char insertLength = 0;
  unsigned int insertPos = 0;
  unsigned char deleteLength = 0;

  while ( readPos < readLength )
  {
    if ( i >= RA_CIGAR_LENGTH )
    {
      printf ( "reupdateDupSnpCounter error: cigar length too long\n" );
      exit ( 1 );
    }

    unsigned int range = ( cigar[i] >> 2 ) & 0x3F;
    char state = stateMap[cigar[i] & 3];

    if ( range == 0 )
    {
      i++;
      continue;
    }

    while ( i + 1 < RA_CIGAR_LENGTH && cigar[i + 1] != 0 && state == stateMap[cigar[i + 1] & 3] )
    {
      range += ( ( cigar[i + 1] >> 2 ) & 0x3F );
      i++;
    }

    if ( state == CIGAR_MATCHMISMATCH )
    {
      if ( readPos < trimHeadSize )
      {
        int trimRange = ( range < trimHeadSize - trimmedHead ? range : trimHeadSize - trimmedHead );
        genomePos += trimRange;
        readPos += trimRange;
        range -= trimRange;
        trimmedHead += trimRange;
      }

      if ( readPos >= readLength - trimTailSize )
      {
        break;
      }

      if ( insertLength > 0 )
      {
        updateDupSnpCounterForInsert ( snpBundle, vcSetting,
            insertPos, insertLength,
            ( char ) ( totalWeight / insertLength ), insertSeq );

        totalWeight = 0.0f;
        memset ( insertSeq, 0, MAX_READ_LENGTH );
        insertLength = 0;
      }

      if ( deleteLength > 0 )
      {
        char qScore = strand ? qualities[readPos] : qualities[readLength - readPos - 1];
        if ( qScore < 0 || qScore > 40 )
        {
          if ( vcSetting->enableQualityCorrection )
          {
            qScore = qScore < 0 ? 0 : 40;
          }
          else
          {
            fprintf ( stderr, "Error in Base Quality. Please check your input read file.\n" );
            exit ( 1 );
          }
        }
        char recalScore = scores[qScore];
        
        char weight = vcSetting->weightMap[recalScore];

        updateDupSnpCounterForDelete ( snpBundle, vcSetting,
            genomePos - deleteLength, deleteLength, weight );

        deleteLength = 0;
      }

      while ( range > 0 )
      {
        if ( readPos >= readLength - trimTailSize )
        {
          readPos = readLength;
          break;
        }

        char baseBit = strand ? query[readPos] : 3 - query[readLength - readPos - 1];
        char qScore = strand ? qualities[readPos] : qualities[readLength - readPos - 1];
        if ( qScore < 0 || qScore > 40 )
        {
          if ( vcSetting->enableQualityCorrection )
          {
            qScore = qScore < 0 ? 0 : 40;
          }
          else
          {
            fprintf ( stderr, "Error in Base Quality. Please check your input read file.\n" );
            exit ( 1 );
          }
        }
        char recalScore = scores[qScore];
        
        char weight = vcSetting->weightMap[recalScore];

        updateDupSnpCounter ( snpBundle, vcSetting,
            genomePos, baseBit, -weight, strand, -1 );

        readPos++;
        genomePos++;
        range--;
      }
    }
    else if ( state == CIGAR_INSERT || state == CIGAR_SOFT_CLIP )
    {
      if ( state == CIGAR_SOFT_CLIP && range >= 5 )
      {
        unsigned int clipPos = readPos ? genomePos - 1 : genomePos;
#ifdef ENABLE_SOFTCLIP_COUNTER
        if ( !( snpCounter[clipPos].softClipCount == 255
              || snpCounter[clipPos].softClipCount == 0 ) )
        {
          snpCounter[clipPos].softClipCount--;
        }
#endif
#ifdef SNP_DIRECTIONAL_SOFTCLIP
        SnpDirectionalSoftclip bufferSoftclip = { clipPos, ( readPos ? 0 : 1) };
        fwrite ( &bufferSoftclip, sizeof ( SnpDirectionalSoftclip ), 1, softclipFilePtr );
#endif
      }

      if ( readPos < trimHeadSize )
      {
        readPos += range;
        i++;
        continue;
      }
      if ( readPos + range >= readLength )
      {
        break;
      }

      if ( state == CIGAR_INSERT )
      {
        char baseBit;
        char recalScore;

        for ( int k = 0; k < range; k++ )
        {
          baseBit = strand ? query[readPos] : 3 - query[readLength - readPos - 1];
          char qScore = strand ? qualities[readPos] : qualities[readLength - readPos - 1];
          if ( qScore < 0 || qScore > 40 )
          {
            if ( vcSetting->enableQualityCorrection )
            {
              qScore = qScore < 0 ? 0 : 40;
            }
            else
            {
              fprintf ( stderr, "Error in Base Quality. Please check your input read file.\n" );
              exit ( 1 );
            }
          }
          recalScore = scores[qScore];
          totalWeight += vcSetting->weightMap[recalScore];
          insertSeq[insertLength + ( k >> 2 )] |= ( baseBit << ( ( k & 3 ) << 1 ) );
          readPos++;
        }
        insertLength += range;
        insertPos = genomePos;
        range = 0;
      }
      else
      {
        if ( insertLength > 0 )
        {
          updateDupSnpCounterForInsert ( snpBundle, vcSetting,
              insertPos, insertLength,
              ( char ) ( totalWeight / insertLength ), insertSeq );

          totalWeight = 0.0f;
          memset ( insertSeq, 0, MAX_READ_LENGTH );
          insertLength = 0;
        }

        if ( deleteLength > 0 )
        {
          char qScore = strand ? qualities[readPos] : qualities[readLength - readPos - 1];
          if ( qScore < 0 || qScore > 40 )
          {
            if ( vcSetting->enableQualityCorrection )
            {
              qScore = qScore < 0 ? 0 : 40;
            }
            else
            {
              fprintf ( stderr, "Error in Base Quality. Please check your input read file.\n" );
              exit ( 1 );
            }
          }
          char recalScore = scores[qScore];
          
          char weight = vcSetting->weightMap[recalScore];

          updateDupSnpCounterForDelete ( snpBundle, vcSetting,
              genomePos - deleteLength, deleteLength, weight );

          deleteLength = 0;
        }

        readPos += range;
      }

    }
    else
    {
      if ( readPos < trimHeadSize )
      {
        genomePos += range;
        i++;
        continue;
      }

      if ( readPos + 1 >= readLength )
      {
        break;
      }

      deleteLength += range;

      genomePos += range;
    }
    i++;
  }
}




