#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "readIndex.h"


size_t loadTranslateWithTranslateSize ( const char * inputFileName, unsigned int & dnaLength,
                                     unsigned int ** ambiguityMap, Translate ** translate )
{
    FILE * inputFile;
    inputFile = ( FILE * ) fopen ( inputFileName, "r" );

    if ( inputFile == NULL )
    {
        fprintf ( stderr, "Cannot open file %s\n", inputFileName );
        exit ( 1 );
    }

    unsigned int i;
    unsigned int gridEntries;
    unsigned int removedSegmentCount;
    fscanf ( inputFile, "%u %d %d %u\n", &dnaLength, &i, &removedSegmentCount, &gridEntries );

    unsigned int * ambMap = ( unsigned int * ) malloc ( ( gridEntries ) * sizeof ( unsigned int ) );
    Translate * tran = ( Translate * ) malloc ( ( i + removedSegmentCount + 1 ) * sizeof ( Translate ) );

    unsigned int j = 0;

    while ( !feof ( inputFile ) && j < gridEntries )
    {
        fscanf ( inputFile, "%u\n", & ( ambMap[j] ) );
        j++;
    }

    if ( j < gridEntries )
    {
        fprintf ( stderr, "Translate missing entries!\n" );
        exit ( 1 );
    }

    j = 0;

    while ( !feof ( inputFile ) && j < i + removedSegmentCount )
    {
        fscanf ( inputFile, "%u %hu %u\n", & ( tran[j].startPos ), & ( tran[j].chrID ), & ( tran[j].correction ) );
        j++;
    }

    if ( j < i + removedSegmentCount )
    {
        fprintf ( stderr, "Translate missing entries\n" );
        exit ( 1 );
    }
    
    j = i + removedSegmentCount;
    tran[j].chrID = USHRT_MAX;
    tran[j].startPos = UINT_MAX;
    tran[j].correction = UINT_MAX;

    fclose ( inputFile );

    *ambiguityMap = ambMap;
    *translate = tran;

    return i + removedSegmentCount + 1; 
}

void loadTranslate( const char * inputFileName, unsigned int & dnaLength,
                    unsigned int ** ambiguityMap, Translate ** translate )
{
    loadTranslateWithTranslateSize( inputFileName, dnaLength, ambiguityMap, translate );
}

void freeTranslate ( unsigned int * ambiguityMap, Translate * translate )
{
    free ( ambiguityMap );
    free ( translate );
}

void loadSeqInfo ( const char * inputFileName, unsigned int & dnaLength,
                   Annotation ** annotation, SeqOffset ** seqOffset, unsigned int & numOfSeq )
{
    FILE * inputFile;
    inputFile = ( FILE * ) fopen ( inputFileName, "r" );

    if ( inputFile == NULL )
    {
        fprintf ( stderr, "Cannot open file %s\n", inputFileName );
        exit ( 1 );
    }

    unsigned int randomSeed;

    fscanf ( inputFile, "%u %d %u\n", &dnaLength, &numOfSeq, &randomSeed );

    if ( numOfSeq == 0 )
    {
        fprintf ( stderr, "Annotation empty entry\n" );
        exit ( 1 );
    }

    Annotation * ann = ( Annotation * ) malloc ( ( numOfSeq + 1 ) * sizeof ( Annotation ) );
    memset ( ann, 0, ( numOfSeq + 1 ) * sizeof ( Annotation ) );
    SeqOffset * seq = ( SeqOffset * ) malloc ( ( numOfSeq + 1 ) * sizeof ( SeqOffset ) );
    memset ( seq, 0, ( numOfSeq + 1 ) * sizeof ( SeqOffset ) );

    int i = 0;
    int j, k;
    while ( !feof ( inputFile ) && i < numOfSeq )
    {
        fscanf ( inputFile, "%u ", &( ann[i].gi ) );
        fgets ( ann[i].text, MAX_SEQ_NAME_LENGTH, inputFile );
        fscanf ( inputFile, "%u %u %d\n", &( seq[i].startPos ), &( seq[i].endPos ), &( seq[i + 1].firstAmbiguityIndex ) );
        seq[i].lastAmbiguityIndex = seq[i + 1].firstAmbiguityIndex;
        seq[i].endPos = seq[i].startPos + seq[i].endPos - 1;

        j = 0;
        while ( j < MAX_SEQ_NAME_LENGTH )
        {
            if ( ann[i].text[j] == '\n'
                 || ann[i].text[j] == '\t'
                 || ann[i].text[j] == ' '
                 || ann[i].text[j] == '\0' )
            {
                break;
            }
            j++;
        }
        ann[i].text[j] = '\0';

        j++;
        if ( ann[i].text[j] != '\0' )
        {
            k = 0;
            while ( j < MAX_SEQ_NAME_LENGTH )
            {
                ann[i].decoratedText[k] = ann[i].text[j];
                k++;
                j++;
            }
        }

        i++;
    }

    if ( i < numOfSeq )
    {
        fprintf ( stderr, "Annotation missing entries\n" );
        exit ( 1 );
    }
    ann[i].gi = 0;
    seq[0].firstAmbiguityIndex = 1;

    if ( annotation )
    {
        *annotation = ann;
    }
    else
    {
        free ( ann );
    }
    if ( seqOffset )
    {
        *seqOffset = seq;
    }
    else
    {
        free ( seq );
    }
}

void freeSeqInfo ( Annotation * annotation, SeqOffset * seqOffset )
{
    free ( annotation );
    free ( seqOffset );
}
