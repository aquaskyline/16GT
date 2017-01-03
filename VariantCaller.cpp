#include "VariantCaller.h"
#include "FisherExactTest.h"
#include "SNP_Caller.h"
#include "SNPFunctions.h"

#define chkPos -1

char getSnpMetrics(unsigned int *packedSeq,
                   SnpBundle snpBundle, unsigned int position, unsigned int dnaLength,
                   InputStat &iStat, unsigned int &weightedCount,
                   unsigned int &posStrandCount, unsigned int &negStrandCount,
                   unsigned char &gcCount, unsigned char &gCount,
                   unsigned int &indelHqCount, unsigned int &indelLqCount,
                   IndelInfo *indelInfo) {
  memset(&iStat, 0, sizeof(InputStat));
  weightedCount = 0;
  posStrandCount = 0;
  negStrandCount = 0;
  gcCount = 0;
  gCount = 0;
  indelHqCount = 0;
  indelLqCount = 0;

  unsigned int *invalidPos = snpBundle.invalidSnpCounterPos;

  if (position == chkPos) {
    fprintf(stderr, "[%s-%d] %u, %u, %u\n%u, %u, %u\n%u, %u, %u\n%u, %u, %u\n%u, %u, %u\n%u, %u, %u\n%u, %u, %u\n\n",
            __FILE__, __LINE__,
            position, indelHqCount, indelLqCount,
            iStat.W[0], iStat.F[0], iStat.R[0],
            iStat.W[1], iStat.F[1], iStat.R[1],
            iStat.W[2], iStat.F[2], iStat.R[2],
            iStat.W[3], iStat.F[3], iStat.R[3],
            iStat.W[4], iStat.F[4], iStat.R[4],
            iStat.W[5], iStat.F[5], iStat.R[5]
    );
  }

  if ((invalidPos[position >> 5] & (1 << (31 - (position & 31)))) || position >= dnaLength) {
    return 0;
  }

  SnpCounter *snpCounter = snpBundle.snpCounter;
  SnpOverflowCounterArray *snpOverflowCounterArray = snpBundle.snpOverflowCounterArray;
  unsigned int numOfCPUThreads = snpBundle.numOfCPUThreads;
  MemoryPool *pool = snpBundle.snpMemoryPool;
  unsigned int region = findRegionByPosition(position, dnaLength, numOfCPUThreads);

  if (((snpCounter[position].weightedCount1 & snpCounter[position].weightedCount2) & 0xC000) != 0xC000) {
    if (snpCounter[position].data.strandCounts.posStrandCount1 |
        snpCounter[position].data.strandCounts.negStrandCount1) {
      unsigned char base = (snpCounter[position].weightedCount1 >> 14) & 3;
      iStat.W[base] = snpCounter[position].weightedCount1 & 0x3FFF;
      iStat.F[base] = snpCounter[position].data.strandCounts.posStrandCount1;
      iStat.R[base] = snpCounter[position].data.strandCounts.negStrandCount1;

      weightedCount += iStat.W[base];
      posStrandCount += iStat.F[base];
      negStrandCount += iStat.R[base];
      if (base == 1 || base == 2) 
      {
        if (base == 2) {
          gCount = 1;
        }
        gcCount = 1;
      }
    }

    if (snpCounter[position].data.strandCounts.posStrandCount2 |
        snpCounter[position].data.strandCounts.negStrandCount2) {
      unsigned char base = (snpCounter[position].weightedCount2 >> 14) & 3;
      iStat.W[base] = snpCounter[position].weightedCount2 & 0x3FFF;
      iStat.F[base] = snpCounter[position].data.strandCounts.posStrandCount2;
      iStat.R[base] = snpCounter[position].data.strandCounts.negStrandCount2;

      weightedCount += iStat.W[base];
      posStrandCount += iStat.F[base];
      negStrandCount += iStat.R[base];
      if (base == 1 || base == 2) 
      {
        if (base == 2) {
          gCount = 1;
        }
        gcCount = 1;
      }
    }

    if (position == chkPos) {
      fprintf(stderr, "[%s-%d] %u, %u, %u\n%u, %u, %u\n%u, %u, %u\n%u, %u, %u\n%u, %u, %u\n%u, %u, %u\n%u, %u, %u\n\n",
              __FILE__, __LINE__,
              position, indelHqCount, indelLqCount,
              iStat.W[0], iStat.F[0], iStat.R[0],
              iStat.W[1], iStat.F[1], iStat.R[1],
              iStat.W[2], iStat.F[2], iStat.R[2],
              iStat.W[3], iStat.F[3], iStat.R[3],
              iStat.W[4], iStat.F[4], iStat.R[4],
              iStat.W[5], iStat.F[5], iStat.R[5]
      );
    }
  } else {
    unsigned int arrayIndex = snpCounter[position].data.arrayIndex;

    unsigned int i;
    for (i = 0; i < ALPHABET_SIZE; ++i) {
      if (snpOverflowCounterArray[region].counters[arrayIndex].posStrandCount[i] |
          snpOverflowCounterArray[region].counters[arrayIndex].negStrandCount[i]) {
        iStat.W[i] = snpOverflowCounterArray[region].counters[arrayIndex].weightedCount[i];
        iStat.F[i] = snpOverflowCounterArray[region].counters[arrayIndex].posStrandCount[i];
        iStat.R[i] = snpOverflowCounterArray[region].counters[arrayIndex].negStrandCount[i];

        weightedCount += iStat.W[i];
        posStrandCount += iStat.F[i];
        negStrandCount += iStat.R[i];
        if (i == 1 || i == 2) 
        {
          if (i == 2) {
            gCount = 1;
          }
          gcCount = 1;
        }
      }
    }
  }
#ifdef GENOTYPE_16
  unsigned int positionInsertion = position + 1;
  positionInsertion = positionInsertion >= dnaLength ? dnaLength - 1 : positionInsertion;
  unsigned int positionDeletion = position + 1;
  positionDeletion = positionDeletion >= dnaLength ? dnaLength - 1 : positionDeletion;
  if (!(((snpCounter[positionInsertion].weightedCount1 & snpCounter[positionInsertion].weightedCount2) & 0xC000) !=
        0xC000)) {
    unsigned int insertionIndex = snpCounter[positionInsertion].data.arrayIndex;
    unsigned int deletionIndex = snpCounter[positionDeletion].data.arrayIndex;

    unsigned int deepestIndel[3] = {0}; 
    unsigned int secondDeepestIndel[3] = {0}; 
    unsigned int _indelInfoCurrIndex = 0;
    unsigned int _indelInfoIndex0 = 0;
    unsigned int _indelInfoIndex1 = 0;
    unsigned int _indelInfoLength0 = 0;
    unsigned int _indelInfoLength1 = 0;

    if (indelInfo) {
      _indelInfoCurrIndex = indelInfo->size;
    }

    if (snpCounter[positionInsertion].weightedCount1 & 0x2000) {
      if (snpCounter[positionInsertion].weightedCount1 & 0x1000) {
        unsigned short arraySize = (snpCounter[positionInsertion].weightedCount1 & 0xFFF);
        unsigned char *insertInfo =
            (unsigned char *) getAddress(pool, snpOverflowCounterArray[region].counters[insertionIndex].insertion.ptr);
        unsigned int infoIndex = 0;
        unsigned int _index;
        while (infoIndex < arraySize) {
          unsigned char insertLength = insertInfo[infoIndex++];
          _index = infoIndex;
          infoIndex += (insertLength + 3) / 4;
          unsigned short hqCount;
          hqCount = insertInfo[infoIndex++];
          hqCount |= insertInfo[infoIndex++] << 8;
          unsigned short lqCount;
          lqCount = insertInfo[infoIndex++];
          lqCount |= insertInfo[infoIndex++] << 8;

          indelHqCount += hqCount;
          indelLqCount += lqCount;
          if (hqCount + lqCount > deepestIndel[1] + deepestIndel[2]) {
            if (indelInfo) {
              _indelInfoIndex1 = _indelInfoIndex0;
              _indelInfoLength1 = _indelInfoLength0;
              _indelInfoIndex0 = indelInfo->size;
              checkIndelInfoSpace(indelInfo, (4 + (insertLength + 3) / 4));
              indelInfo->array[indelInfo->size++] = 'I';
              indelInfo->array[indelInfo->size++] = insertLength;
              for (unsigned int i = _index; i < _index + (insertLength + 3) / 4; ++i) {
                indelInfo->array[indelInfo->size++] = insertInfo[i];
              }
              indelInfo->array[indelInfo->size++] = hqCount & 0xff;
              indelInfo->array[indelInfo->size++] = (hqCount >> 8) & 0xff;
              indelInfo->array[indelInfo->size++] = lqCount & 0xff;
              indelInfo->array[indelInfo->size++] = (lqCount >> 8) & 0xff;
              _indelInfoLength0 = indelInfo->size - _indelInfoIndex0;
            }

            secondDeepestIndel[0] = deepestIndel[0];
            secondDeepestIndel[1] = deepestIndel[1];
            secondDeepestIndel[2] = deepestIndel[2];
            deepestIndel[0] = hqCount * 4; 
            deepestIndel[1] = hqCount;
            deepestIndel[2] = lqCount;
          } else if (hqCount + lqCount > secondDeepestIndel[1] + secondDeepestIndel[2]) {
            if (indelInfo) {
              _indelInfoIndex1 = indelInfo->size;
              checkIndelInfoSpace(indelInfo, (4 + (insertLength + 3) / 4));
              indelInfo->array[indelInfo->size++] = 'I';
              indelInfo->array[indelInfo->size++] = insertLength;
              for (unsigned int i = _index; i < _index + (insertLength + 3) / 4; ++i) {
                indelInfo->array[indelInfo->size++] = insertInfo[i];
              }
              indelInfo->array[indelInfo->size++] = hqCount & 0xff;
              indelInfo->array[indelInfo->size++] = (hqCount >> 8) & 0xff;
              indelInfo->array[indelInfo->size++] = lqCount & 0xff;
              indelInfo->array[indelInfo->size++] = (lqCount >> 8) & 0xff;
              _indelInfoLength1 = indelInfo->size - _indelInfoIndex1;
            }

            secondDeepestIndel[0] = hqCount * 4; 
            secondDeepestIndel[1] = hqCount;
            secondDeepestIndel[2] = lqCount;
          }
        }
      } else {
        unsigned char insertLength = ((snpCounter[positionInsertion].weightedCount1 >> 8) & 0xF);
        unsigned char hqCount = (snpCounter[positionInsertion].weightedCount1 & 0xFF);
        unsigned char lqCount = snpOverflowCounterArray[region].counters[insertionIndex].insertion.counters.lqCount;
        indelHqCount += hqCount;
        indelLqCount += lqCount;
        if (hqCount + lqCount > deepestIndel[1] + deepestIndel[2]) {
          if (indelInfo) {
            _indelInfoIndex1 = _indelInfoIndex0;
            _indelInfoLength1 = _indelInfoLength0;
            _indelInfoIndex0 = indelInfo->size;
            checkIndelInfoSpace(indelInfo, (10 + (insertLength + 3) / 4));
            indelInfo->array[indelInfo->size++] = 'I';
            indelInfo->array[indelInfo->size++] = insertLength;
            for (unsigned int i = 0; i < (insertLength + 3) / 4; ++i) {
              indelInfo->array[indelInfo->size++] = snpOverflowCounterArray[region].counters[insertionIndex].insertion.counters.insertSeq[i];
            }
            indelInfo->array[indelInfo->size++] = hqCount & 0xff;
            indelInfo->array[indelInfo->size++] = (hqCount >> 8) & 0xff;
            indelInfo->array[indelInfo->size++] = lqCount & 0xff;
            indelInfo->array[indelInfo->size++] = (lqCount >> 8) & 0xff;
            _indelInfoLength0 = indelInfo->size - _indelInfoIndex0;
          }

          secondDeepestIndel[0] = deepestIndel[0];
          secondDeepestIndel[1] = deepestIndel[1];
          secondDeepestIndel[2] = deepestIndel[2];
          deepestIndel[0] = hqCount * 4; 
          deepestIndel[1] = hqCount;
          deepestIndel[2] = lqCount;
        } else if (hqCount + lqCount > secondDeepestIndel[1] + secondDeepestIndel[2]) {
          if (indelInfo) {
            _indelInfoIndex1 = indelInfo->size;
            checkIndelInfoSpace(indelInfo, (8 + (insertLength + 3) / 4));
            indelInfo->array[indelInfo->size++] = 'I';
            indelInfo->array[indelInfo->size++] = insertLength;
            for (unsigned int i = 0; i < (insertLength + 3) / 4; ++i) {
              indelInfo->array[indelInfo->size++] = snpOverflowCounterArray[region].counters[insertionIndex].insertion.counters.insertSeq[i];
            }
            indelInfo->array[indelInfo->size++] = hqCount & 0xff;
            indelInfo->array[indelInfo->size++] = (hqCount >> 8) & 0xff;
            indelInfo->array[indelInfo->size++] = lqCount & 0xff;
            indelInfo->array[indelInfo->size++] = (lqCount >> 8) & 0xff;
            _indelInfoLength1 = indelInfo->size - _indelInfoIndex1;
          }

          secondDeepestIndel[0] = hqCount * 4; 
          secondDeepestIndel[1] = hqCount;
          secondDeepestIndel[2] = lqCount;
        }
      }
    }

    if (snpCounter[positionDeletion].weightedCount2 & 0x2000) {
      if (snpCounter[positionDeletion].weightedCount2 & 0x1000) {
        unsigned short arraySize = (snpCounter[positionDeletion].weightedCount2 & 0xFFF);
        unsigned char *deleteInfo =
            (unsigned char *) getAddress(pool, snpOverflowCounterArray[region].counters[deletionIndex].deletion.ptr);
        unsigned int infoIndex = 0;
        while (infoIndex < arraySize) {
          unsigned char deleteLength = deleteInfo[infoIndex++];
          unsigned short hqCount;
          hqCount = deleteInfo[infoIndex++];
          hqCount |= deleteInfo[infoIndex++] << 8;
          unsigned short lqCount;
          lqCount = deleteInfo[infoIndex++];
          lqCount |= deleteInfo[infoIndex++] << 8;

          indelHqCount += hqCount;
          indelLqCount += lqCount;
          if (hqCount + lqCount > deepestIndel[1] + deepestIndel[2]) {
            if (indelInfo) {
              _indelInfoIndex1 = _indelInfoIndex0;
              _indelInfoLength1 = _indelInfoLength0;
              _indelInfoIndex0 = indelInfo->size;
              checkIndelInfoSpace(indelInfo, (4 + (deleteLength + 3) / 4));
              indelInfo->array[indelInfo->size++] = 'D';
              indelInfo->array[indelInfo->size++] = deleteLength;
              memset(&(indelInfo->array[indelInfo->size]), 0, (deleteLength + 3) / 4 * sizeof(unsigned char));
              for (unsigned int i = 0; i < deleteLength; ++i) {
                unsigned int _arrayIdx = i >> 2;
                unsigned int _arrayPos = i & 3;
                unsigned int _deletePos = positionDeletion + i;
                unsigned char base = ((packedSeq[_deletePos >> 4] >> (30 - ((_deletePos & 15) << 1))) & 3);
                indelInfo->array[indelInfo->size + _arrayIdx] |= (base << (_arrayPos << 1));
              }
              indelInfo->size += (deleteLength + 3) / 4;
              indelInfo->array[indelInfo->size++] = hqCount & 0xff;
              indelInfo->array[indelInfo->size++] = (hqCount >> 8) & 0xff;
              indelInfo->array[indelInfo->size++] = lqCount & 0xff;
              indelInfo->array[indelInfo->size++] = (lqCount >> 8) & 0xff;
              _indelInfoLength0 = indelInfo->size - _indelInfoIndex0;
            }

            secondDeepestIndel[0] = deepestIndel[0];
            secondDeepestIndel[1] = deepestIndel[1];
            secondDeepestIndel[2] = deepestIndel[2];
            deepestIndel[0] = hqCount * 4; 
            deepestIndel[1] = hqCount;
            deepestIndel[2] = lqCount;
          } else if (hqCount + lqCount > secondDeepestIndel[1] + secondDeepestIndel[2]) {
            if (indelInfo) {
              _indelInfoIndex1 = indelInfo->size;
              checkIndelInfoSpace(indelInfo, (4 + (deleteLength + 3) / 4));
              indelInfo->array[indelInfo->size++] = 'D';
              indelInfo->array[indelInfo->size++] = deleteLength;
              memset(&(indelInfo->array[indelInfo->size]), 0, (deleteLength + 3) / 4 * sizeof(unsigned char));
              for (unsigned int i = 0; i < deleteLength; ++i) {
                unsigned int _arrayIdx = i >> 2;
                unsigned int _arrayPos = i & 3;
                unsigned int _deletePos = positionDeletion + i;
                unsigned char base = ((packedSeq[_deletePos >> 4] >> (30 - ((_deletePos & 15) << 1))) & 3);
                indelInfo->array[indelInfo->size + _arrayIdx] |= (base << (_arrayPos << 1));
              }
              indelInfo->size += (deleteLength + 3) / 4;
              indelInfo->array[indelInfo->size++] = hqCount & 0xff;
              indelInfo->array[indelInfo->size++] = (hqCount >> 8) & 0xff;
              indelInfo->array[indelInfo->size++] = lqCount & 0xff;
              indelInfo->array[indelInfo->size++] = (lqCount >> 8) & 0xff;
              _indelInfoLength1 = indelInfo->size - _indelInfoIndex1;
            }

            secondDeepestIndel[0] = hqCount * 4; 
            secondDeepestIndel[1] = hqCount;
            secondDeepestIndel[2] = lqCount;
          }
        }
      } else {
        unsigned char deleteLength = snpOverflowCounterArray[region].counters[deletionIndex].deletion.counters.deletionLength;
        unsigned char hqCount = snpOverflowCounterArray[region].counters[deletionIndex].deletion.counters.hqCount;
        unsigned char lqCount = snpOverflowCounterArray[region].counters[deletionIndex].deletion.counters.lqCount;
        indelHqCount += hqCount;
        indelLqCount += lqCount;
        if (hqCount + lqCount > deepestIndel[1] + deepestIndel[2]) {
          if (indelInfo) {
            _indelInfoIndex1 = _indelInfoIndex0;
            _indelInfoLength1 = _indelInfoLength0;
            _indelInfoIndex0 = indelInfo->size;
            checkIndelInfoSpace(indelInfo, (10 + (deleteLength + 3) / 4));
            indelInfo->array[indelInfo->size++] = 'D';
            indelInfo->array[indelInfo->size++] = deleteLength;
            memset(&(indelInfo->array[indelInfo->size]), 0, (deleteLength + 3) / 4 * sizeof(unsigned char));
            for (unsigned int i = 0; i < deleteLength; ++i) {
              unsigned int _arrayIdx = i >> 2;
              unsigned int _arrayPos = i & 3;
              unsigned int _deletePos = positionDeletion + i;
              unsigned char base = ((packedSeq[_deletePos >> 4] >> (30 - ((_deletePos & 15) << 1))) & 3);
              indelInfo->array[indelInfo->size + _arrayIdx] |= (base << (_arrayPos << 1));
            }
            indelInfo->size += (deleteLength + 3) / 4;
            indelInfo->array[indelInfo->size++] = hqCount & 0xff;
            indelInfo->array[indelInfo->size++] = (hqCount >> 8) & 0xff;
            indelInfo->array[indelInfo->size++] = lqCount & 0xff;
            indelInfo->array[indelInfo->size++] = (lqCount >> 8) & 0xff;
            _indelInfoLength0 = indelInfo->size - _indelInfoIndex0;
          }

          secondDeepestIndel[0] = deepestIndel[0];
          secondDeepestIndel[1] = deepestIndel[1];
          secondDeepestIndel[2] = deepestIndel[2];
          deepestIndel[0] = hqCount * 4; 
          deepestIndel[1] = hqCount;
          deepestIndel[2] = lqCount;
        } else if (hqCount + lqCount > secondDeepestIndel[1] + secondDeepestIndel[2]) {
          if (indelInfo) {
            _indelInfoIndex1 = indelInfo->size;
            checkIndelInfoSpace(indelInfo, (10 + (deleteLength + 3) / 4));
            indelInfo->array[indelInfo->size++] = 'D';
            indelInfo->array[indelInfo->size++] = deleteLength;
            memset(&(indelInfo->array[indelInfo->size]), 0, (deleteLength + 3) / 4 * sizeof(unsigned char));
            for (unsigned int i = 0; i < deleteLength; ++i) {
              unsigned int _arrayIdx = i >> 2;
              unsigned int _arrayPos = i & 3;
              unsigned int _deletePos = positionDeletion + i;
              unsigned char base = ((packedSeq[_deletePos >> 4] >> (30 - ((_deletePos & 15) << 1))) & 3);
              indelInfo->array[indelInfo->size + _arrayIdx] |= (base << (_arrayPos << 1));
            }
            indelInfo->size += (deleteLength + 3) / 4;
            indelInfo->array[indelInfo->size++] = hqCount & 0xff;
            indelInfo->array[indelInfo->size++] = (hqCount >> 8) & 0xff;
            indelInfo->array[indelInfo->size++] = lqCount & 0xff;
            indelInfo->array[indelInfo->size++] = (lqCount >> 8) & 0xff;
            _indelInfoLength1 = indelInfo->size - _indelInfoIndex1;
          }

          secondDeepestIndel[0] = hqCount * 4; 
          secondDeepestIndel[1] = hqCount;
          secondDeepestIndel[2] = lqCount;
        }
      }
    }

    iStat.W[ALPHABET_SIZE] = deepestIndel[0];
    iStat.F[ALPHABET_SIZE] = deepestIndel[1];
    iStat.R[ALPHABET_SIZE] = deepestIndel[2];
    iStat.W[ALPHABET_SIZE + 1] = secondDeepestIndel[0];
    iStat.F[ALPHABET_SIZE + 1] = secondDeepestIndel[1];
    iStat.R[ALPHABET_SIZE + 1] = secondDeepestIndel[2];

    if (position == chkPos) {
      fprintf(stderr, "[%s-%d] %u, %u, %u\n%u, %u, %u\n%u, %u, %u\n%u, %u, %u\n%u, %u, %u\n%u, %u, %u\n%u, %u, %u\n\n",
              __FILE__, __LINE__,
              position, indelHqCount, indelLqCount,
              iStat.W[0], iStat.F[0], iStat.R[0],
              iStat.W[1], iStat.F[1], iStat.R[1],
              iStat.W[2], iStat.F[2], iStat.R[2],
              iStat.W[3], iStat.F[3], iStat.R[3],
              iStat.W[4], iStat.F[4], iStat.R[4],
              iStat.W[5], iStat.F[5], iStat.R[5]
      );
    }

    if (position == chkPos) {
      fprintf(stderr, "[%s-%d] %u, %u, %u\n%u, %u, %u\n%u, %u, %u\n%u, %u, %u\n%u, %u, %u\n%u, %u, %u\n%u, %u, %u\n\n",
              __FILE__, __LINE__,
              position, indelHqCount, indelLqCount,
              iStat.W[0], iStat.F[0], iStat.R[0],
              iStat.W[1], iStat.F[1], iStat.R[1],
              iStat.W[2], iStat.F[2], iStat.R[2],
              iStat.W[3], iStat.F[3], iStat.R[3],
              iStat.W[4], iStat.F[4], iStat.R[4],
              iStat.W[5], iStat.F[5], iStat.R[5]
      );
    }

    if (secondDeepestIndel[1] | secondDeepestIndel[2]) {
      if (indelInfo && _indelInfoCurrIndex == indelInfo->size) {
        fprintf(stderr, "Error in getting Metrics\n");
        exit(1);
      }
      if (indelInfo) {
        
        
        if (_indelInfoIndex0 == _indelInfoCurrIndex) {
          _indelInfoCurrIndex = _indelInfoIndex0 + _indelInfoLength0;
          if (_indelInfoIndex1 == _indelInfoCurrIndex) {
            _indelInfoCurrIndex = _indelInfoIndex1 + _indelInfoLength1;
          } else {
            for (unsigned int i = 0; i < _indelInfoLength1; ++i) {
              indelInfo->array[_indelInfoCurrIndex++] = indelInfo->array[_indelInfoIndex1++];
            }
          }
        } else {
          if (_indelInfoCurrIndex + _indelInfoLength0 > _indelInfoIndex1) {
            unsigned char *tmpArray = (unsigned char *) malloc(_indelInfoLength1 * sizeof(unsigned char));
            for (unsigned int i = 0; i < _indelInfoLength1; ++i) {
              tmpArray[i] = indelInfo->array[_indelInfoIndex1++];
            }

            for (unsigned int i = 0; i < _indelInfoLength0; ++i) {
              indelInfo->array[_indelInfoCurrIndex++] = indelInfo->array[_indelInfoIndex0++];
            }
            for (unsigned int i = 0; i < _indelInfoLength1; ++i) {
              indelInfo->array[_indelInfoCurrIndex++] = tmpArray[i];
            }
            free(tmpArray);
          } else {
            for (unsigned int i = 0; i < _indelInfoLength0; ++i) {
              indelInfo->array[_indelInfoCurrIndex++] = indelInfo->array[_indelInfoIndex0++];
            }
            for (unsigned int i = 0; i < _indelInfoLength1; ++i) {
              indelInfo->array[_indelInfoCurrIndex++] = indelInfo->array[_indelInfoIndex1++];
            }
          }
        }
        indelInfo->size = _indelInfoCurrIndex;
      }
      return 3; 
    }

    if (deepestIndel[1] | deepestIndel[2]) {
      if (indelInfo && _indelInfoCurrIndex == indelInfo->size) {
        fprintf(stderr, "Error in getting Metrics\n");
        exit(1);
      }
      return 2; 
    }
#endif
  }
  return 1;
}

void selectPossibleSNPs(SnpBundle snpBundle,
                        unsigned int *packedSeq, unsigned int dnaLength,
                        Annotation *annotation,
                        unsigned int *attriSize,
                        char isExome, ExomeRegion *exomeRegion, unsigned int numExomeRegion,
                        LikelihoodCache *likelihood_cache,
                        InputOptions &input_options, IniParams ini_params,
                        unsigned int *ambMap, Translate *translate,
                        const char *snp_noRF_filename,
                        double startTime, double &lastEventTime) {
  
  
  
  unsigned int genotype2int[16] =
      {
          0, 1, 2, 3, 9, 10, 11, 18, 19, 27, 4, 12, 20, 28, 36, 37
      };
  copyGenotypeTableToGPU(genotype2int);

  double *snpPrior = createSnpPriorArray();
  copySnpPriorToGPU(snpPrior);

  FILE *snpOutput = (FILE *) fopen(snp_noRF_filename, "wb");

  double cpu2gpuFisherTime;
  double totalCpu2gpuFisherTime = 0.0;
  double fisherTime;
  double totalFisherTime = 0.0;
  double gpu2cpuFisherTime;
  double totalGpu2cpuFisherTime = 0.0;
  double calMetaTime;
  double totalCalMetaTime = 0.0;

  unsigned int threadsPerBlock = 256; 
  unsigned int numBlock = 32768; 
  unsigned int batchSize = threadsPerBlock * numBlock;

  
  InputStat *iStat0 = (InputStat *) malloc(batchSize * sizeof(InputStat));
  memset(iStat0, 0, batchSize * sizeof(InputStat));
  InputStat *iStat1 = (InputStat *) malloc(batchSize * sizeof(InputStat));
  memset(iStat1, 0, batchSize * sizeof(InputStat));
  InputStat *iStatNext = (InputStat *) malloc(batchSize * sizeof(InputStat));
  memset(iStatNext, 0, batchSize * sizeof(InputStat));
  MetaReference *reference0 = (MetaReference *) malloc(batchSize * sizeof(MetaReference));
  MetaReference *reference1 = (MetaReference *) malloc(batchSize * sizeof(MetaReference));
  MetaReference *referenceNext = (MetaReference *) malloc(batchSize * sizeof(MetaReference));
  unsigned char *refChar = (unsigned char *) malloc(batchSize * sizeof(unsigned char));
  GenotypeLikelihood *genotypeLikelihood0 = (GenotypeLikelihood *) malloc(batchSize * sizeof(GenotypeLikelihood));
  GenotypeLikelihood *genotypeLikelihood1 = (GenotypeLikelihood *) malloc(batchSize * sizeof(GenotypeLikelihood));
  GenotypeLikelihood *genotypeLikelihoodNext = (GenotypeLikelihood *) malloc(batchSize * sizeof(GenotypeLikelihood));
  SnpCallingInfo *scInfo0 = (SnpCallingInfo *) malloc(batchSize * sizeof(SnpCallingInfo));
  SnpCallingInfo *scInfo1 = (SnpCallingInfo *) malloc(batchSize * sizeof(SnpCallingInfo));
  SnpCallingInfo *scInfoNext = (SnpCallingInfo *) malloc(batchSize * sizeof(SnpCallingInfo));
  MetaWindowInfo *mwInfo0 = (MetaWindowInfo *) malloc((batchSize + SNP_META_WINDOW_SIZE * 2) * sizeof(MetaWindowInfo));
  MetaWindowInfo *mwInfo1 = (MetaWindowInfo *) malloc((batchSize + SNP_META_WINDOW_SIZE * 2) * sizeof(MetaWindowInfo));
  MetaWindowInfo *mwInfoNext = (MetaWindowInfo *) malloc(
      (batchSize + SNP_META_WINDOW_SIZE * 2) * sizeof(MetaWindowInfo));
  IndelInfo *indelInfo0 = (IndelInfo *) malloc(sizeof(IndelInfo));
  indelInfo0->array = (unsigned char *) malloc(batchSize * sizeof(unsigned char));
  memset(indelInfo0->array, 0, batchSize * sizeof(unsigned char));
  indelInfo0->size = 0;
  indelInfo0->limit = batchSize;
  IndelInfo *indelInfo1 = (IndelInfo *) malloc(sizeof(IndelInfo));
  indelInfo1->array = (unsigned char *) malloc(batchSize * sizeof(unsigned char));
  memset(indelInfo1->array, 0, batchSize * sizeof(unsigned char));
  indelInfo1->size = 0;
  indelInfo1->limit = batchSize;
  IndelInfo *indelInfoNext = (IndelInfo *) malloc(sizeof(IndelInfo));
  indelInfoNext->array = (unsigned char *) malloc(batchSize * sizeof(unsigned char));
  memset(indelInfoNext->array, 0, batchSize * sizeof(unsigned char));
  indelInfoNext->size = 0;
  indelInfoNext->limit = batchSize;

  SNPMetaBuffer *buffer0 = createSNPMetaBuffer(reference0, iStat0, genotypeLikelihood0,
                                               scInfo0, mwInfo0, indelInfo0);
  SNPMetaBuffer *buffer1 = createSNPMetaBuffer(reference1, iStat1, genotypeLikelihood1,
                                               scInfo1, mwInfo1, indelInfo1);
  SNPMetaBuffer *bufferNext = createSNPMetaBuffer(referenceNext, iStatNext, genotypeLikelihoodNext,
                                                  scInfoNext, mwInfoNext, indelInfoNext);

  SNPMetaBuffer *buffer;
  MetaReference *reference;
  InputStat *iStat;
  GenotypeLikelihood *genotypeLikelihood;
  SnpCallingInfo *scInfo;
  MetaWindowInfo *mwInfo;
  IndelInfo *indelInfo;
  pthread_t filterSnpThread;

  startFilterSNPThread(buffer0, buffer1, snpOutput, attriSize, filterSnpThread);

  buffer = getIdleSNPMetaBuffer(buffer0, buffer1);
  reference = buffer->reference;
  iStat = buffer->snpCounter;
  genotypeLikelihood = buffer->genotypeLikelihood;
  scInfo = buffer->scInfo;
  mwInfo = buffer->mwInfo;
  indelInfo = buffer->indelInfo;
  memset(mwInfo, 0, (batchSize + SNP_META_WINDOW_SIZE * 2) * sizeof(MetaWindowInfo));

  InputStat inputStat;
  unsigned int weightedCount;
  unsigned int posStrandCount;
  unsigned int negStrandCount;
  unsigned char gcStrandCount;
  unsigned char gStrandCount;
  unsigned int indelHqCount;
  unsigned int indelLqCount;
  char isValid;

  unsigned int i, j;
  unsigned int exomeGroup = 0;
  unsigned int batchIndex = 0;
  unsigned int position;

  if (!isExome) {
    exomeRegion = (ExomeRegion *) malloc(sizeof(ExomeRegion));
    exomeRegion[0].startPos = 0;
    exomeRegion[0].endPos = dnaLength;
    numExomeRegion = 1;
  }

  if (isExome && input_options.verbose) {
    printf("[BALSA] Number of exome groups : %u\n", numExomeRegion);
  }

  for (exomeGroup = 0; exomeGroup < numExomeRegion; ++exomeGroup) {
    unsigned int start = exomeRegion[exomeGroup].startPos;
    unsigned int end = exomeRegion[exomeGroup].endPos;
    if (start >= dnaLength || end > dnaLength) {
      fprintf(stderr,
              "[Error] Positions in exome region are incompatible with the index, please generate the exome list with the correct index.\n");
      exit(1);
    }
    for (position = start; position < end; ++position) {
      char refBase = ((packedSeq[position >> 4] >> (30 - ((position & 15) << 1))) & 3);
      isValid = getSnpMetrics(packedSeq, snpBundle, position, dnaLength, inputStat,
                              weightedCount, posStrandCount, negStrandCount, gcStrandCount,
                              gStrandCount, indelHqCount, indelLqCount, indelInfo);

      reference[batchIndex].amb = position;
      reference[batchIndex].refChar = refBase;
      unsigned short chrID;
      getChrAndPos(ambMap, translate, position, &(reference[batchIndex].chrPos), &(chrID));
      strcpy(reference[batchIndex].chrName, annotation[chrID - 1].text);

      iStat[batchIndex] = inputStat;

      mwInfo[SNP_META_WINDOW_SIZE + batchIndex].weightedCount = weightedCount;
      mwInfo[SNP_META_WINDOW_SIZE + batchIndex].posStrandCount = posStrandCount;
      mwInfo[SNP_META_WINDOW_SIZE + batchIndex].negStrandCount = negStrandCount;
      mwInfo[SNP_META_WINDOW_SIZE + batchIndex].gcStrandCount = gcStrandCount;
      mwInfo[SNP_META_WINDOW_SIZE + batchIndex].gStrandCount = gStrandCount;
      mwInfo[SNP_META_WINDOW_SIZE + batchIndex].indelHqCount = indelHqCount;
      mwInfo[SNP_META_WINDOW_SIZE + batchIndex].indelLqCount = indelLqCount;
      mwInfo[SNP_META_WINDOW_SIZE + batchIndex].isValid = isValid;
      if ((batchIndex == 0 && position > 0) || (position == start && position > 0)) {
        char preBase = ((packedSeq[(position - 1) >> 4] >> (30 - (((position - 1) & 15) << 1))) & 3);
        mwInfo[SNP_META_WINDOW_SIZE + batchIndex].preBase = preBase;
      } else if (batchIndex == 0) {
        mwInfo[SNP_META_WINDOW_SIZE + batchIndex].preBase = 4; 
      } else {
        mwInfo[SNP_META_WINDOW_SIZE + batchIndex].preBase = reference[batchIndex - 1].refChar;
      }

      unsigned short polyrun = 1;
      unsigned int testPos = position;
      char testBase;
      while (testPos > 0) {
        testBase = ((packedSeq[(testPos - 1) >> 4] >> (30 - (((testPos - 1) & 15) << 1))) & 3);
        if (testBase == refBase) {
          if (polyrun < 0xFFFF) {
            polyrun++;
          }
          testPos--;
        } else {
          break;
        }
      }
      testPos = position + 1;
      while (testPos < dnaLength) {
        testBase = ((packedSeq[testPos >> 4] >> (30 - ((testPos & 15) << 1))) & 3);
        if (testBase == refBase) {
          if (polyrun < 0xFFFF) {
            polyrun++;
          }
          testPos++;
        } else {
          break;
        }
      }
      mwInfo[SNP_META_WINDOW_SIZE + batchIndex].polyrun = polyrun;

      batchIndex++;

      if (batchIndex == batchSize) {
        for (j = 0; j < batchIndex; ++j) {
          refChar[j] = reference[j].refChar;
        }
        
        cpu2gpuFisherTime = getElapsedTime(startTime);
        totalCpu2gpuFisherTime += cpu2gpuFisherTime - lastEventTime;
        lastEventTime = cpu2gpuFisherTime;

        
        
        computeGenotypeLikelihoodWrapper(iStat, reference, genotypeLikelihood, ini_params, batchIndex,
                                         likelihood_cache);
        computeGenotypeWrapper(refChar, genotypeLikelihood, scInfo, ini_params, batchIndex);

        fisherTime = getElapsedTime(startTime);
        totalFisherTime += fisherTime - lastEventTime;
        lastEventTime = fisherTime;

        gpu2cpuFisherTime = getElapsedTime(startTime);
        totalGpu2cpuFisherTime += gpu2cpuFisherTime - lastEventTime;
        lastEventTime = gpu2cpuFisherTime;

        unsigned int mwIdx = 0;
        if (reference[0].amb > SNP_META_WINDOW_SIZE) {
          for (i = reference[0].amb - SNP_META_WINDOW_SIZE; i < reference[0].amb; ++i) {
            isValid = getSnpMetrics(packedSeq, snpBundle, i, dnaLength, inputStat,
                                    weightedCount, posStrandCount, negStrandCount, gcStrandCount,
                                    gStrandCount, indelHqCount, indelLqCount, NULL);

            mwInfo[mwIdx].weightedCount = weightedCount;
            mwInfo[mwIdx].posStrandCount = posStrandCount;
            mwInfo[mwIdx].negStrandCount = negStrandCount;
            mwInfo[mwIdx].gcStrandCount = gcStrandCount;
            mwInfo[mwIdx].gStrandCount = gStrandCount;
            mwInfo[mwIdx].indelHqCount = indelHqCount;
            mwInfo[mwIdx].indelLqCount = indelLqCount;
            mwInfo[mwIdx].isValid = isValid;
            mwIdx++;
          }
        }

        mwIdx = SNP_META_WINDOW_SIZE + batchIndex;
        for (i = reference[batchIndex - 1].amb + 1;
             i < dnaLength && i < reference[batchIndex - 1].amb + SNP_META_WINDOW_SIZE; ++i) {
          isValid = getSnpMetrics(packedSeq, snpBundle, i, dnaLength, inputStat,
                                  weightedCount, posStrandCount, negStrandCount, gcStrandCount,
                                  gStrandCount, indelHqCount, indelLqCount, NULL);

          mwInfo[mwIdx].weightedCount = weightedCount;
          mwInfo[mwIdx].posStrandCount = posStrandCount;
          mwInfo[mwIdx].negStrandCount = negStrandCount;
          mwInfo[mwIdx].gcStrandCount = gcStrandCount;
          mwInfo[mwIdx].gStrandCount = gStrandCount;
          mwInfo[mwIdx].indelHqCount = indelHqCount;
          mwInfo[mwIdx].indelLqCount = indelLqCount;
          mwInfo[mwIdx].isValid = isValid;
          mwIdx++;
        }

        setReadySNPMetaBufferStatus(buffer, batchIndex, 0);
        buffer = getIdleSNPMetaBuffer(buffer0, buffer1);
        reference = buffer->reference;
        iStat = buffer->snpCounter;
        genotypeLikelihood = buffer->genotypeLikelihood;
        scInfo = buffer->scInfo;
        mwInfo = buffer->mwInfo;
        indelInfo = buffer->indelInfo;
        memset(mwInfo, 0, (batchSize + SNP_META_WINDOW_SIZE * 2) * sizeof(MetaWindowInfo));

        calMetaTime = getElapsedTime(startTime);
        totalCalMetaTime += calMetaTime - lastEventTime;
        lastEventTime = calMetaTime;

        memset(iStat, 0, batchSize * sizeof(InputStat));
        indelInfo->size = 0;
        memset(indelInfo->array, 0, indelInfo->limit * sizeof(unsigned char));
        batchIndex = 0;
      }
    }

    if (batchIndex > 0 && exomeGroup == numExomeRegion - 1) {
      for (j = 0; j < batchIndex; ++j) {
        refChar[j] = reference[j].refChar;
      }

      
      cpu2gpuFisherTime = getElapsedTime(startTime);
      totalCpu2gpuFisherTime += cpu2gpuFisherTime - lastEventTime;
      lastEventTime = cpu2gpuFisherTime;

      
      
      computeGenotypeLikelihoodWrapper(iStat, reference, genotypeLikelihood, ini_params, batchIndex, likelihood_cache);
      computeGenotypeWrapper(refChar, genotypeLikelihood, scInfo, ini_params, batchIndex);

      fisherTime = getElapsedTime(startTime);
      totalFisherTime += fisherTime - lastEventTime;
      lastEventTime = fisherTime;

      
      gpu2cpuFisherTime = getElapsedTime(startTime);
      totalGpu2cpuFisherTime += gpu2cpuFisherTime - lastEventTime;
      lastEventTime = gpu2cpuFisherTime;

      unsigned int mwIdx = 0;
      if (reference[0].amb > SNP_META_WINDOW_SIZE) {
        for (i = reference[0].amb - SNP_META_WINDOW_SIZE; i < reference[0].amb; ++i) {
          isValid = getSnpMetrics(packedSeq, snpBundle, i, dnaLength, inputStat,
                                  weightedCount, posStrandCount, negStrandCount, gcStrandCount,
                                  gStrandCount, indelHqCount, indelLqCount, NULL);

          mwInfo[mwIdx].weightedCount = weightedCount;
          mwInfo[mwIdx].posStrandCount = posStrandCount;
          mwInfo[mwIdx].negStrandCount = negStrandCount;
          mwInfo[mwIdx].gcStrandCount = gcStrandCount;
          mwInfo[mwIdx].gStrandCount = gStrandCount;
          mwInfo[mwIdx].indelHqCount = indelHqCount;
          mwInfo[mwIdx].indelLqCount = indelLqCount;
          mwInfo[mwIdx].isValid = isValid;
          mwIdx++;
        }
      }

      mwIdx = SNP_META_WINDOW_SIZE + batchIndex;
      for (i = reference[batchIndex - 1].amb + 1;
           i < dnaLength && i < reference[batchIndex - 1].amb + SNP_META_WINDOW_SIZE; ++i) {
        isValid = getSnpMetrics(packedSeq, snpBundle, i, dnaLength, inputStat,
                                weightedCount, posStrandCount, negStrandCount, gcStrandCount,
                                gStrandCount, indelHqCount, indelLqCount, NULL);

        mwInfo[mwIdx].weightedCount = weightedCount;
        mwInfo[mwIdx].posStrandCount = posStrandCount;
        mwInfo[mwIdx].negStrandCount = negStrandCount;
        mwInfo[mwIdx].gcStrandCount = gcStrandCount;
        mwInfo[mwIdx].gStrandCount = gStrandCount;
        mwInfo[mwIdx].indelHqCount = indelHqCount;
        mwInfo[mwIdx].indelLqCount = indelLqCount;
        mwInfo[mwIdx].isValid = isValid;
        mwIdx++;
      }

      setReadySNPMetaBufferStatus(buffer, batchIndex, 0);
      waitFinishSNPMetaBuffer(buffer);
      setFinishSNPMetaBufferStatus(buffer0);
      setFinishSNPMetaBufferStatus(buffer1);

      calMetaTime = getElapsedTime(startTime);
      totalCalMetaTime += calMetaTime - lastEventTime;
      lastEventTime = calMetaTime;
    }
  }

  if (!isExome) {
    free(exomeRegion);
  }

  closeFilterSNPThread(filterSnpThread);
  freeSNPMetaBuffer(buffer0);
  freeSNPMetaBuffer(buffer1);
  freeSNPMetaBuffer(bufferNext);
  free(iStat0);
  free(iStat1);
  free(iStatNext);
  free(reference0);
  free(reference1);
  free(referenceNext);
  free(refChar);
  free(genotypeLikelihood0);
  free(genotypeLikelihood1);
  free(genotypeLikelihoodNext);
  free(scInfo0);
  free(scInfo1);
  free(scInfoNext);
  free(mwInfo0);
  free(mwInfo1);
  free(mwInfoNext);
  free(indelInfo0->array);
  free(indelInfo0);
  free(indelInfo1->array);
  free(indelInfo1);
  free(indelInfoNext->array);
  free(indelInfoNext);

  freeSnpPriorArray(snpPrior);

  fclose(snpOutput);

  if (input_options.verbose) {
    printf("Total Likelihood and Bias computing time : %9.4f seconds.\n", totalFisherTime);
    printf("Total time of SNP calling module : %9.4f seconds.\n", totalCalMetaTime);
  }
}



void calRFpedictProb(const char *snp_noRF_filename,
                     const char *snp_filename,
                     unsigned int attriSize, char verbose) {
  calRandomForestpPedictProb(snp_noRF_filename, snp_filename, attriSize, verbose);
}
