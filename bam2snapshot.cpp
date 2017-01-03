#include "bam2snapshot.h"

static void printTimeElapsed(const time_t startTime) {
  if (startTime) {
    printf("Time elapsed by now is %.0f seconds.\n", difftime(time(NULL), startTime));
    fflush(stdout);
  }
}


void initAlignment(Alignment *align, int maxReadLength, int maxCigarLength) {
  align->startPos = 0;
  align->strand = 0;
  align->readLength = 0;
  align->read = (unsigned char *) calloc(maxReadLength + 1, sizeof(unsigned char));
  align->qualities = (char *) calloc(maxReadLength + 1, sizeof(char));
  align->cigar = (char *) calloc(maxCigarLength + 1, sizeof(char));
}


void freeAlignment(Alignment *align) {
  free(align->read);
  free(align->qualities);
  free(align->cigar);
}



void clearSnapshot(unsigned dnaLength, unsigned numOfCPUThreads, SnpBundle *snpBundle) {
  memset(snpBundle, 0, sizeof(snpBundle));

  snpBundle->snpCounter = (SnpCounter *) malloc(dnaLength * sizeof(SnpCounter));
  assert(snpBundle->snpCounter != NULL);
  memset(snpBundle->snpCounter, 0, dnaLength * sizeof(SnpCounter));

  snpBundle->numOfCPUThreads = numOfCPUThreads;

  constructSnpOverflowCounter(&(snpBundle->snpOverflowCounterArray), numOfCPUThreads);


  snpBundle->invalidSnpCounterPos = (unsigned *) calloc((dnaLength + 31) >> 5, sizeof(unsigned));

  snpBundle->snpMemoryPool = createPool(SNP_POOL_SIZE); 
}



void freeSnapshot(SnpBundle *snpBundle) {
  free(snpBundle->snpCounter);
  destroySnpOverflowCounter(snpBundle->snpOverflowCounterArray, snpBundle->numOfCPUThreads);
  free(snpBundle->invalidSnpCounterPos);
  destroyPool(snpBundle->snpMemoryPool);
}

int getLenFromCigar(char *cigarString) {
  char *cigar = cigarString;
  int cigarLength = strlen(cigar);
  int cigarInt;
  char cigarChar;
  int lengthScanned;
  int lenFromCigar = 0;

  if (cigarLength == 0) return 0;
  while (1) {
    sscanf(cigar, "%d%c%n", &cigarInt, &cigarChar, &lengthScanned);
    if (cigarChar == 'M' || cigarChar == 'I' || cigarChar == 'S') {
      lenFromCigar += cigarInt;
    }

    cigarLength -= lengthScanned;
    cigar += lengthScanned;
    if (cigarLength <= 0) {
      break;
    }
  }
  return lenFromCigar;
}

inline void bamToAlignment(bam1_t *bam, const bam_header_t *header, Alignment *align,
                           unsigned dnaLength, unsigned mapqThreshold,
                           Annotation *annotations, unsigned numOfSeq,
                           unsigned int *ambMap, Translate *translate, int chrIdMap[]) {

  if (bam->core.qual < mapqThreshold) {
    align->startPos = dnaLength + 1;    
    return;
  }

  int32_t tid = bam->core.tid;

  if (tid < 0 || tid >= header->n_targets || bam->core.pos < 0 || bam->core.pos == INT32_MAX) {
    align->startPos = dnaLength + 2;    
    return;
  }

  if (chrIdMap[tid] == -1) {
    chrIdMap[tid] = getChrIDFromName(annotations, numOfSeq, header->target_name[tid]);
    assert(chrIdMap[tid] >= 1);
    assert(chrIdMap[tid] <= numOfSeq);
  }

  unsigned posInChr = bam->core.pos + 1;    

  align->startPos = getAmbPos(chrIdMap[tid], posInChr, ambMap, translate, dnaLength);

#ifdef DEBUG
  assert(align->startPos < dnaLength);
#endif

  align->strand = ((bam1_strand(bam)) ? 1 : 0);

  align->readLength = bam->core.l_qseq;

  if (align->readLength <= 0) return;

  uint8_t *read = bam1_seq(bam);
  uint8_t *qual = bam1_qual(bam);

  if (!align->strand) {
    for (int i = 0; i < align->readLength; i++) {
      align->qualities[i] = qual[i];
      align->read[i] = bamcode2snpcode[(bam1_seqi(read, i))];
    }
  } else {
    for (int i = 0; i < align->readLength; i++) {
      int j = align->readLength - 1 - i;
      align->qualities[j] = qual[i];
      align->read[j] = 3 - bamcode2snpcode[(bam1_seqi(read, i))];
    }
  }

  align->qualities[align->readLength] = '\0';

  char *cur = align->cigar;
  int n_cigar = bam->core.n_cigar;
  for (int i = 0; i < n_cigar; i++) {
    uint32_t cigar_base = bam1_cigar(bam)[i];
    int print = sprintf(cur, "%d%c", bam_cigar_oplen(cigar_base), bam_cigar_opchr(cigar_base));
    cur += print;
  }

  sprintf(cur, "\0");
}

void updateSnapshotForFile(SnpBundle *snpBundle, VcSetting *vcSetting,
                           unsigned char recalScore[],
                           Annotation *annotations, int numOfSeq,
                           unsigned int *ambMap, Translate *translate,
                           FILE *overflowFile, samfile_t *samFilePtr, FILE *softclipFilePtr,
                           AlignmentQC &alignmentQC, unsigned int *exomePosBitVec) {

  bam1_t *bam = bam_init1();
  Alignment align;

  initAlignment(&align, Max_Read_Length, 4 * Max_Read_Length);

  int chrIdMap[Max_Num_Of_Chrs];
  memset(chrIdMap, 0xFF, sizeof(chrIdMap));

  long long int count = 0;

  while (samread(samFilePtr, bam) >= 0) {

    count++;


    bamToAlignment(bam, samFilePtr->header, &align, vcSetting->dnaLength, vcSetting->mapqThreshold, annotations,
                   numOfSeq, ambMap, translate, chrIdMap);

    if ((bam->core.flag & 0x400) == 0x400) {
      
      continue;
    }

    
    updateAlignmentQCForReads(alignmentQC, align.readLength, align.cigar, align.startPos,
                              align.startPos == vcSetting->dnaLength + 2 ? 0 : 1, exomePosBitVec);

    if (align.readLength <= 0 || align.readLength != getLenFromCigar(align.cigar)) {

      continue;
    }

    if (align.startPos >= vcSetting->dnaLength || align.cigar[0] == '\0') {
      continue;
    }

    int trimHead, trimTail;
    getTrimSize(NULL, align.read, align.readLength, align.strand == 1,
                align.cigar, align.startPos,
                &trimHead, &trimTail, vcSetting->trimSize);

    updateSnpCounterForReads(*snpBundle, vcSetting, trimHead, trimTail,
                             recalScore, NULL, 0, align.startPos, 1 - align.strand, align.read, align.qualities,
                             align.readLength, align.cigar, overflowFile, softclipFilePtr);


  }

  printf("Total number of alignments: %lld\n", count);

  freeAlignment(&align);

  bam_destroy1(bam);

}


void setEnvironment(SnpBundle *snpBundle) {
  
  constructSnpOverflowBufferArray(snpBundle->numOfCPUThreads);
}


void freeEnvironment(SnpBundle *snpBundle) {
  
  destroySnpOverflowBufferArray(snpBundle->numOfCPUThreads);
}


void bamFilesToSnapshot(SnpBundle *snpBundle, unsigned dnaLength,
                        char *inputFilenames[], int inputFileTypes[], int n_inputFiles, char *outputFilenamePrefix,
                        unsigned char recalScore[], unsigned char weightMap[], unsigned snpTrimSize,
                        unsigned softClipThreshold,
                        unsigned mapqThreshold, unsigned indelWeightThreshold,
                        Annotation *annotations, int numOfSeq, unsigned int *ambMap, Translate *translate,
                        AlignmentQC &alignmentQC, unsigned int *exomePosBitVec,
                        const time_t startTime) {

  
  setEnvironment(snpBundle);


  char *tmpOverFlowFilename = (char *) malloc(strlen(outputFilenamePrefix) + 20);
  sprintf(tmpOverFlowFilename, "%s.tmp", outputFilenamePrefix);
  
  FILE *overflowFile = fopen(tmpOverFlowFilename, "wb");
  assert(overflowFile);

  FILE *softclipFilePtr = tmpfile();
  assert(softclipFilePtr);


  
  printf("Updating snapshot...\n");
  fflush(stdout);

  VcSetting vcSetting;
  memset(&vcSetting, 0, sizeof(VcSetting));
  vcSetting.dnaLength = dnaLength;
  vcSetting.mapqThreshold = mapqThreshold;
  vcSetting.trimSize = snpTrimSize;
  vcSetting.softClipThreshold = softClipThreshold;
  vcSetting.weightMap = weightMap;
  vcSetting.indelWeightThreshold = indelWeightThreshold;
  vcSetting.enableQualityCorrection = 1;

  for (int i = 0; i < n_inputFiles; i++) {

    samfile_t *samFilePtr;
    if (inputFileTypes[i] == 1) {
      samFilePtr = samopen(inputFilenames[i], "rb", NULL);
    } else {
      samFilePtr = samopen(inputFilenames[i], "r", NULL);
    }

    if (!samFilePtr) {
      fprintf(stderr, "[ERROR] Failed accessing %s.\n",
              inputFilenames[i]);
      exit(1);
    }

    printf("Processing %s...\n", inputFilenames[i]);
    fflush(stdout);

    if (numOfSeq != samFilePtr->header->n_targets) {
      printf("[WARNING] %d chromosomes in reference, but %s is %d.\n", numOfSeq, inputFilenames[i],
             samFilePtr->header->n_targets);
      fflush(stdout);
    }

    updateSnapshotForFile(snpBundle, &vcSetting, recalScore,
                          annotations, numOfSeq, ambMap, translate, overflowFile, samFilePtr, softclipFilePtr,
                          alignmentQC, exomePosBitVec
    );

    printf("Done with %s.\n", inputFilenames[i]);
    fflush(stdout);

    samclose(samFilePtr);

    printTimeElapsed(startTime);
  }

  fclose(overflowFile);

  printf("Done with all files.\n");
  fflush(stdout);

  
  printf("Dealing with overflows...\n");
  fflush(stdout);

  readFilesForOverflow(*snpBundle, &vcSetting, tmpOverFlowFilename);

  printf("Done.\n");
  fflush(stdout);

  
  printf("Dealing with directional soft-clips...\n");

  FILE *tmpFile[5];
  for (int i = 0; i < 5; i++)
    tmpFile[i] = tmpfile();

  constructDirectionalSoftclipCounterFromFiles(1,
                                               &softclipFilePtr, tmpFile[0], tmpFile[1],
                                               &tmpFile[2], tmpFile[3], tmpFile[4],
                                               snpBundle->snpDirectionalSoftclipCounter,
                                               snpBundle->snpDirectionalSoftclipCounterSize,
                                               snpBundle->snpDirectionalSoftclipCounterCapacity);

  printf("Done.\n");
  fflush(stdout);

  fclose(softclipFilePtr);

  for (int i = 0; i < 5; i++)
    fclose(tmpFile[i]);

  
  freeEnvironment(snpBundle);
  remove(tmpOverFlowFilename);
}

void loadBam2snapshotParas(const char *iniFileName, unsigned *numOfCPUThreads, unsigned char recalScore[],
                           unsigned char weightMap[],
                           unsigned *snpTrimSize, unsigned *softClipThreshold,
                           unsigned *mapqThreshold, unsigned *indelWeightThreshold) {
  if (!fileExists(iniFileName)) {
    fprintf(stderr, "[ERROR] Configuration file %s doesn't exists.\n", iniFileName);
    exit(1);
  }

  IniParams iniParams;

  ParseIniFile((char *) iniFileName, iniParams);

  *numOfCPUThreads = iniParams.Ini_NumOfCpuThreads;

  
  for (int i = 0; i < Q_SCORE_ARRAY_SIZE; i++)
    recalScore[i] = i;

  
  for (int i = iniParams.Ini_weightBound0 + 1; i <= iniParams.Ini_weightBound1; i++)
    weightMap[i] = 1;
  for (int i = iniParams.Ini_weightBound1 + 1; i <= iniParams.Ini_weightBound2; i++)
    weightMap[i] = 2;
  for (int i = iniParams.Ini_weightBound2 + 1; i <= iniParams.Ini_weightBound3; i++)
    weightMap[i] = 3;
  for (int i = iniParams.Ini_weightBound3 + 1; i < Q_SCORE_ARRAY_SIZE; i++)
    weightMap[i] = 4;

  *snpTrimSize = iniParams.Ini_statTrimSize;
  *softClipThreshold = iniParams.Ini_statSoftClipThreshold;
  *mapqThreshold = iniParams.Ini_statMapqThreshold;

  *indelWeightThreshold = 0;
  if (iniParams.Ini_IndelQualityThreshold == 40) {
    *indelWeightThreshold = 4;
  } else if (iniParams.Ini_IndelQualityThreshold >= iniParams.Ini_weightBound3) {
    *indelWeightThreshold = 3;
  } else if (iniParams.Ini_IndelQualityThreshold >= iniParams.Ini_weightBound2) {
    *indelWeightThreshold = 2;
  } else if (iniParams.Ini_IndelQualityThreshold >= iniParams.Ini_weightBound1) {
    *indelWeightThreshold = 1;
  }
}

void printCommandUsage(char *programName) {
  fprintf(stderr, "Usage:\n");
  fprintf(stderr,
          "    %s -i <Reference Index Prefix> -o <Output Prefix> -b|-s <BAM/SAM> ... [-e regionListFile]\n",
          programName);
  fprintf(stderr, "    -o: Output Prefix\n");
  fprintf(stderr, "    -b|-s <BAM/SAM>: Can be specified multiple times\n");
  fprintf(stderr, "    -e <Exome Region Index>: Exon Region Index generated by RegionIndexBuilder\n");
}


void parseCommandlineArguments(int argn, char *args[], char *&indexPrefix, char *&outputFilenamePrefix,
                               char *inputFilenames[], int inputFileTypes[], int &n_inputFiles, char *&exomeFilename) {
  n_inputFiles = 0;
  indexPrefix = NULL;
  outputFilenamePrefix = NULL;

	if(argn <= 1) {
		printCommandUsage(args[0]);
		exit(1);
	}

  for (int i = 1; i < argn; i += 2) {
    if (i + 1 >= argn) {
      printCommandUsage(args[0]);
      exit(1);
    }

    char *flag = args[i];
    if (strcmp(flag, "-i") == 0) {
      if (indexPrefix) {
        fprintf(stderr, "Usage Error: Multiple indexes spcified.\n\n");
        printCommandUsage(args[0]);
        exit(1);
      }
      indexPrefix = args[i + 1];
    } else if (strcmp(flag, "-b") == 0 || strcmp(flag, "-s") == 0) {

      if (n_inputFiles >= Max_Num_Of_Input_Files) {
        fprintf(stderr, "Usage Error: Exceed maximum %d BAM/SAM files.\n", Max_Num_Of_Input_Files);
        exit(1);
      }

      inputFilenames[n_inputFiles] = args[i + 1];
      inputFileTypes[n_inputFiles++] = (flag[1] == 'b' ? 1 : 0);
    } else if (strcmp(flag, "-o") == 0) {
      if (outputFilenamePrefix) {
        fprintf(stderr, "Usage Error: Multiple Output Prefixes specified.\n");
        printCommandUsage(args[0]);
        exit(1);
      }
      outputFilenamePrefix = args[i + 1];
    } else if (strcmp(flag, "-e") == 0) {
      if (exomeFilename) {
        fprintf(stderr, "Usage Error: Multiple Exome Region Index specified.\n");
        printCommandUsage(args[0]);
        exit(1);
      }
      exomeFilename = args[i + 1];
    } else {
      fprintf(stderr, "Usage Error: Invalid arguments.\n");
      printCommandUsage(args[0]);
      exit(1);
    }
  }

  if (!indexPrefix) {
    fprintf(stderr, "Usage Error: Missing index.\n");
    printCommandUsage(args[0]);
    exit(1);
  }

  if (!outputFilenamePrefix) {
    fprintf(stderr, "Usage Error: Missing Output Prefix.\n");
    printCommandUsage(args[0]);
    exit(1);
  }

  if (n_inputFiles == 0) {
    fprintf(stderr, "Usage Error: Missing BAM/SAM.\n");
    printCommandUsage(args[0]);
    exit(1);
  }

  if (!dirOfPrefixExists(indexPrefix)) {
    fprintf(stderr, "Error: The directory of Reference Index doesn't exists.\n");
    exit(1);
  }

  if (!dirOfPrefixExists(outputFilenamePrefix)) {
    fprintf(stderr, "Error: The directory of Output Prefix doesn't exists.\n");
    exit(1);
  }

  for (int i = 0; i < n_inputFiles; i++) {
    if (!fileExists(inputFilenames[i])) {
      fprintf(stderr, "Error: The BAM/SAM file %s doesn't exists.\n", inputFilenames[i]);
      exit(1);
    }
  }
}


void
loadTranslateAndAnnotations(const char *indexPrefix, unsigned *dnaLength, Annotation **annotations, unsigned *numOfSeq,
                            unsigned int **ambMap, Translate **translate) {
  printf("Loading reference...\n");
  fflush(stdout);

  char *traFilename = (char *) malloc(strlen(indexPrefix) + 5);
  sprintf(traFilename, "%s.tra", indexPrefix);

  loadTranslate(traFilename, *dnaLength, ambMap, translate);

  free(traFilename);

  char *annFilename = (char *) malloc(strlen(indexPrefix) + 5);
  sprintf(annFilename, "%s.ann", indexPrefix);
  unsigned dnaLength2 = 0;

  loadSeqInfo(annFilename, dnaLength2, annotations, NULL, *numOfSeq);

  free(annFilename);

  printf("Done.\n");
  fflush(stdout);

  if (*dnaLength != dnaLength2) {
    fprintf(stderr,
            "Reference corrupted. Please use another one.\n");
    exit(1);
  }
}

unsigned int setExomeRegion(ExomeRegion *region, unsigned int textLength, unsigned int *bamshotBitVec) {
  unsigned int regionSize = 0;
  for (unsigned int i = 0; i < textLength; ++i) {
    if (bamshotBitVec[i >> 5] & (1 << (31 - (i & 31)))) {
      region[regionSize].startPos = i;
      region[regionSize].endPos = i;
      ++regionSize;
    }
  }
  return regionSize;
}

bool
isLowPopulation(SnpCounter *snpCounter, unsigned int textLength, unsigned int *bamshotBitVec, unsigned int &bitSetCount,
                unsigned int &maxArraySize) {
  bitSetCount = 0;
  maxArraySize = 0;
  for (unsigned int i = 0; i < textLength; ++i) {
    unsigned short w1 = snpCounter[i].weightedCount1;
    unsigned short w2 = snpCounter[i].weightedCount2;
    if (((w1 & w2) & 0xC000) == 0xC000) {
      bamshotBitVec[i >> 5] |= (1 << (31 - (i & 31)));
      bitSetCount++;
      maxArraySize = ((w1 & 0xFFF) < maxArraySize) ? maxArraySize : w1 & 0xFFF;
      maxArraySize = ((w2 & 0xFFF) < maxArraySize) ? maxArraySize : w2 & 0xFFF;
    } else {
      unsigned int strandCount = snpCounter[i].data.strandCounts.posStrandCount1 +
                                 snpCounter[i].data.strandCounts.negStrandCount1 +
                                 snpCounter[i].data.strandCounts.posStrandCount2 +
                                 snpCounter[i].data.strandCounts.negStrandCount2;
#           ifdef ENABLE_SOFTCLIP_COUNTER
      unsigned int scCount = snpCounter[i].softClipCount;
#           else
#           endif
      if (snpCounter[i].weightedCount1 != 0 || snpCounter[i].weightedCount2 != 0 || strandCount != 0 ||
          false) {
        bamshotBitVec[i >> 5] |= (1 << (31 - (i & 31)));
        bitSetCount++;
      }
    }
  }

  if ((double) bitSetCount / textLength < 0.1) {
    return true;
  }

  return false;
}

void writeSnapshotToFiles(SnpBundle *snpBundle, const char *outputFilenamePrefix, unsigned int *bamshotBitVec) {
  int prefix_len = strlen(outputFilenamePrefix);

  unsigned int bitSetCount;
  unsigned int maxArraySize;
  unsigned char isExome = 0;
  ExomeRegion *region = NULL;
  if (isLowPopulation(snpBundle->snpCounter, snpBundle->textLength, bamshotBitVec, bitSetCount, maxArraySize)) {
    region = (ExomeRegion *) malloc(bitSetCount * sizeof(ExomeRegion));
    setExomeRegion(region, snpBundle->textLength, bamshotBitVec);
    isExome = 1;

    snpBundle->snpMemoryPool->capacity = snpBundle->snpMemoryPool->curPtr + maxArraySize;
  }

  char *snapshotFilename = (char *) malloc(prefix_len + 20);
  sprintf(snapshotFilename, "%s.snapshot", outputFilenamePrefix);

  writeSnpInfoSnapshot(snpBundle, region, bitSetCount, isExome, snapshotFilename);

  if (region != NULL) {
    free(region);
  }

  free(snapshotFilename);
}

int main(int argn, char *args[]) {

  time_t startTime = time(NULL);

  
  char *indexPrefix;
  char *outputFilenamePrefix;
  char *inputFilenames[Max_Num_Of_Input_Files];
  int inputFileTypes[Max_Num_Of_Input_Files];
  int n_inputFiles;
  char *exomeFilename = NULL;

  parseCommandlineArguments(argn, args, indexPrefix, outputFilenamePrefix, inputFilenames, inputFileTypes,
                            n_inputFiles, exomeFilename);

  
  Annotation *annotations;
  unsigned numOfSeq = 0, dnaLength = 0;

  unsigned int *ambMap;
  Translate *translate;

  loadTranslateAndAnnotations(indexPrefix, &dnaLength, &annotations, &numOfSeq, &ambMap, &translate);

  printf("%d chromosomes, %ubp in length.\n", numOfSeq, dnaLength);
  fflush(stdout);


  
  unsigned numOfCPUThreads = 0;
  unsigned char recalScore[Q_SCORE_ARRAY_SIZE] = {0};
  unsigned char weightMap[Q_SCORE_ARRAY_SIZE] = {0};
  unsigned snpTrimSize = 0;
  unsigned softClipThreshold = 0;
  unsigned mapqThreshold = 0;
  unsigned indelWeightThreshold = 0;

  char *iniFile = (char *) malloc(strlen(args[0]) + 40);
  sprintf(iniFile, "%s.ini", args[0]);

  loadBam2snapshotParas(iniFile, &numOfCPUThreads, recalScore, weightMap, &snpTrimSize, &softClipThreshold,
                        &mapqThreshold, &indelWeightThreshold);

  free(iniFile);

  printf("Parameters:\n");
  printf("    #CPUThreads=%d\n"
             "    TrimSize=%d\n"
             "    SoftClipThreshold=%d\n"
             "    MQThreshold=%d\n"
             "    indelWeightThreshold=%d\n",
         numOfCPUThreads, snpTrimSize, softClipThreshold, mapqThreshold, indelWeightThreshold);

  
  SnpBundle snpBundle;
  snpBundle.textLength = dnaLength;

  printf("Allocating memory...\n");
  fflush(stdout);

  clearSnapshot(dnaLength, numOfCPUThreads, &snpBundle);

  printf("Done.\n");
  fflush(stdout);

  printTimeElapsed(startTime);


  
  AlignmentQC alignmentQC;
  resetAlignmentQC(alignmentQC);
  unsigned int *exomePosBitVec = NULL;
  if (exomeFilename != NULL) {
    printf("Reading Exome Index File %s...\n", exomeFilename);
    fflush(stdout);
    unsigned int numExomeRegion = 0;
    FILE *exomeFile = fopen(exomeFilename, "rb");
    fread(&numExomeRegion, sizeof(unsigned int), 1, exomeFile);
    exomePosBitVec = (unsigned int *) malloc(((dnaLength + 31) >> 5) * sizeof(unsigned int));
    memset(exomePosBitVec, 0, ((dnaLength + 31) >> 5) * sizeof(unsigned int));
    for (int i = 0; i < numExomeRegion; ++i) {
      unsigned int startPos, endPos;
      fread(&startPos, sizeof(unsigned int), 1, exomeFile);
      fread(&endPos, sizeof(unsigned int), 1, exomeFile);
      for (unsigned int j = startPos + 1; j <= endPos; ++j) {
        exomePosBitVec[j >> 5] |= (1 << (31 - (j & 31)));
      }
    }
    fclose(exomeFile);
    printf("Done.\n");
    fflush(stdout);
  }

  bamFilesToSnapshot(&snpBundle, dnaLength,
                     inputFilenames, inputFileTypes, n_inputFiles, outputFilenamePrefix,
                     recalScore, weightMap, snpTrimSize, softClipThreshold,
                     mapqThreshold, indelWeightThreshold,
                     annotations, numOfSeq, ambMap, translate,
                     alignmentQC, exomePosBitVec,
                     startTime);


  printTimeElapsed(startTime);

  
  printf("Writing snapshot...\n");
  fflush(stdout);

  unsigned int *bamshotBitVec = (unsigned int *) malloc(((snpBundle.textLength + 31) >> 5) * sizeof(unsigned int));
  writeSnapshotToFiles(&snpBundle, outputFilenamePrefix, bamshotBitVec);

  writeAlignmentQCToFile(alignmentQC, outputFilenamePrefix);

  printf("Done.\n");
  fflush(stdout);

  printTimeElapsed(startTime);

  
  freeSnapshot(&snpBundle);
  free(bamshotBitVec);
  free(exomePosBitVec);
  printTimeElapsed(startTime);

  return 0;
}


