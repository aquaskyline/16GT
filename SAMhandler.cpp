

#include "SAMhandler.h"

int SAMIUint8ConcatUint8 ( uint8_t * data, int * curSize,
                           uint8_t key )
{
	data[ ( *curSize ) ++] = key;
}
int SAMIUint8ConcatUint32 ( uint8_t * data, int * curSize,
                            uint32_t key )
{
	int i;
	int len = sizeof ( uint32_t ) / sizeof ( uint8_t );
	uint8_t * key_8 = ( uint8_t * ) &key;

	for ( i = 0; i < len; i++ )
	{
		data[ ( *curSize ) ++] = key_8[i];
	}
}
int SAMIUint8ConcatString ( uint8_t * data, int * curSize,
                            char * key, int len )
{
	int i;

	for ( i = 0; i < len; i++ )
	{
		data[ ( *curSize ) ++] = key[i];
	}
}

void SAMOutputHeaderConstruct(bam_header_t *sheader, HSP *hsp)
{
	int i, j;
	
	sheader->n_targets = hsp->numOfSeq;
	sheader->target_name = ( char ** ) malloc ( sizeof ( char * ) * hsp->numOfSeq );
	sheader->target_len = ( uint32_t * ) malloc ( sizeof ( uint32_t ) * hsp->numOfSeq );

	for ( i = 0; i < hsp->numOfSeq; i++ )
	{
		for ( j = 0; j < 255; j++ )
		{
			if ( hsp->annotation[i].text[j] == '\0' ||
			        hsp->annotation[i].text[j] == ' ' ||
			        hsp->annotation[i].text[j] == '\t' ||
			        hsp->annotation[i].text[j] == '\r' ||
			        hsp->annotation[i].text[j] == '\n' )
			{
				break;
			}
		}

		hsp->annotation[i].text[j] = '\0';
		
		
		sheader->target_name[i] = hsp->annotation[i].text;
		
		sheader->target_len[i] = hsp->seqActualOffset[i].endPos - hsp->seqActualOffset[i].startPos;
	}

	char * program_info = PROGRAM_INFO_SAM;
	int textLen = strlen ( hsp->readGroup ) + strlen ( hsp->sampleName ) + strlen ( hsp->readGrpOption ) + strlen ( program_info ) + 50;
	sheader->text = ( char * ) malloc ( textLen * sizeof ( char ) );
	sprintf ( sheader->text, "@RG\tID:%s\tSM:%s\tPL:ILLUMINA\t%s\n%s", hsp->readGroup, hsp->sampleName, hsp->readGrpOption, program_info );
	sheader->l_text = strlen ( sheader->text );
	
	
	sheader->hash = NULL;
	sheader->rg2lib = NULL;
}

void SAMOutputHeaderDestruct ( bam_header_t * sheader )
{
	free ( sheader->target_name );
	free ( sheader->target_len );
	free ( sheader->text );
}


void SAMOccurrenceConstruct ( OCC * occ )
{
	occ->SAMOutBuffer.data = ( uint8_t * ) malloc ( sizeof ( uint8_t ) * SAM_MDATA_SIZE );
}

void SAMOccurrenceDestruct ( OCC * occ )
{
	free ( occ->SAMOutBuffer.data );
}
