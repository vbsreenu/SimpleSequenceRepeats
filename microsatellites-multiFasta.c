/*
Program for finding simple sequence repeats in a nucleotide fasta file.

Sreenu Vattipally
MRC-University of Glasgow Centre for Virus Research
University of Glasgow
Glasgow, G61 1QH

sreenu.vattipally@glasgow.ac.uk
*/

#include "microsatellites.h"

int main(int argc, char **argv) {

	FILE *input;
	char c, *genome;
	int i=0, maxGenomeSize=100000, flag=1, genomeLen=0;
	int minMotifSize=2, maxMotifSize=50, minRepeatNumber=2, flankSeq=0; 
	char *inFile = NULL;

	while ((c = getopt (argc, argv, "m:M:r:i:f:h")) != -1) 
	switch (c) { 
	case 'm':
        minMotifSize = atoi(optarg);
        break;

      	case 'M':
        maxMotifSize = atoi(optarg);
	if(maxMotifSize < minMotifSize){
		printf("\n\nWrong -M option. Maximum motif (-M) size should not be less than minimum motif (-m) size\n\n");
		exit(0);
	}
        break;

      	case 'r':
	minRepeatNumber = atoi(optarg);
        break;

      	case 'f':
	flankSeq = atoi(optarg); 
        break;

      	case 'i':
	inFile = optarg;
        break;

      	case 'h':
		printf("Usage: SimpleSequenceRepeats -i inFile -m minMotifSize -M maxMotifSize -r minRepeatNumber -f flankingSequence\n");
        exit (0);

      	default:
		printf("Usage: SimpleSequenceRepeats -i inFile -m minMotifSize -M maxMotifSize -r minRepeatNumber -f flankingSequence\n");
        abort ();
      }


// printf ("%d\t%d\t%d\t%s\t%d\n",minMotifSize, maxMotifSize, minRepeatNumber, inFile, flankSeq);


	// Read the genome(s) file
	if((input = fopen(inFile,"r"))!=NULL){
	genome=(char*) calloc(maxGenomeSize,sizeof(char));


	while((c=fgetc(input))){
		if(feof(input)) {
			// End of the sequence
			genome[genomeLen]='\0';
			if(minMotifSize == 1) checkMonoRepeats(genome, minRepeatNumber, flankSeq);
			//if(minMotifSize > 1) 
			checkRepeats(genome, minMotifSize, maxMotifSize, minRepeatNumber, flankSeq);
			break;
		}

		if(c=='>') {
			if(strlen(genome)){
				genome[genomeLen]='\0';
				if(minMotifSize == 1) checkMonoRepeats(genome, minRepeatNumber, flankSeq);
				//if(minMotifSize > 1) 
				checkRepeats(genome, minMotifSize, maxMotifSize, minRepeatNumber,flankSeq);
			}
			genome[0]='\0';
			flag=genomeLen=0;
		}

		if(flag==0 && c!='>') printf("%c",c);
		if(c=='\n') flag=1;
		if(flag==1 && c!='\n') genome[genomeLen++]=c;

		// Reallocate the memory if needed
		if(genomeLen > maxGenomeSize-10){
			maxGenomeSize+=10000;
			genome=(char*) realloc(genome, maxGenomeSize*sizeof(char));
		}

	} // End of while
	fclose(input);
	}// End of fopen if

	else printf("Usage: SimpleSequenceRepeats  -i inFile -m minMotifSize -M maxMotifSize -r minRepeatNumber -f flankingSequence\n");
	
} // End of Main
