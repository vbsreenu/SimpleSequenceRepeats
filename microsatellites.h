/*
Program for finding simple sequence repeats in a nucleotide fasta file.

Sreenu Vattipally
MRC-University of Glasgow Centre for Virus Research
University of Glasgow
Glasgow, G61 1QH

sreenu.vattipally@glasgow.ac.uk

*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>

#define MAX_SSR_LEN 50

// Function to check internal repeats in a microsatellite motif
int checkIntRepeats(char *repeat){
	int i=0, j=0, noIntRepeats=0;
	char leftMotif[MAX_SSR_LEN], rightMotif[MAX_SSR_LEN];


	for(i=1;i<strlen(repeat);i++){
		if(repeat[0]!=repeat[i]) {
			noIntRepeats=1;
			break;
		}
	}
	if(noIntRepeats==0){ return(1); exit(0); }


	if(strlen(repeat)==2){ if (repeat[0]==repeat[1]) return(1); else return(0); exit(0); }


	if(strlen(repeat)==3){ return(0); exit(0); }
	
	// SSR motif 3 cannot have internal repeats so, starting from 4

	if(strlen(repeat)>3){
		for(j=2;j<=(int)(strlen(repeat)/2); j++){
			if((strlen(repeat))%j==0){
				noIntRepeats=0;

				for (i=0;i<(strlen(repeat)-j);i+=j){

					leftMotif[0]='\0';
					rightMotif[0]='\0';

					strncpy(leftMotif, repeat+i,j);
					strncpy(rightMotif, repeat+(i+j),j);

					leftMotif[j]='\0';
					rightMotif[j]='\0';

        				if (strcmp (leftMotif,rightMotif)!=0) noIntRepeats=1;
				}
			}
			else noIntRepeats=1;

			if(noIntRepeats==0){
				return(1);
				break;
			}
		}
		if(noIntRepeats) return(0);
		else return(1);
	}
}
/* ======================================================= */ 

// Check mono nucleotide repeats
void checkMonoRepeats(char *genome, int minRepeatNumber, int flankSeq){


int i=0, ssrRepeats=0, ssrLocation=0, j=0;

for(i=0;i<strlen(genome);i++){

        if(genome[i]==genome[i+1]){
                ssrRepeats++;
                if(ssrRepeats==1) ssrLocation=i+1;
        }

       if(genome[i]!=genome[i+1] && ssrRepeats >= minRepeatNumber-1){
                printf("%c\t1\t%d\t%d\t%d",genome[i],ssrRepeats+1,ssrLocation, i+1);

		if(flankSeq){
				printf("\t");

				// printing left flanking sequence
				for(j=(ssrLocation-flankSeq)-1;j<ssrLocation-1; j++) 
					if(j>=0 && j<=strlen(genome)) printf("%c",genome[j]);

				printf("\t");

				// printing right flanking sequence
				for(j=i+1; j< (i+1)+flankSeq; j++) 
				if(j>0 && j<=strlen(genome)) printf("%c",genome[j]);
		}
		printf("\n");
                ssrRepeats=0;
        }

	//printf("minRepeatNumber\t%d\tssrRepeats %d\n",minRepeatNumber, ssrRepeats);

       if(genome[i]!=genome[i+1]) ssrRepeats=0;

	}
}
// End of checkMonoRepeats function

/* ======================================================= */ 

void checkRepeats (char *genome, int minMotifSize, int maxMotifSize, int minRepeatNumber, int flankSeq){

int i=0, ssrFlag=0,ssrRepeats=0, ssrLocation=0, motifLen, genomeLen, ssrEnd=0, j=0;
char leftMotif[MAX_SSR_LEN], rightMotif[MAX_SSR_LEN];

genomeLen=strlen(genome);

// Change the minimum motif size length if it set to 1 (this function is for non-mononucleotide repeats)
if(minMotifSize==1) minMotifSize=2;

for(motifLen=minMotifSize;motifLen<=maxMotifSize;motifLen++){  // 1st for loop

ssrRepeats=0;

        for (i=0;;i++){  // 2nd for loop
                // break the loop if reached the end
                if((genomeLen-i) < (motifLen*2)) {
                        if(ssrFlag!=0 && checkIntRepeats(leftMotif)==0)
				if(ssrRepeats >= minRepeatNumber-1)
                                printf("%s\t%d\t%d\t%d\n",leftMotif,ssrRepeats+1,ssrLocation+1,(motifLen*(ssrRepeats+1))+ssrLocation);
                        break;
                }

                leftMotif[0]='\0';
                rightMotif[0]='\0';
                strncpy(leftMotif, genome+i,motifLen);
                strncpy(rightMotif, genome+(i+motifLen),motifLen);

                leftMotif[motifLen]='\0';
                rightMotif[motifLen]='\0';


                if (strcmp (leftMotif,rightMotif)==0){
                        if (ssrFlag==0) ssrLocation=i;

                        i+=(motifLen-1);
                        ssrRepeats+=1;
                        ssrFlag++;
                }
                else{
                        if(ssrFlag!=0 && checkIntRepeats(leftMotif)==0){
				if(ssrRepeats >= minRepeatNumber-1){

				ssrEnd=(motifLen*(ssrRepeats+1))+ssrLocation;

                                printf("%s\t%d\t%d\t%d\t%d",leftMotif,motifLen,ssrRepeats+1,ssrLocation+1,ssrEnd);

					if(flankSeq){
						printf("\t");

					// printing left flanking sequence
					for(j=(ssrLocation-flankSeq);j<ssrLocation; j++) 
						if(j>=0 && j<=strlen(genome)) printf("%c",genome[j]);

					printf("\t");

					// printing right flanking sequence
					for(j=ssrEnd; j<(ssrEnd)+flankSeq; j++) 
						if(j>0 && j<=strlen(genome)) printf("%c",genome[j]);
					}
					printf("\n");
                        	}
			}
                        ssrFlag=0;
                        ssrRepeats=0;
                }
        }// End of 2nd for
} // End of 1st for
}
