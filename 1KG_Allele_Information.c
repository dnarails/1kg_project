/*
Developer:	Yu-Cheng Li
Description:	Calculate the AF, AC, and some information.
		For 1000 Genome Project PhaseI version.
Date:	2013.12.01
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define	LINE_LEN	30000

void calculate_allele_information (int start, int end,int table_count[14][3], int table_total[14]);

int main(int argc, char *argv[]){
	FILE	*file_csv;
	FILE	*file_table;
	char	line[LINE_LEN];
	int	i=0,j=0;
	char	subset[14][4];
	char	sample[1092][8];
	int	sample_index[1092];
	int	table_index[1092];
	int	table_count[14][3];
	int	table_total[14];

	char	*tok;

	int	subset_num=0;
	int	sample_num=0;
	int	index;
	int	ac_check;
	int	ac_subset;
	int	ac_set;

	int	count[3];
	int	total;

	char	chr[10];
	int	pos;
	char	rs[20];
	char	ref;
	char	alt;
	char	vt[10];
	int	ac;
	int	an;
	float	af;
	char	asn[10];
	char	amr[10];
	char	afr[10];
	char	eur[10];
	char	rsq[10];
	char	ldaf[10];
	char	genotype[3];

	/* Sample's Group */
	file_csv = fopen(argv[1],"r");
	/* Sample's Genotype */
	file_table = fopen(argv[2],"r");



	/* Read Sample's Group */
	while(fgets(line, LINE_LEN, file_csv) != NULL){
		tok = NULL;
		tok=strtok(line,",");
		index=0;
		while (tok !=NULL){
			if (index == 0){
				strncpy(subset[subset_num],tok,3);
				subset_num++;
				index++;
			}else if (strlen(tok)>2){
				strncpy(sample[sample_num],tok,7);
				sample_index[sample_num]=subset_num-1;
				sample_num++;
				index++;
			}
			tok = strtok(NULL, ",");
		}
		table_total[subset_num-1]=index-1;
	}
	



	/* Basic Column's Name */	
	for (i=0;i < 15;i++){
		fscanf(file_table,"%s", line);
		printf("%s\t",line);
	}

	/* Sample's Name */	
	for (i=0;i < 1092;i++){
		fscanf(file_table,"%s", rs);
		j=0;
		while (strcmp(rs+4, sample[j])!=0){
			j++;
		}
		table_index[i]=sample_index[j];
		//printf("%s\t%s\n",rs+4,subset[sample_index[j]]);
		//printf("%s\t",rs);
	}

	printf("ASN_00\tASN_01\tASN_11\tASN_AF\tASN_COUNT\t");
	printf("AMR_00\tAMR_01\tAMR_11\tAMR_AF\tAMR_COUNT\t");
	printf("AFR_00\tAFR_01\tAFR_11\tAFR_AF\tAFR_COUNT\t");
	printf("EUR_00\tEUR_01\tEUR_11\tEUR_AF\tEUR_COUNT\t");

	for (i = 0;i < 14;i++){
		printf("%s_00\t%s_01\t%s_11\t%s_AF\t%s_COUNT\t",subset[i], subset[i], subset[i], subset[i], subset[i]);
	}
	printf("\n");
	index = 0;	


	/* Sample's Genotype */
	while(fscanf(file_table,"%s", chr)!=EOF){
		/* 1000G's information */
		fscanf(file_table,"%d", &pos);
		fscanf(file_table,"%s", rs);
		fscanf(file_table,"\t%c", &ref);
		fscanf(file_table,"\t%c", &alt);
		fscanf(file_table,"%s", vt);
		fscanf(file_table,"%d", &ac);
		fscanf(file_table,"%d", &an);
		fscanf(file_table,"%f", &af);
		fscanf(file_table,"%s", asn);
		fscanf(file_table,"%s", amr);
		fscanf(file_table,"%s", afr);
		fscanf(file_table,"%s", eur);
		fscanf(file_table,"%s", rsq);
		fscanf(file_table,"%s", ldaf);
		printf("%s\t%d\t%s\t%c\t%c\t%s\t%d\t%d\t%f\t%s\t%s\t%s\t%s\t%s\t%s", chr, pos, rs, ref, alt, vt, ac, an, af, asn, amr, afr, eur, rsq, ldaf);

		/* 1092 Sample's Genotype */
		for (i=0;i < 1092;i++){
			fscanf(file_table,"%s", genotype);
			if (genotype[0] == genotype[1]){
				if (genotype[0] == ref){
					table_count[table_index[i]][0]++;
				}else if (genotype[0] == alt){
					table_count[table_index[i]][2]++;
				}
			}else {
				table_count[table_index[i]][1]++;
			}
		}
		/* 4 main Group */
		calculate_allele_information ( 0, 2, table_count, table_total);
		calculate_allele_information ( 3, 5, table_count, table_total);
		calculate_allele_information ( 6, 8, table_count, table_total);
		calculate_allele_information ( 9, 14, table_count, table_total);

		/* 14 sub-group and double check the result*/
		ac_check=0;
		for (i = 0;i < 14;i++){
			ac_subset = (table_count[i][1]+table_count[i][2]+table_count[i][2]);
			ac_check += ac_subset;
			printf("\t%.5f\t%.5f\t%.5f\t%.5f", (float)table_count[i][0] /table_total[i], (float)table_count[i][1] /table_total[i], (float)table_count[i][2] /table_total[i], (float)ac_subset /(table_total[i]*2));
			printf("\t%d,%d,%d,%d", table_count[i][0], table_count[i][1], table_count[i][2], ac_subset);
			table_count[i][0]	= 0;
			table_count[i][1]	= 0;
			table_count[i][2]	= 0;
		}
		printf("\n");
		if (ac != ac_check){
			printf("Error\n");
			return -1;
		}
		index++;
	}

	fclose(file_csv);
	fclose(file_table);
	return 0;
}

void calculate_allele_information (int start, int end,int table_count[14][3], int table_total[14]){
	int	i;
	int	count[3];
	int	total;
	int	ac_set;

	count[0] = 0;	count[1] = 0;	count[2] = 0;	total = 0;
	for (i = start; i<=end;i++){
		count[0] += table_count[i][0];
		count[1] += table_count[i][1];
		count[2] += table_count[i][2];
		total	+= table_total[i];
	}
	ac_set = count[1]+count[2]+count[2];
	printf("\t%.5f\t%.5f\t%.5f\t%.5f", (float)count[0] /total, (float)count[1] /total, (float)count[2] /total, (float)ac_set /(total*2));
	printf("\t%d,%d,%d,%d", count[0], count[1], count[2], ac_set);
}
