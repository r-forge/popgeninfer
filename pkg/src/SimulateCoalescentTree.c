#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>

#define COALESCENCE		0
#define MIGRATION		1
#define OUTPUT           "tmp.dat"


//	Tests:
//	10 demes sampled out of 10, with N = 1000, m = 0.001, mu = 0.0001, n = 100
//	OBS:	0.4660676	0.2451498
//	EXP:	0.465482	0.24519

//	10 demes sampled out of 10, with N = 1000, m = 0.01, mu = 0.0001, n = 100
//	OBS:	0.3510464	0.3221066 
//	EXP:	0.350394	0.321906

//	20 demes sampled out of 20, with N = 200, m = 0.01, mu = 0.0001, n = 100
//	OBS:	0.5971124	0.5022783
//	EXP:	0.595106	0.501354


struct Node {
	double time;                                                                      // age of the node (time = 0 at sample time)
	int allele;                                                                    // allelic state of the node
	struct Node *ancestor;                                                         // address of the ancestor of the node
	struct Node *descendant[2];                                              // addresses of the descendants of the nodes
};

struct Tree {
	int n;
	int nbr_ancestors;
	double N;
	double m;
	struct Node *tree;
	struct Node **list;
};

typedef struct {
	double mu;
	double N;
	int n;
	int nt;
	int nd;
	double m;
	int nloci;
} Parameters;

typedef struct {
	int nbr;
	int *lineage;
} Deme;

void AddMutations(Parameters P,struct Node *node);
struct Node *BuildTreeHudson(Parameters P);
#ifdef TEST
void ComputeIdentityProbabilities(Parameters P);
#endif
int Mutation(int allele);
void WriteOutputs(Parameters P);

struct Tree *population;
FILE *outputfile;

void SimulateCoalescentTree(int *total_number_of_demes,int *effective_size,double *migration_rate,double *mutation_rate,int *number_of_sampled_demes,int *sample_size,int *number_of_loci) {

	int i,j,k;
	int sampled_lineages,max_lineages;
	Parameters P;  
	struct Node *mrca;
	
	GetRNGstate();	
	outputfile = fopen(OUTPUT, "w");
	P.nt = *total_number_of_demes;
	P.mu = *mutation_rate;
	P.N = *effective_size;
	P.n = *sample_size;
	P.nd = *number_of_sampled_demes;
	P.m = *migration_rate;
	P.nloci = *number_of_loci;
	population = (struct Tree *) malloc (P.nt * sizeof (struct Tree));
	sampled_lineages = P.n * P.nd;
	for (i = 0; i < P.nt; ++i) {
		population[i].N = P.N;
		population[i].m = P.m;
		if (i < P.nd) {
			population[i].n = P.n;
		} else {
			population[i].n = 0;			
		}
	}	
	for (i = 0; i < P.nloci; ++i) {                                                // Loop over loci
		max_lineages = 2 * sampled_lineages - 1;
		for (j = 0; j < P.nt; ++j) {                                           // Initialization of terminal nodes in population 0
			population[j].list = (struct Node **) malloc (max_lineages * sizeof (struct Node *));// These are the pointers to the lineages left in the population
			population[j].tree =(struct Node *) malloc(max_lineages * sizeof(struct Node));      // This is the genealogy
		}
		for (j = 0; j < P.nt; ++j) {                                           // Initialization of terminal nodes in population 0
			for (k = 0; k < max_lineages; ++k) {                                           // Initialization of terminal nodes in population 0
				(population[j].tree[k]).ancestor = NULL;
				(population[j].tree[k]).descendant[0] = NULL;
				(population[j].tree[k]).descendant[1] = NULL;
			}
		}
		mrca = BuildTreeHudson(P);                                                            // Give a sample of genes from a (generation-by-generation) coalescent alogrithm with mutations
		mrca -> allele = 0;
		AddMutations(P,mrca);
#ifndef TEST
		WriteOutputs(P);
#else                                                                                                               // ... only if in 'non-test' mode
		ComputeIdentityProbabilities(P);
#endif
		for (j = 0; j < P.nt; ++j) {                                           // Initialization of terminal nodes in population 0
			free(population[j].list);
			free(population[j].tree);
		}
	}
	free(population);
	fclose(outputfile);                                                            // Close output files
	PutRNGstate();	
}

struct Node *BuildTreeHudson(Parameters P)

{
	int *last_node;																				// last_node[i] is the index of the last created node in population i (former nbr_nodes[i])
	int *nbr_lineages;																			// nbr_lineages[i] is the number of surviving lineages in population i
	int total_nbr_lineages = 0;																	// total_nbr_lineages is the total number of survivng lienages across all demes
	double *p;
	struct Node *topnode = NULL;
	int i,j,number1,number2,event;
	double sum,x;
	double time = 0.0;
	
	for (i = 0; i < P.nt; ++i) {
		for (j = 0; j < population[i].n; ++j) {
			population[i].tree[j].time = time;														// tree[j] for j = 0, 1, ..., (n - 1) are the sampled (terminal) nodes
			population[i].list[j] = population[i].tree + j;													// list points to the sampled nodes of population i
		}
	}
	last_node = (int *) malloc(P.nt * sizeof(int));
	nbr_lineages = (int *) malloc(P.nt * sizeof(int));
	p = (double *) malloc((2 * P.nt) * sizeof(double));
	for (i = 0; i < P.nt; ++i) {
		last_node[i] = nbr_lineages[i] = population[i].n;   // !!!!!
		total_nbr_lineages += nbr_lineages[i];
	}
	while (total_nbr_lineages > 1) {
		p[0] = (double) nbr_lineages[0] * (nbr_lineages[0] - 1) / 2 / population[0].N;					// Probability of coalescence in population 0
		p[1] = (double) p[0] + nbr_lineages[0] * population[0].m;												// Probability of migration in population 0
		for (i = 2; i < (2 * P.nt); i += 2) {
			j = (int) i / 2;																		// Index for the demes
			p[i] = p[i - 1] + (double) nbr_lineages[j] * (nbr_lineages[j] - 1) / 2 / population[j].N;	// Probability of coalescence in population i
			p[i + 1] = p[i] + (double) nbr_lineages[j] * population[j].m;										// Probability of migration in population i
		}
		sum = p[2 * P.nt - 1];																// Sum over all probabilities
		for (i = 0; i < (2 * P.nt); ++i) {
			p[i] /= sum;
		}
		p[(2 * P.nt - 1)] = 1.0;
		time += - log(1 - unif_rand()) / sum;												// Time of the next event		
		x = unif_rand();
		i = 0;
		while (x > p[i]) i += 1;		
		event = i%2;
		i /= 2;
		switch (event) {
			case COALESCENCE : {															// Coalescence in population i
				number1 = (int) (unif_rand() * nbr_lineages[i]);						// Draw two distinct lineages at random (number1 and number2) in population[i]
				do {
					number2 = (int) (unif_rand() * nbr_lineages[i]);						// Draw two distinct lineages at random (number1 and number2) in population[i]
				}
				while (number2 == number1);
				population[i].tree[last_node[i]].time = time;									// population[i].tree[last_node[i]] is the next internal node
				population[i].list[number1] -> ancestor = population[i].tree + last_node[i];		// The ancestor of population[i].tree[number1] is population[i].tree[last_node[i]] (the internal node now considered)
				population[i].list[number2] -> ancestor = population[i].tree + last_node[i];		// The ancestor of population[i].tree[number2] is population[i].tree[last_node[i]] (the internal node now considered)
				population[i].tree[last_node[i]].descendant[0] = population[i].list[number1];		// The descendant_1 of population[i].tree[last_node[i]] (the internal node now considered) is population[i].list[number1]
				population[i].tree[last_node[i]].descendant[1] = population[i].list[number2];		// The descendant_2 of population[i].tree[last_node[i]] (the internal node now considered) is population[i].list[number2]
				population[i].list[number1] = population[i].tree + last_node[i];					// population[i].list[number1] now points to the next internal node
				population[i].list[number2] = population[i].list[--nbr_lineages[i]];				// The number of lineages in population[i] (nbr_lineages[i]) decreases (coalescence) AND population[i].list[number2] points to the last active lineage in population[i]
				--total_nbr_lineages;														// The total number of lineages decreases (coalescence)
				++last_node[i];																// The index (last_node[i]) of the next internal node to create is increased
				break;					
			}
			case MIGRATION : {																// Migration from population j in population i (backward in time)
				do {
					j = (int) (unif_rand() * P.nt);
				}
				while (j == i);
				number1 = (int) (unif_rand() * nbr_lineages[i]);						// Draw one lineage at random (number1) in population[i]
				population[j].list[nbr_lineages[j]++] = population[i].list[number1];				//
				population[i].list[number1] = population[i].list[--nbr_lineages[i]];				//
				break;
			}
		}
	}
	for (i = 0; i < P.nt; ++i) {
		population[i].nbr_ancestors = nbr_lineages[i];
	}
	if (total_nbr_lineages == 1) {															    		    	    // If there remains a single lineage (e.g., if this is the mrca of the metapopulation considered...
		for (i = 0; i < P.nt; ++i) {															    // ... or, more insterestingly, if this is the mrca of the full sample)
			if (nbr_lineages[i] > 0) break;															    		    // then find the population where this lineage is
		}
		topnode = &population[i].tree[(last_node[i] - 1)];															// 'topnode' is the address of this last lineage 
	}
	free(last_node);
	free(nbr_lineages);
	free(p);
	return(topnode);															    		    	                // 'topnode' is passed to the function MakeSamples 
}

void AddMutations(Parameters P,
                  struct Node *node)

{
	int i,j;
	int nbr_mut,tmp;
	int time;
	
	for (i = 0; i < 2; ++i) {                                // This loop is over the total number of descendants of the current node
		if (node -> descendant[i] != NULL) {
			(node -> descendant[i]) -> allele = node -> allele;                        // Copy the ancestral allelic state into the current node
			time = node -> time - ((node -> descendant[i]) -> time);                   // Calculate the time elapsed since ancestor
			nbr_mut = (int) rpois(time * P.mu);                                           // Calculate the number of mutations along the branch
			for (j = 0; j < nbr_mut; ++j) {                                            // This loop is over all the mutations that occurred
				do {                                                                     // Draw a new allelic state
					tmp = Mutation((node -> descendant[i]) -> allele);                                                      // bug fixed 17-01-2007 : was 'tmp = Mutation(node -> allele)'
				}                                                                        // which caused the number of mutations to be underestimated... (mutations different from the ancestral state, not from the
				while (tmp == (node -> descendant[i]) -> allele);                        // the new allele must be different from the current state
				(node -> descendant[i]) -> allele = tmp;                                 // Copy the mutation into the current node
			}
			AddMutations(P,node -> descendant[i]);                                // Play it again (recursively, throughout the tree) starting from the descendants of the current node
		}
	}
}

int Mutation(int allele)

{
	int new;
 
	new = 1 - allele;
	return(new);
}

void WriteOutputs(Parameters P)

{
	int i,j,cpt;
	
	for (i = 0; i < P.nd; ++i) {
		cpt = 0;
		for (j = 0; j < population[i].n; ++j) {
			if (population[i].tree[j].allele == 0) ++cpt;
		}
		fprintf(outputfile,"%8d\t",cpt);
	}
	fprintf(outputfile,"\n");
}

#ifdef TEST
void ComputeIdentityProbabilities(Parameters P)

{
	int i,j,k,l;
	int npairs[2];
	double Q[2];
	
	npairs[0] = npairs[1] = 0;
	Q[0] = Q[1] = 0.0;
	for (i = 0; i < P.nd; ++i) {
		for (j = 0; j < (population[i].n - 1); ++j) {
			for (k = (j + 1); k < population[i].n; ++k) {
				if (population[i].tree[j].allele == population[i].tree[k].allele) Q[0] += 1;
				npairs[0] += 1;
			}
		}
	}
	for (i = 0; i < (P.nd - 1); ++i) {
		for (j = (i + 1); j < P.nd; ++j) {
			for (k = 0; k < population[i].n; ++k) {
				for (l = 0; l < population[j].n; ++l) {
					if (population[i].tree[k].allele == population[j].tree[l].allele) Q[1] += 1;
					npairs[1] += 1;
				}
			}
		}
	}
	for (i = 0; i < 2; ++i) {
		Q[i] /= (double) npairs[i];
	}
	fprintf(outputfile,"%10.6f %10.6f\n",Q[0],Q[1]);
}
#endif
