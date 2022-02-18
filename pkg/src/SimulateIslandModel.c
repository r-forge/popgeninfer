#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>

#define WITHIN               0
#define BETWEEN              1
#define COALESCENCE          0
#define MIGRATION            1
#define IAM                  0
#define KAM                  1
#define SMM                  2

#define OUTPUT               "sim_stats.out"
#define PRIOR                "prior_val.in"

struct node_struct {
    int allele;                                       // allelic state of the node
  double time;                                      // age of the node (time = 0 at sample time)
    struct node_struct *ancestor;                     // address of the ancestor of the node
    struct node_struct *descendant[2];                // addresses of the descendants of the nodes
};

struct tree_struct {
    int nbr_sampled_lineages;
    int nbr_ancestors;
  double M;
    struct node_struct *tree;
    struct node_struct **list;
};

typedef struct {
  int mutation_model;
  int nbr_allelic_states;
    int nbr_demes;
  int nbr_loci;
    int nbr_sampled_demes;
  int nbr_simulations;
  int sample_size;
  double M;         // = 2Nm
  double theta;     // = 2Nu
} parameters_struct;

typedef struct {
  int min;
  int max;
  double total_nbr_alleles;
  double avrge_nbr_alleles;
  double total_heterozygosity;
  double avrge_heterozygosity;
  double total_shannon;
  double avrge_shannon;
  double total_var_allele_length;
  double avrge_var_allele_length;
  double F_ST;
  double F_ST_num;
  double F_ST_den;
#ifdef DEBUG
  double Q[2];
#endif
} statistics_struct;

struct node_struct *build_tree_Hudson(parameters_struct P);
void add_mutations(parameters_struct P,struct node_struct *node);
int mutation(parameters_struct P,int allele);
void write_outputs(statistics_struct multi_locus);
void compute_locus_specific_statistics(parameters_struct P,statistics_struct *sum_stat);
void compute_multilocus_statistics(parameters_struct P,statistics_struct *multi_locus,statistics_struct *locus_specific);

#ifdef DEBUG
void compute_locus_specific_identity_probabilities(parameters_struct P,statistics_struct *sum_stat);
void compute_multilocus_identity_probabilities(parameters_struct P,statistics_struct *multi_locus,statistics_struct *locus_specific);
void write_identity_probabilities(statistics_struct multi_locus);
#endif

int k_iam;
struct tree_struct *deme;
FILE *outfile,*priorfile;

void SimulateIslandModel(int *number_of_simulations,
                         int *mutation_model,
                         int *total_number_of_demes,
                         int *number_of_loci,
                         int *number_of_sampled_demes,
                         int *sample_size) {

    int i,j,k,sim;
    int sampled_lineages,max_lineages;
  int dummy,ncols,nlines;
  int error;
  char X;
  double val;
    parameters_struct P;
  statistics_struct *locus_specific;
  statistics_struct multi_locus;
    struct node_struct *mrca;
    
    GetRNGstate();
  P.mutation_model = *mutation_model;
  P.nbr_demes = *total_number_of_demes;
    P.nbr_loci = *number_of_loci;
    P.nbr_sampled_demes = *number_of_sampled_demes;
  P.nbr_simulations = *number_of_simulations;
    P.sample_size = *sample_size;
  if ((priorfile = fopen(PRIOR,"r")) == NULL) {
    Rprintf("%s: file not found... Run the generate.prior() function first...\n",PRIOR);
    goto end;
  }
  nlines = 0;
  do {
    ncols = 0;
    do {
      dummy = -1;
      if ((X = fscanf(priorfile,"%lf%n",&val,&dummy)) != EOF) {
        if (dummy == -1) {
          Rprintf("%s: unexpected character at line no. %d... Please check the prior file\n",PRIOR,(nlines + 1));
          goto end;
        }
      }
      ncols++;
      do {
        X = getc(priorfile);
        if (X == '\n') {
          nlines++;
          if (ncols > 2) {
            Rprintf("%s: the number of columns at line no. %d is larger than expected... Please check the prior file\n",PRIOR,nlines);
            goto end;
          }
          if (ncols < 2) {
            Rprintf("%s: the number of columns at line no. %d is lower than expected... Please check the prior file\n",PRIOR,nlines);
            goto end;
          }
          if (ncols == 2) {
            break;
          }
        }
      } while (X  == '\t' || X == ' ');
      fseek(priorfile, -1, SEEK_CUR);
    } while ((X != '\n') & (X != EOF));
  } while (X != EOF);
  if (nlines > P.nbr_simulations) {
    Rprintf("%s: the number of lines (%d) is larger than expected... Please check the prior file\n",PRIOR,nlines);
    goto end;
  }
  if (nlines < P.nbr_simulations) {
    Rprintf("%s: the number of lines (%d) is lower than expected... Please check the prior file\n",PRIOR,nlines);
    goto end;
  }
  if ((outfile = fopen(OUTPUT, "w")) == NULL) {
    Rprintf("Can't open %s for output!\n",OUTPUT);
    goto end;
  }
#ifndef DEBUG
  fprintf(outfile,"tot_nb_all  avg_nb_all     tot_het     avg_het        F_ST tot_var_sze avg_var_sze tot_Shannon avg_Shannon\n");
#else
  fprintf(outfile,"  Q_within   Q_between\n");
#endif
    sampled_lineages = P.sample_size * P.nbr_sampled_demes;
  max_lineages = 2 * sampled_lineages - 1;
  deme = (struct tree_struct *) malloc (P.nbr_demes * sizeof (struct tree_struct));
  for (i = 0; i < P.nbr_demes; ++i) {
        if (i < P.nbr_sampled_demes) {
            deme[i].nbr_sampled_lineages = P.sample_size;
        } else {
            deme[i].nbr_sampled_lineages = 0;
        }
    }
  for (j = 0; j < P.nbr_demes; ++j) {
    deme[j].list = (struct node_struct **) malloc (max_lineages * sizeof (struct node_struct *)); // These are the pointers to the lineages left in the deme
    deme[j].tree =(struct node_struct *) malloc(max_lineages * sizeof(struct node_struct)); // This is the genealogy
  }
  locus_specific = (statistics_struct *) malloc (P.nbr_loci * sizeof(statistics_struct));
  rewind(priorfile);
  for (sim = 0; sim < P.nbr_simulations; ++sim) {
      if ((error = fscanf(priorfile,"%lf",&P.theta)) == 0 || error == EOF) {
          Rprintf("STOPPED: problem reading prior file...\n");
      }
      if ((error = fscanf(priorfile,"%lf",&P.M)) == 0 || error == EOF) {
          Rprintf("STOPPED: problem reading prior file...\n");
      }
    for (j = 0; j < P.nbr_demes; ++j) {
      deme[j].M = P.M;
    }
    for (i = 0; i < P.nbr_loci; ++i) {
      for (j = 0; j < P.nbr_demes; ++j) {
        for (k = 0; k < max_lineages; ++k) {
          (deme[j].tree[k]).ancestor = NULL;
          (deme[j].tree[k]).descendant[0] = NULL;
          (deme[j].tree[k]).descendant[1] = NULL;
        }
      }
      mrca = build_tree_Hudson(P);
      mrca -> allele = k_iam = 0;
      add_mutations(P,mrca);
#ifndef DEBUG
      compute_locus_specific_statistics(P,&locus_specific[i]);
#else
      compute_locus_specific_identity_probabilities(P,&locus_specific[i]);
#endif
    }
#ifndef DEBUG
    compute_multilocus_statistics(P,&multi_locus,locus_specific);
    write_outputs(multi_locus);
#else
    compute_multilocus_identity_probabilities(P,&multi_locus,locus_specific);
    write_identity_probabilities(multi_locus);
#endif
  }
  for (j = 0; j < P.nbr_demes; ++j) {
    free(deme[j].list);
    free(deme[j].tree);
  }
    free(deme);
  free(locus_specific);
  fclose(outfile);
end :
  fclose(priorfile);
    PutRNGstate();
}

struct node_struct *build_tree_Hudson(parameters_struct P)

{
    int *last_node;                                                                        // last_node[i] is the index of the last created node in deme i (former nbr_nodes[i])
    int *nbr_lineages;                                                                // nbr_lineages[i] is the number of surviving lineages in deme i
    int total_nbr_lineages = 0;                                                // total_nbr_lineages is the total number of survivng lienages across all demes
    double *p;
    struct node_struct *topnode = NULL;
    int i,j,number1,number2,event;
    double sum,x;
    double time = 0.0;
    
    for (i = 0; i < P.nbr_demes; ++i) {
        for (j = 0; j < deme[i].nbr_sampled_lineages; ++j) {
            deme[i].tree[j].time = time;                                    // tree[j] for j = 0, 1, ..., (n - 1) are the sampled (terminal) nodes
            deme[i].list[j] = deme[i].tree + j;                        // list points to the sampled nodes of deme i
        }
    }
    last_node = (int *) malloc(P.nbr_demes * sizeof(int));
    nbr_lineages = (int *) malloc(P.nbr_demes * sizeof(int));
    p = (double *) malloc((2 * P.nbr_demes) * sizeof(double));
    for (i = 0; i < P.nbr_demes; ++i) {
        last_node[i] = nbr_lineages[i] = deme[i].nbr_sampled_lineages;
        total_nbr_lineages += nbr_lineages[i];
    }
    while (total_nbr_lineages > 1) {
        p[0] = (double) nbr_lineages[0] * (nbr_lineages[0] - 1)/* / 2 / deme[0].N*/; // Probability of coalescence in deme 0
        p[1] = (double) p[0] + nbr_lineages[0] * deme[0].M/*deme[0].m*/;    // Probability of migration in deme 0
        for (i = 2; i < (2 * P.nbr_demes); i += 2) {
            j = (int) i / 2;                                                            // Index for the demes
            p[i] = p[i - 1] + (double) nbr_lineages[j] * (nbr_lineages[j] - 1)/* / 2 / deme[j].N*/;    // Probability of coalescence in deme i
            p[i + 1] = p[i] + (double) nbr_lineages[j] * deme[j].M/*deme[j].m*/;    // Probability of migration in deme i
        }
        sum = p[2 * P.nbr_demes - 1];                                        // Sum over all probabilities
        for (i = 0; i < (2 * P.nbr_demes); ++i) {
            p[i] /= sum;
        }
        p[(2 * P.nbr_demes - 1)] = 1.0;
        time += - log(1 - unif_rand()) / sum;                        // Time of the next event
        x = unif_rand();
        i = 0;
        while (x > p[i]) i += 1;
        event = i%2;
        i /= 2;
        switch (event) {
            case COALESCENCE : {                                                    // Coalescence in population i
                number1 = (int) (unif_rand() * nbr_lineages[i]); // Draw two distinct lineages at random (number1 and number2) in deme[i]
                do {
                    number2 = (int) (unif_rand() * nbr_lineages[i]); // Draw two distinct lineages at random (number1 and number2) in deme[i]
                }
                while (number2 == number1);
                deme[i].tree[last_node[i]].time = time;            // deme[i].tree[last_node[i]] is the next internal node
                deme[i].list[number1] -> ancestor = deme[i].tree + last_node[i]; // The ancestor of deme[i].tree[number1] is deme[i].tree[last_node[i]] (the internal node now considered)
                deme[i].list[number2] -> ancestor = deme[i].tree + last_node[i]; // The ancestor of deme[i].tree[number2] is deme[i].tree[last_node[i]] (the internal node now considered)
                deme[i].tree[last_node[i]].descendant[0] = deme[i].list[number1]; // The descendant_1 of deme[i].tree[last_node[i]] (the internal node now considered) is deme[i].list[number1]
                deme[i].tree[last_node[i]].descendant[1] = deme[i].list[number2];    // The descendant_2 of deme[i].tree[last_node[i]] (the internal node now considered) is deme[i].list[number2]
                deme[i].list[number1] = deme[i].tree + last_node[i]; // deme[i].list[number1] now points to the next internal node
                deme[i].list[number2] = deme[i].list[--nbr_lineages[i]]; // The number of lineages in deme[i] (nbr_lineages[i]) decreases (coalescence) AND deme[i].list[number2] points to the last active lineage in deme[i]
                --total_nbr_lineages;                                                // The total number of lineages decreases (coalescence)
                ++last_node[i];                                                            // The index (last_node[i]) of the next internal node to create is increased
                break;
            }
            case MIGRATION : {                                                        // Migration from deme j in deme i (backward in time)
                do {
                    j = (int) (unif_rand() * P.nbr_demes);
                }
                while (j == i);
                number1 = (int) (unif_rand() * nbr_lineages[i]); // Draw one lineage at random (number1) in deme[i]
                deme[j].list[nbr_lineages[j]++] = deme[i].list[number1];
                deme[i].list[number1] = deme[i].list[--nbr_lineages[i]];
                break;
            }
        }
    }
    for (i = 0; i < P.nbr_demes; ++i) {
        deme[i].nbr_ancestors = nbr_lineages[i];
    }
    if (total_nbr_lineages == 1) {                                        // If there remains a single lineage (e.g., if this is the mrca of the metapopulation considered...
        for (i = 0; i < P.nbr_demes; ++i) {                            // ... or, more insterestingly, if this is the mrca of the full sample)
            if (nbr_lineages[i] > 0) break;                                // then find the deme where this lineage is
        }
        topnode = &deme[i].tree[(last_node[i] - 1)];        // 'topnode' is the address of this last lineage
    }
    free(last_node);
    free(nbr_lineages);
    free(p);
    return(topnode);
}

void add_mutations(parameters_struct P,
                   struct node_struct *node)

{
    int i,j;
    int nbr_mut,tmp;
    int time;
    
    for (i = 0; i < 2; ++i) {                         // This loop is over the total number of descendants of the current node
        if (node -> descendant[i] != NULL) {
            (node -> descendant[i]) -> allele = node -> allele; // Copy the ancestral allelic state into the current node
            time = node -> time - ((node -> descendant[i]) -> time); // Calculate the time elapsed since ancestor
      nbr_mut = (int) rpois(time * P.theta);        // Calculate the number of mutations along the branch
            for (j = 0; j < nbr_mut; ++j) {               // This loop is over all the mutations that occurred
                do {                                        // Draw a new allelic state
                    tmp = mutation(P,(node -> descendant[i]) -> allele); // bug fixed 17-01-2007 : was 'tmp = Mutation(node -> allele)'
                }                                           // which caused the number of mutations to be underestimated... (mutations different from the ancestral state, not from the
                while (tmp == (node -> descendant[i]) -> allele); // the new allele must be different from the current state
                (node -> descendant[i]) -> allele = tmp;    // Copy the mutation into the current node
            }
            add_mutations(P,node -> descendant[i]);       // Play it again (recursively, throughout the tree) starting from the descendants of the current node
        }
    }
}

int mutation(parameters_struct P,
             int allele)

{
    int new_allele = 0;
 
  switch (P.mutation_model) {
    
    case IAM :
      new_allele = ++k_iam;
      break;
      
    case KAM :
      if (P.nbr_allelic_states == 2) {
        new_allele = 1 - allele;
      }
      else {
        do {
          new_allele = (int) (P.nbr_allelic_states * unif_rand()); // draw an allelic state within a range of possible states ('model_characteristics' allelic states)
        }
        while (new_allele == allele);
      }
      break;
    
    case SMM :
      new_allele = allele + 2 * (int) (2 * unif_rand()) - 1; // calculate the sign (+1 or -1) of the mutation step, and add it to the current allelic state
      break;
  }
  return(new_allele);
}

void compute_locus_specific_statistics(parameters_struct P,
                                       statistics_struct *sum_stat)

{
  int i,j,k,n;
  int min,max,tmp;
  int *ss;
  double SSI,SSP,MSI,MSP,ss2,nc;
  double freq,homozygosity,meanlength,meansqrlength;
  int **count;
  
  min = 2147483647L;
  max = -2147483647L;
  for (i = 0; i < P.nbr_sampled_demes; ++i) {
    for (j = 0; j < deme[i].nbr_sampled_lineages; ++j) {
      if (deme[i].tree[j].allele < min) {
        min = deme[i].tree[j].allele;
      }
      if (deme[i].tree[j].allele > max) {
        max = deme[i].tree[j].allele;
      }
    }
  }
  for (i = 0; i < P.nbr_sampled_demes; ++i) {
    for (j = 0; j < deme[i].nbr_sampled_lineages; ++j) {
      deme[i].tree[j].allele -= min;
    }
  }
  sum_stat -> min = 0;
  sum_stat -> max = (max - min);
  count = (int **) malloc((P.nbr_sampled_demes + 1) * sizeof(int *)); // These are the allele frequencies
  for (i = 0; i < (P.nbr_sampled_demes + 1); ++i) {
    count[i] = (int *) malloc((sum_stat -> max + 1) * sizeof(int));
  }
  ss = (int *) malloc (P.nbr_sampled_demes * sizeof(int));
  for (i = 0; i < (P.nbr_sampled_demes + 1); ++i) {
    for (k = 0; k < (sum_stat -> max + 1); ++k) {
      count[i][k] = 0.0;
    }
  }
  n = ss2 = 0;                                      // n is a counter, should be equal to the total number of genes at the end of the loop
  for (i = 0; i < P.nbr_sampled_demes; ++i) {
    ss[i] = 0;                                      // ss[j] is the sample size (number of genes, twice the number of inividuals) of the jth deme
    for (j = 0; j < deme[i].nbr_sampled_lineages; ++j) { // deme[0][i][j].n is a number of genes (twice the number of individuals)
      ++count[i][deme[i].tree[j].allele];           // Calculate the allele frequencies
      ++count[P.nbr_sampled_demes][deme[i].tree[j].allele]; // Calculate the marginal allele frequencies
      ++ss[i];                                      // Calculate the sample sizes (number of sampled genes)
      ++n;
    }
    ss2 += ss[i] * ss[i];                           // ss2 is the sum, over demes, of squared sample sizes (number of genes, twice the number of inividuals)
  }
  nc = ((double) n - (double) ss2 / n) / ((double) P.nbr_sampled_demes - 1.0);
  SSI = SSP = 0.0;
  for (i = 0; i < P.nbr_sampled_demes; ++i) {
    for(k = 0; k < (sum_stat -> max + 1); ++k) {
      SSI += (double) count[i][k] - pow(count[i][k],2) / ss[i];
      SSP += ss[i] * pow(((double) count[i][k] / ss[i] - count[P.nbr_sampled_demes][k] / n),2);
    }
  }
  MSI = (double) SSI / (n - P.nbr_sampled_demes);
  MSP = (double) SSP / (P.nbr_sampled_demes - 1.0);
  sum_stat -> F_ST_num = (MSP - MSI);
  sum_stat -> F_ST_den = (MSP + (nc - 1) * MSI);
  sum_stat -> F_ST = (MSP - MSI)  / (MSP + (nc - 1) * MSI);
  sum_stat -> avrge_heterozygosity = 0.0;
  sum_stat -> avrge_nbr_alleles = 0.0;
  sum_stat -> avrge_shannon = 0.0;
  for (i = 0; i < P.nbr_sampled_demes; ++i) {
    homozygosity = 0.0;
    meanlength = meansqrlength = 0.0;
    tmp = 0;
    for (k = 0; k < (sum_stat -> max + 1); ++k) {
      if (count[i][k] > 0) {
        tmp += 1;
        freq = (double) count[i][k] / ss[i];
        sum_stat -> avrge_shannon += freq * log(freq);
      }
      homozygosity += pow((double) count[i][k] / ss[i],2);
    }
    sum_stat -> avrge_heterozygosity += (1.0 - homozygosity) * ss[i] / (ss[i] - 1.0);
    sum_stat -> avrge_nbr_alleles += (double) tmp;
    sum_stat -> avrge_shannon *= (-1.0);
  }
  sum_stat -> avrge_heterozygosity /= P.nbr_sampled_demes;
  sum_stat -> avrge_nbr_alleles /= P.nbr_sampled_demes;
  sum_stat -> avrge_var_allele_length = 0.0;
  for (i = 0; i < P.nbr_sampled_demes; ++i) {
    meanlength = meansqrlength = 0.0;
    for (j = 0; j < deme[i].nbr_sampled_lineages; ++j) {
      meanlength += (double) deme[i].tree[j].allele;
      meansqrlength += (double) pow(deme[i].tree[j].allele,2);
    }
    meanlength /= (double) ss[i];
    meansqrlength /= (double) ss[i];
    sum_stat -> avrge_var_allele_length += (meansqrlength - pow(meanlength,2)) * ss[i] / (ss[i] - 1.0);
  }
  sum_stat -> avrge_var_allele_length /= P.nbr_sampled_demes;
  sum_stat -> total_nbr_alleles = 0.0;
  sum_stat -> total_heterozygosity = 0.0;
  sum_stat -> total_shannon = 0.0;
  homozygosity = 0.0;
  for (k = 0; k < (sum_stat -> max + 1); ++k) {
    if (count[P.nbr_sampled_demes][k] > 0) {
      sum_stat -> total_nbr_alleles += 1.0;
      freq = (double) count[P.nbr_sampled_demes][k] / n;
      sum_stat -> total_shannon += freq * log(freq);
    }
    homozygosity += pow((double) count[P.nbr_sampled_demes][k] / n,2);
  }
  sum_stat -> total_heterozygosity = (1.0 - homozygosity) * n / (n - 1.0);
  sum_stat -> total_shannon *= (-1.0);
  meanlength = meansqrlength = 0.0;
  for (i = 0; i < P.nbr_sampled_demes; ++i) {
    for (j = 0; j < deme[i].nbr_sampled_lineages; ++j) {
      meanlength += (double) deme[i].tree[j].allele;
      meansqrlength += (double) pow(deme[i].tree[j].allele,2);
    }
  }
  meanlength /= (double) n;
  meansqrlength /= (double) n;
  sum_stat -> total_var_allele_length = (meansqrlength - pow(meanlength,2)) * n / (n - 1.0);
  for (i = 0; i < (P.nbr_sampled_demes + 1); ++i) {
    free(count[i]);
  }
  free(count);
  free(ss);
}

void compute_multilocus_statistics(parameters_struct P,
                                   statistics_struct *multi_locus,
                                   statistics_struct *locus_specific)

{
  int locus;
  
  multi_locus -> avrge_heterozygosity = 0.0;
  multi_locus -> total_heterozygosity = 0.0;
  multi_locus -> avrge_nbr_alleles = 0.0;
  multi_locus -> total_nbr_alleles = 0.0;
  multi_locus -> avrge_var_allele_length = 0.0;
  multi_locus -> total_var_allele_length = 0.0;
  multi_locus -> F_ST_num = 0.0;
  multi_locus -> F_ST_den = 0.0;
  multi_locus -> avrge_shannon = 0.0;
  multi_locus -> total_shannon = 0.0;
  for (locus = 0; locus < P.nbr_loci; ++locus) {
    multi_locus -> avrge_heterozygosity += locus_specific[locus].avrge_heterozygosity;
    multi_locus -> total_heterozygosity += locus_specific[locus].total_heterozygosity;
    multi_locus -> avrge_nbr_alleles += locus_specific[locus].avrge_nbr_alleles;
    multi_locus -> total_nbr_alleles += locus_specific[locus].total_nbr_alleles;
    multi_locus -> avrge_var_allele_length += locus_specific[locus].avrge_var_allele_length;
    multi_locus -> total_var_allele_length += locus_specific[locus].total_var_allele_length;
    multi_locus -> F_ST_num += locus_specific[locus].F_ST_num;
    multi_locus -> F_ST_den += locus_specific[locus].F_ST_den;
    multi_locus -> avrge_shannon += locus_specific[locus].avrge_shannon;
    multi_locus -> total_shannon += locus_specific[locus].total_shannon;
  }
  multi_locus -> avrge_heterozygosity /= P.nbr_loci;
  multi_locus -> total_heterozygosity /= P.nbr_loci;
  multi_locus -> avrge_nbr_alleles /= P.nbr_loci;
  multi_locus -> total_nbr_alleles /= P.nbr_loci;
  multi_locus -> avrge_var_allele_length /= P.nbr_loci;
  multi_locus -> total_var_allele_length /= P.nbr_loci;
  multi_locus -> avrge_shannon /= P.nbr_loci;
  multi_locus -> total_shannon /= P.nbr_loci;
  multi_locus -> F_ST = multi_locus -> F_ST_num / multi_locus -> F_ST_den;
}

void write_outputs(statistics_struct multi_locus)

{
  fprintf(outfile,"%10.6f\t",multi_locus.total_nbr_alleles);
  fprintf(outfile,"%10.6f\t",multi_locus.avrge_nbr_alleles);
  fprintf(outfile,"%10.6f\t",multi_locus.total_heterozygosity);
  fprintf(outfile,"%10.6f\t",multi_locus.avrge_heterozygosity);
  fprintf(outfile,"%10.6f\t",multi_locus.F_ST);
  fprintf(outfile,"%10.6f\t",multi_locus.total_var_allele_length);
  fprintf(outfile,"%10.6f\t",multi_locus.avrge_var_allele_length);
  fprintf(outfile,"%10.6f\t",multi_locus.total_shannon);
  fprintf(outfile,"%10.6f\n",multi_locus.avrge_shannon);
}

#ifdef DEBUG
void compute_locus_specific_identity_probabilities(parameters_struct P,
                                                   statistics_struct *sum_stat)

{
    int i,j,k,l;
    int npairs[2];
    double Q[2];
    
    npairs[WITHIN] = npairs[BETWEEN] = 0;
    sum_stat -> Q[WITHIN] = sum_stat -> Q[BETWEEN] = 0.0;
    for (i = 0; i < P.nbr_sampled_demes; ++i) {
        for (j = 0; j < (deme[i].nbr_sampled_lineages - 1); ++j) {
            for (k = (j + 1); k < deme[i].nbr_sampled_lineages; ++k) {
                if (deme[i].tree[j].allele == deme[i].tree[k].allele) sum_stat -> Q[WITHIN] += 1.0;
                npairs[WITHIN] += 1;
            }
        }
    }
    for (i = 0; i < (P.nbr_sampled_demes - 1); ++i) {
        for (j = (i + 1); j < P.nbr_sampled_demes; ++j) {
            for (k = 0; k < deme[i].nbr_sampled_lineages; ++k) {
                for (l = 0; l < deme[j].nbr_sampled_lineages; ++l) {
                    if (deme[i].tree[k].allele == deme[j].tree[l].allele) sum_stat -> Q[BETWEEN] += 1.0;
                    npairs[BETWEEN] += 1;
                }
            }
        }
    }
  sum_stat -> Q[WITHIN] /= (double) npairs[WITHIN];
  sum_stat -> Q[BETWEEN] /= (double) npairs[BETWEEN];
}

void compute_multilocus_identity_probabilities(parameters_struct P,
                                               statistics_struct *multi_locus,
                                               statistics_struct *locus_specific)

{
  int locus;

  multi_locus -> Q[WITHIN] = 0.0;
  multi_locus -> Q[BETWEEN] = 0.0;
  for (locus = 0; locus < P.nbr_loci; ++locus) {
    multi_locus -> Q[WITHIN] += locus_specific[locus].Q[WITHIN];
    multi_locus -> Q[BETWEEN] += locus_specific[locus].Q[BETWEEN];
  }
  multi_locus -> Q[WITHIN] /= P.nbr_loci;
  multi_locus -> Q[BETWEEN] /= P.nbr_loci;
}

void write_identity_probabilities(statistics_struct multi_locus)

{
  fprintf(outfile,"%10.6f\t",multi_locus.Q[WITHIN]);
  fprintf(outfile,"%10.6f\n",multi_locus.Q[BETWEEN]);
}
#endif

