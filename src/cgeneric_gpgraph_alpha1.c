#include "cgeneric_defs.h"
#include "stdio.h"

// This version uses 'padded' matrices with zeroes
double *inla_cgeneric_gpgraph_alpha1_model(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data) {

  double *ret = NULL;

  double lkappa, lsigma, kappa, sigma;

  double c1, c2, c_1, c_2, one_m_c2, l_e;

  int N, M, i, j, k;
  
  // the size of the model
  assert(data->n_ints == 5);

  // the number of doubles
  assert(data->n_doubles == 8);

  assert(!strcasecmp(data->ints[0]->name, "n"));       // this will always be the case
  N = data->ints[0]->ints[0];			       // this will always be the case
  assert(N > 0);

  assert(!strcasecmp(data->ints[1]->name, "debug"));    // this will always be the case
  int debug = data->ints[1]->ints[0];	        // this will always be the case

  if(debug == 1){
    debug = 1;
  }

  assert(!strcasecmp(data->ints[2]->name, "prec_graph_i"));
  inla_cgeneric_vec_tp *graph_i = data->ints[2];
  M = graph_i->len;

  assert(!strcasecmp(data->ints[3]->name, "prec_graph_j"));
  inla_cgeneric_vec_tp *graph_j = data->ints[3];
  assert(M == graph_j->len);

  assert(!strcasecmp(data->doubles[0]->name, "EtV2"));
  inla_cgeneric_vec_tp *EtV2 = data->doubles[0];
  
  int nE = EtV2 -> len;

  assert(!strcasecmp(data->doubles[1]->name, "EtV3"));
  inla_cgeneric_vec_tp *EtV3 = data->doubles[1];

  assert(!strcasecmp(data->doubles[2]->name, "El"));
  inla_cgeneric_vec_tp *El = data->doubles[2];

  assert(!strcasecmp(data->ints[4]->name, "prec_graph_i"));
  inla_cgeneric_vec_tp *stationary_endpoints = data->ints[4];

  // prior parameters
  assert(!strcasecmp(data->doubles[3]->name, "start_kappa"));
  double start_kappa = data->doubles[3]->doubles[0];

  assert(!strcasecmp(data->doubles[4]->name, "start_sigma"));
  double start_sigma = data->doubles[4]->doubles[0];

  double start_lkappa = log(start_kappa);
  double start_lsigma = log(start_sigma);


  if (theta) {
    // interpretable parameters 
    lkappa = theta[0];
    lsigma = theta[1];
    kappa = exp(lkappa);
    sigma = exp(lsigma);
  }
  else {   
    lsigma = lkappa = sigma = kappa = NAN;
  }
  
  switch (cmd) {
  case INLA_CGENERIC_VOID:
    {
      assert(!(cmd == INLA_CGENERIC_VOID));
      break;
    }
    
  case INLA_CGENERIC_GRAPH:
    {
      k=2;
      ret = Calloc(k + 2 * M, double);
      ret[0] = N;        	       /* dimension */
      ret[1] = M;		   /* number of (i <= j) */
      for (i = 0; i < M; i++) {
	      ret[k++] = graph_i->ints[i];
      }
      for (i = 0; i < M; i++) {
	      ret[k++] = graph_j->ints[i];
      }
      break;
    }
    
  case INLA_CGENERIC_Q:
    {
      k = 2;
      ret = Calloc(k + M, double);
      ret[0] = -1;		/* REQUIRED */
      ret[1] = M;		/* REQUIRED */

     

      break;
    }
    
  case INLA_CGENERIC_MU:
    {
      ret = Calloc(1, double);
      ret[0] = 0.0;
      break;
    }
    
  case INLA_CGENERIC_INITIAL:
    {
      // return c(P, initials)
      // where P is the number of hyperparameters      
      ret = Calloc(3, double);
      ret[0] = 2;
      ret[1] = start_lkappa;
      ret[2] = start_lsigma;
      break;
    }
    
  case INLA_CGENERIC_LOG_NORM_CONST:
    {
      break;
    }
    
  case INLA_CGENERIC_LOG_PRIOR:
    {
      ret = Calloc(1, double);

      ret[0] = 0.0;

      ret[0] += -0.5 * SQR(theta[0] - start_lkappa)/(10.0) - 
      0.5*log(10) - 0.5 * log(2.0 * M_PI);

      ret[0] += -0.5 * SQR(theta[1] - start_lsigma)/(10) - 
      0.5*log(10) - 0.5 * log(2.0 * M_PI);

	    break;
    }
    
  case INLA_CGENERIC_QUIT:
  default:
    break;
  }
  
  return (ret);
}