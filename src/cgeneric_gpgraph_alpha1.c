#include "cgeneric_defs.h"
#include "stdio.h"

// This version uses 'padded' matrices with zeroes
double *inla_cgeneric_gpgraph_alpha1_model(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data) {

  double *ret = NULL;

  double lkappa, lsigma, kappa, sigma;

  double c1, c2, c_1, c_2, one_m_c2, l_e;

  int N, M, i, j, k;

  char *parameterization;
  
  // the size of the model
  assert(data->n_ints == 7);

  // the number of doubles
  assert(data->n_doubles == 9);

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

  assert(!strcasecmp(data->ints[4]->name, "index_graph"));
  inla_cgeneric_vec_tp *idx_ij = data->ints[4];
  int M_full = idx_ij->len;

  assert(!strcasecmp(data->ints[5]->name, "count_idx"));
  inla_cgeneric_vec_tp *count_idx = data->ints[5];
  assert(M == count_idx->len);

  assert(!strcasecmp(data->ints[6]->name, "stationary_endpoints"));
  inla_cgeneric_vec_tp *stationary_endpoints = data->ints[6];

  assert(!strcasecmp(data->doubles[0]->name, "EtV2"));
  inla_cgeneric_vec_tp *EtV2 = data->doubles[0];
  
  int nE = EtV2 -> len;

  assert(!strcasecmp(data->doubles[1]->name, "EtV3"));
  inla_cgeneric_vec_tp *EtV3 = data->doubles[1];

  assert(nE == EtV3 -> len);

  assert(!strcasecmp(data->doubles[2]->name, "El"));
  inla_cgeneric_vec_tp *El = data->doubles[2];

  // prior parameters
  assert(!strcasecmp(data->doubles[3]->name, "start_theta"));
  double start_theta = data->doubles[3]->doubles[0];

  assert(!strcasecmp(data->doubles[4]->name, "start_lsigma"));
  double start_lsigma = data->doubles[4]->doubles[0];

  assert(!strcasecmp(data->doubles[5]->name, "prior_theta_meanlog"));
  double prior_theta_meanlog = data->doubles[5]->doubles[0];

  assert(!strcasecmp(data->doubles[6]->name, "prior_theta_sdlog"));
  double prior_theta_sdlog = data->doubles[6]->doubles[0];

  assert(!strcasecmp(data->doubles[7]->name, "prior_sigma_meanlog"));
  double prior_sigma_meanlog = data->doubles[7]->doubles[0];

  assert(!strcasecmp(data->doubles[8]->name, "prior_sigma_sdlog"));
  double prior_sigma_sdlog = data->doubles[8]->doubles[0];

  assert(!strcasecmp(data->chars[2]->name, "parameterization"));
  parameterization = &data->chars[2]->chars[0];

  if (theta) {
    // interpretable parameters 

    if(!strcasecmp(parameterization, "matern")){
      lkappa = log(2.0) - theta[1];
    } else {
      lkappa = theta[1];
    }
    lsigma = theta[0];
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
      ret[0] = N;       /* dimension */
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

      double *raw_entries;
      raw_entries = Calloc(M_full, double);

      // int count=0;
      //     for(i = 0; i < nE; i++){
      //       l_e = El->doubles[i];
      //       c1 = exp(-kappa*l_e);
      //       c2 = SQR(c1);
      //       one_m_c2 = 1-c2;
      //       c_1 = 0.5 + c2/one_m_c2;
      //       c_2 = -c1/one_m_c2;

      //       if(EtV2->doubles[i] != EtV3->doubles[i]){
            
      //       ret[k + idx_ij->ints[count]] = c_1;
      //       ret[k + idx_ij->ints[count + 1]] = c_1;
      //       ret[k + idx_ij->ints[count + 2]] = c_2;
      //       count += 3;
      //     }else{
      //       ret[k + idx_ij->ints[count]] = tanh(0.5 * kappa * l_e);
      //       count++;
      //     }
      //   }

      //   if(stationary_endpoints->ints[0]>=0){
      //     int stat_endpt_len = stationary_endpoints->len;
      //     for(i = 0; i < stat_endpt_len; i++){
      //       ret[k+idx_ij->ints[count+i]] = 0.5;
      //     }
      //     count += stat_endpt_len;
      //   }

      //   assert(count == M);

      //   double fact = 2*kappa / (pow(sigma,2));

      //   int one=1;
      //   dscal_(&M, &fact, &ret[k], &one);
    
       int count=0;
          for(i = 0; i < nE; i++){
            l_e = El->doubles[i];
            c1 = exp(-kappa*l_e);
            c2 = SQR(c1);
            one_m_c2 = 1-c2;
            c_1 = 0.5 + c2/one_m_c2;
            c_2 = -c1/one_m_c2;

            if(EtV2->doubles[i] != EtV3->doubles[i]){
            
            raw_entries[idx_ij->ints[count]] = c_1;
            raw_entries[idx_ij->ints[count + 1]] = c_1;
            raw_entries[idx_ij->ints[count + 2]] = c_2;
            count += 3;
          }else{
            raw_entries[idx_ij->ints[count]] = tanh(0.5 * kappa * l_e);
            count++;
          }
        }

        if(stationary_endpoints->ints[0]>=0){
          int stat_endpt_len = stationary_endpoints->len;
          for(i = 0; i < stat_endpt_len; i++){
            raw_entries[idx_ij->ints[count+i]] = 0.5;
          }
          count += stat_endpt_len;
        }

        assert(count == M_full);

        double fact = 2*kappa / (pow(sigma,2));

        int one=1;
        dscal_(&M_full, &fact, &raw_entries[0], &one);

        count = 0;
        for(i = 0; i < M; i++){
          for(j = 0; j < count_idx->ints[i]; j++){
            // ret[k + i] += raw_entries[count];
            ret[k + i] += raw_entries[count];
            count++;
          }
        }
        assert(M_full == count);

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
      ret[1] = start_lsigma;
      ret[2] = start_theta;
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

      ret[0] += -0.5 * SQR(theta[1] - prior_theta_meanlog)/(SQR(prior_theta_sdlog)) - 
      log(prior_theta_sdlog) - 0.5 * log(2.0 * M_PI); 

      ret[0] += -0.5 * SQR(lsigma - prior_sigma_meanlog)/(SQR(prior_sigma_sdlog)) - 
      log(prior_sigma_sdlog) - 0.5 * log(2.0 * M_PI);
	    break;
    }
    
  case INLA_CGENERIC_QUIT:
  default:
    break;
  }
  
  return (ret);
}