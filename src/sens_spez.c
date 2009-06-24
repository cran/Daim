#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>


SEXP one_roc(SEXP prob_oob, SEXP label, SEXP thres, SEXP labn)
{
	int i, j, N_pos, N_neg;
	double FP, TP;
	int N_thres = LENGTH(thres);
	double *th = REAL(thres);
	double *prob = REAL(prob_oob);
	double *class = REAL(label);
	
	SEXP FP_Rate, TP_Rate, result, names_result, row_names;
	
	PROTECT(FP_Rate = allocVector(REALSXP,N_thres));
	PROTECT(TP_Rate = allocVector(REALSXP,N_thres));
	
	for(i=0; i < N_thres; i++){
		TP=0.0, FP=0.0, N_pos=0;
		for(j=0; j < LENGTH(prob_oob); j++){
			TP += (th[i] <= prob[j]) && (class[j] != REAL(labn)[0]);
			FP += (th[i] <= prob[j]) && (class[j] == REAL(labn)[0]);
			N_pos += class[j] != REAL(labn)[0];
		}
		N_neg = LENGTH(prob_oob) - N_pos;
		REAL(FP_Rate)[i] = FP / N_neg;
		REAL(TP_Rate)[i] = TP / N_pos;
	}
	PROTECT(result = allocVector(VECSXP,2));
	PROTECT(names_result = allocVector(STRSXP, 2));
	SET_STRING_ELT(names_result, 0, mkChar("FPR"));
	SET_STRING_ELT(names_result, 1, mkChar("TPR"));
	setAttrib(result, R_NamesSymbol, names_result);
	UNPROTECT(1);
	
	PROTECT(row_names = allocVector(INTSXP, 2));
	INTEGER(row_names)[0] = NA_INTEGER;
	INTEGER(row_names)[1] = N_thres;
	setAttrib(result, R_RowNamesSymbol, row_names);
	UNPROTECT(1);
	
	SET_VECTOR_ELT(result, 0, FP_Rate);
	SET_VECTOR_ELT(result, 1, TP_Rate);
	UNPROTECT(3);
	setAttrib(result, R_ClassSymbol, mkString("data.frame"));
	return(result);
}

SEXP roc_value( SEXP prob, SEXP labels, SEXP thres, SEXP labn)
{
	int i;
	int N_obj = LENGTH(prob);
	
	SEXP result = PROTECT(allocVector(VECSXP,N_obj));
	for(i=0; i < N_obj; i++){
		SET_VECTOR_ELT(result, i, one_roc(VECTOR_ELT(prob,i), VECTOR_ELT(labels,i), thres, labn));
	}
	UNPROTECT(1);
	return(result);
}

double error_loob( SEXP loob_i, SEXP classt, double rsp_i)
{
	int i;
	double ans=0.0;
	for(i=0; i < LENGTH(loob_i); i++){
		ans += (REAL(loob_i)[i] >= REAL(classt)[0]) != rsp_i;
	}
	ans = ans / i;
	return(ans);
}


double error_loob2( SEXP loob_i, double classt, double rsp_i)
{
	int i;
	double ans=0.0;
	for(i=0; i < LENGTH(loob_i); i++){
		ans += (REAL(loob_i)[i] >= classt) != rsp_i;
	}
	ans = ans / i;
	return(ans);
}


SEXP sens_spez_obs (SEXP loob, SEXP app, SEXP thres, SEXP lab_neg, SEXP classt, SEXP rsp)
{
	int i,j,t,sum_tlab=0;
	double sum_loob=0.0, err_app=0.0, err_loob=0.0;
	SEXP temp_loob;
	
	int N_lab = LENGTH(rsp);
	int N_thres = LENGTH(thres);
	double *thr = REAL(thres);
	double *prob_app = REAL(app);
	double *RSP = REAL(rsp);
	double neg_lab = REAL(lab_neg)[0]; 
	
	double sensspezloob[N_lab][N_thres];
	int *templab;
	templab = Calloc(N_lab,int);
	
	if(N_lab > LENGTH(loob))
		 error("\n Increase nboot and try again! \n");
	for(i=0; i < N_lab; i++){
		err_app += (prob_app[i] >= REAL(classt)[0]) != (RSP[i] - 1.0);
		templab[i] = (RSP[i] - 1.0) == neg_lab;
		sum_tlab += templab[i];
		temp_loob = VECTOR_ELT(loob,i);
		err_loob += error_loob(temp_loob,classt,(RSP[i]-1.0));
		for(j=0; j < N_thres; j++){
			sensspezloob[i][j] = 0.0;
			sum_loob = 0.0;
			if(templab[i]){
				for(t=0; t < LENGTH(temp_loob); t++){
					sum_loob += (REAL(temp_loob)[t] >= thr[j]) == neg_lab; 
				}
			}
			else{
				for(t=0; t < LENGTH(temp_loob); t++){
					sum_loob += (REAL(temp_loob)[t] >= thr[j]) != neg_lab; 
				}
			}
			sensspezloob[i][j] = sum_loob / LENGTH(temp_loob);
		}
	}

	double *spez_loob, *sens_loob;
	double *spez_app, *sens_app, *my_q;
	spez_loob = Calloc(N_thres,double);
	sens_loob = Calloc(N_thres,double);
	spez_app = Calloc(N_thres,double);
	sens_app = Calloc(N_thres,double);
	my_q = Calloc(N_thres,double);
	
	int neg=0, pos=0;
	
	for(j=0; j < N_thres; j++){
		spez_loob[j] = 0.0;
		sens_loob[j] = 0.0;
		spez_app[j] = 0.0;
		sens_app[j] = 0.0;
		my_q[j] = 0.0;
		for(i=0; i < N_lab; i++){
			my_q[j] += prob_app[i] >= thr[j];
			if(templab[i]){
				spez_loob[j] += sensspezloob[i][j];
				spez_app[j] += prob_app[i] < thr[j];
				neg++;
			}
			else{
				sens_loob[j] += sensspezloob[i][j];
				sens_app[j] += prob_app[i] >= thr[j];
				pos++;				
			}
		}
	}
	
	SEXP result, names_result, app_err, loob_err;
	SEXP sensloob, spezloob, sensapp, spezapp, myq, myp;
	
	PROTECT(result = allocVector(VECSXP,8));
	PROTECT(names_result = allocVector(STRSXP, 8));
	SET_STRING_ELT(names_result, 0, mkChar("sensloob"));
	SET_STRING_ELT(names_result, 1, mkChar("spezloob"));
	SET_STRING_ELT(names_result, 2, mkChar("sensapp"));
	SET_STRING_ELT(names_result, 3, mkChar("spezapp"));
	SET_STRING_ELT(names_result, 4, mkChar("myq"));
	SET_STRING_ELT(names_result, 5, mkChar("myp"));
	SET_STRING_ELT(names_result, 6, mkChar("errapp"));
	SET_STRING_ELT(names_result, 7, mkChar("errloob"));
	setAttrib(result, R_NamesSymbol, names_result);
	UNPROTECT(1);
	PROTECT(sensloob = allocVector(REALSXP,N_thres));
	PROTECT(spezloob = allocVector(REALSXP,N_thres));
	PROTECT(sensapp = allocVector(REALSXP,N_thres));
	PROTECT(spezapp = allocVector(REALSXP,N_thres));
	PROTECT(myq = allocVector(REALSXP,N_thres));
	PROTECT(myp = allocVector(REALSXP,1));
	PROTECT(app_err = allocVector(REALSXP,1));
	PROTECT(loob_err = allocVector(REALSXP,1));
	for(j=0; j < N_thres; j++){
		REAL(sensloob)[j] = sens_loob[j]/pos*N_thres;
		REAL(spezloob)[j] = spez_loob[j]/neg*N_thres;
		REAL(sensapp)[j] = sens_app[j] / (double) (N_lab - sum_tlab);
		REAL(spezapp)[j] = spez_app[j] / (double) sum_tlab;
		REAL(myq)[j] = my_q[j] / (double) N_lab;
	}
	UNPROTECT(5);
	REAL(myp)[0] = (double) (N_lab - sum_tlab) / (double) N_lab;
	REAL(app_err)[0] = err_app / N_lab;
	REAL(loob_err)[0] = err_loob / N_lab;
	SET_VECTOR_ELT(result, 0, sensloob);
	SET_VECTOR_ELT(result, 1, spezloob);
	SET_VECTOR_ELT(result, 2, sensapp);
	SET_VECTOR_ELT(result, 3, spezapp);
	SET_VECTOR_ELT(result, 4, myq);
	SET_VECTOR_ELT(result, 5, myp);
	SET_VECTOR_ELT(result, 6, app_err);
	SET_VECTOR_ELT(result, 7, loob_err);
	UNPROTECT(4);
	Free(spez_loob);Free(sens_loob);Free(spez_app);Free(sens_app);
	Free(my_q); Free(templab);
	return(result);
}

SEXP sens_spez_none (SEXP loob, SEXP app, SEXP thres, SEXP lab_neg, SEXP classt, SEXP rsp)
{
	int i,j,sum_tlab=0;
	double err_app=0.0, err_loob=0.0;
	SEXP temp_loob;
	
	int N_lab = LENGTH(rsp);
	int N_thres = LENGTH(thres);
	double *thr = REAL(thres);
	double *prob_app = REAL(app);
	double *RSP = REAL(rsp);
	double neg_lab = REAL(lab_neg)[0]; 
	
	double *mean_loob;
	int *templab;
	mean_loob = Calloc(N_lab,double);
	templab = Calloc(N_lab,int);
	
	if(N_lab > LENGTH(loob))
		 error("\n Increase nboot and try again! \n");
	for(i=0; i < N_lab; i++){
		err_app += (prob_app[i] >= REAL(classt)[0]) != (RSP[i] - 1.0);
		templab[i] = (RSP[i] - 1.0) == neg_lab;
		sum_tlab += templab[i];
		temp_loob = VECTOR_ELT(loob,i);
		err_loob += error_loob(temp_loob,classt,(RSP[i]-1.0));
		mean_loob[i] = 0.0;
		for(j=0; j < LENGTH(temp_loob); j++){
			mean_loob[i] += REAL(temp_loob)[j];
		}
		mean_loob[i] = mean_loob[i] / LENGTH(temp_loob);
	}

	double *spez_loob, *sens_loob;
	double *spez_app, *sens_app, *my_q;
	spez_loob = Calloc(N_thres,double);
	sens_loob = Calloc(N_thres,double);
	spez_app = Calloc(N_thres,double);
	sens_app = Calloc(N_thres,double);
	my_q = Calloc(N_thres,double);
	int neg=0, pos=0;

	for(j=0; j < N_thres; j++){
		spez_loob[j] = 0.0;
		sens_loob[j] = 0.0;
		spez_app[j] = 0.0;
		sens_app[j] = 0.0;
		my_q[j] = 0.0;
		for(i=0; i < N_lab; i++){
			my_q[j] += prob_app[i] >= thr[j];
			if(templab[i]){
				spez_loob[j] += mean_loob[i] < thr[j];
				spez_app[j] += prob_app[i] < thr[j];
				neg++;
			}
			else{
				sens_loob[j] += mean_loob[i] >= thr[j];
				sens_app[j] += prob_app[i] >= thr[j];
				pos++;				
			}
		}
	}

	SEXP result, names_result, app_err, loob_err;
	SEXP sensloob, spezloob, sensapp, spezapp, myq, myp;
	
	PROTECT(result = allocVector(VECSXP,8));
	PROTECT(names_result = allocVector(STRSXP, 8));
	SET_STRING_ELT(names_result, 0, mkChar("sensloob"));
	SET_STRING_ELT(names_result, 1, mkChar("spezloob"));
	SET_STRING_ELT(names_result, 2, mkChar("sensapp"));
	SET_STRING_ELT(names_result, 3, mkChar("spezapp"));
	SET_STRING_ELT(names_result, 4, mkChar("myq"));
	SET_STRING_ELT(names_result, 5, mkChar("myp"));
	SET_STRING_ELT(names_result, 6, mkChar("errapp"));
	SET_STRING_ELT(names_result, 7, mkChar("errloob"));
	setAttrib(result, R_NamesSymbol, names_result);
	UNPROTECT(1);
	PROTECT(sensloob = allocVector(REALSXP,N_thres));
	PROTECT(spezloob = allocVector(REALSXP,N_thres));
	PROTECT(sensapp = allocVector(REALSXP,N_thres));
	PROTECT(spezapp = allocVector(REALSXP,N_thres));
	PROTECT(myq = allocVector(REALSXP,N_thres));
	PROTECT(myp = allocVector(REALSXP,1));
	PROTECT(app_err = allocVector(REALSXP,1));
	PROTECT(loob_err = allocVector(REALSXP,1));
	for(j=0; j < N_thres; j++){
		REAL(sensloob)[j] = sens_loob[j]/pos*N_thres;
		REAL(spezloob)[j] = spez_loob[j]/neg*N_thres;
		REAL(sensapp)[j] = sens_app[j]/(N_lab - sum_tlab);
		REAL(spezapp)[j] = spez_app[j]/sum_tlab;
		REAL(myq)[j] = my_q[j] / (double) N_lab;
	}
	UNPROTECT(5);
	REAL(myp)[0] = (double) (N_lab - sum_tlab) / (double) N_lab;
	REAL(app_err)[0] = err_app / (double) N_lab;
	REAL(loob_err)[0] = err_loob / (double) N_lab;
	SET_VECTOR_ELT(result, 0, sensloob);
	SET_VECTOR_ELT(result, 1, spezloob);
	SET_VECTOR_ELT(result, 2, sensapp);
	SET_VECTOR_ELT(result, 3, spezapp);
	SET_VECTOR_ELT(result, 4, myq);
	SET_VECTOR_ELT(result, 5, myp);
	SET_VECTOR_ELT(result, 6, app_err);
	SET_VECTOR_ELT(result, 7, loob_err);
	Free(spez_loob);Free(sens_loob);Free(spez_app);Free(sens_app);
	Free(my_q); Free(mean_loob); Free(templab);
	
	UNPROTECT(4);
	return(result);
}



SEXP sens_spez_obs_cut (SEXP loob, SEXP app, SEXP thres, SEXP lab_neg, SEXP classt, SEXP rsp)
{
	int i,j,t,sum_tlab=0;
	double sum_loob=0.0, err_app=0.0, err_loob=0.0;
	SEXP temp_loob;
	
	int N_lab = LENGTH(rsp);
	int N_thres = LENGTH(thres);
	double *thr = REAL(thres);
	double *prob_app = REAL(app);
	double *RSP = REAL(rsp);
	double neg_lab = REAL(lab_neg)[0]; 
	
	double sensspezloob[N_lab][N_thres];
	int *templab;
	double *err_app2, *err_loob2;
	templab = Calloc(N_lab,int);
	err_app2 = Calloc(N_thres,double);
	err_loob2 = Calloc(N_thres,double);
	
	if(N_lab > LENGTH(loob))
		 error("\n Increase nboot and try again! \n");
	for(j=0; j < N_thres; j++){
		err_app2[j] = 0.0;
		err_loob2[j] = 0.0;
		for(i=0; i < N_lab; i++){
			temp_loob = VECTOR_ELT(loob,i);
			err_app2[j] += (prob_app[i] >= thr[j]) != (RSP[i] - 1.0);
			err_loob2[j] += error_loob2(temp_loob,thr[j],(RSP[i]-1.0));
		}
	}
	for(i=0; i < N_lab; i++){
//		err_app += (prob_app[i] >= REAL(classt)[0]) != (RSP[i] - 1.0);
		templab[i] = (RSP[i] - 1.0) == neg_lab;
		sum_tlab += templab[i];
		temp_loob = VECTOR_ELT(loob,i);
//		err_loob += error_loob(temp_loob,classt,(RSP[i]-1.0));

		for(j=0; j < N_thres; j++){
			sensspezloob[i][j] = 0.0;
			sum_loob = 0.0;
			if(templab[i]){
				for(t=0; t < LENGTH(temp_loob); t++){
					sum_loob += (REAL(temp_loob)[t] >= thr[j]) == neg_lab; 
				}
			}
			else{
				for(t=0; t < LENGTH(temp_loob); t++){
					sum_loob += (REAL(temp_loob)[t] >= thr[j]) != neg_lab; 
				}
			}
			sensspezloob[i][j] = sum_loob / LENGTH(temp_loob);
		}
	}

	double *spez_loob, *sens_loob;
	double *spez_app, *sens_app, *my_q;
	spez_loob = Calloc(N_thres,double);
	sens_loob = Calloc(N_thres,double);
	spez_app = Calloc(N_thres,double);
	sens_app = Calloc(N_thres,double);
	my_q = Calloc(N_thres,double);
	
	int neg=0, pos=0;
	
	for(j=0; j < N_thres; j++){
		spez_loob[j] = 0.0;
		sens_loob[j] = 0.0;
		spez_app[j] = 0.0;
		sens_app[j] = 0.0;
		my_q[j] = 0.0;
		for(i=0; i < N_lab; i++){
			my_q[j] += prob_app[i] >= thr[j];
			if(templab[i]){
				spez_loob[j] += sensspezloob[i][j];
				spez_app[j] += prob_app[i] < thr[j];
				neg++;
			}
			else{
				sens_loob[j] += sensspezloob[i][j];
				sens_app[j] += prob_app[i] >= thr[j];
				pos++;				
			}
		}
	}
	
	SEXP result, names_result, app_err, loob_err, app_err2, loob_err2;
	SEXP sensloob, spezloob, sensapp, spezapp, myq, myp;
	
	PROTECT(result = allocVector(VECSXP,10));
	PROTECT(names_result = allocVector(STRSXP, 10));
	SET_STRING_ELT(names_result, 0, mkChar("sensloob"));
	SET_STRING_ELT(names_result, 1, mkChar("spezloob"));
	SET_STRING_ELT(names_result, 2, mkChar("sensapp"));
	SET_STRING_ELT(names_result, 3, mkChar("spezapp"));
	SET_STRING_ELT(names_result, 4, mkChar("myq"));
	SET_STRING_ELT(names_result, 5, mkChar("myp"));
	SET_STRING_ELT(names_result, 6, mkChar("errapp"));
	SET_STRING_ELT(names_result, 7, mkChar("errloob"));
	SET_STRING_ELT(names_result, 8, mkChar("errapp2"));
	SET_STRING_ELT(names_result, 9, mkChar("errloob2"));
	
	setAttrib(result, R_NamesSymbol, names_result);
	UNPROTECT(1);
	PROTECT(sensloob = allocVector(REALSXP,N_thres));
	PROTECT(spezloob = allocVector(REALSXP,N_thres));
	PROTECT(sensapp = allocVector(REALSXP,N_thres));
	PROTECT(spezapp = allocVector(REALSXP,N_thres));
	PROTECT(myq = allocVector(REALSXP,N_thres));
	PROTECT(myp = allocVector(REALSXP,1));
	PROTECT(app_err = allocVector(REALSXP,1));
	PROTECT(loob_err = allocVector(REALSXP,1));
	PROTECT(app_err2 = allocVector(REALSXP,N_thres));
	PROTECT(loob_err2 = allocVector(REALSXP,N_thres));
	for(j=0; j < N_thres; j++){
		REAL(sensloob)[j] = sens_loob[j]/pos*N_thres;
		REAL(spezloob)[j] = spez_loob[j]/neg*N_thres;
		REAL(sensapp)[j] = sens_app[j] / (double) (N_lab - sum_tlab);
		REAL(spezapp)[j] = spez_app[j] / (double) sum_tlab;
		REAL(myq)[j] = my_q[j] / (double) N_lab;
		REAL(app_err2)[j] = err_app2[j] / N_lab;
		REAL(loob_err2)[j] = err_loob2[j] / N_lab;
	}
	UNPROTECT(7);
	REAL(myp)[0] = (double) (N_lab - sum_tlab) / (double) N_lab;
	REAL(app_err)[0] = err_app / N_lab;
	REAL(loob_err)[0] = err_loob / N_lab;
	SET_VECTOR_ELT(result, 0, sensloob);
	SET_VECTOR_ELT(result, 1, spezloob);
	SET_VECTOR_ELT(result, 2, sensapp);
	SET_VECTOR_ELT(result, 3, spezapp);
	SET_VECTOR_ELT(result, 4, myq);
	SET_VECTOR_ELT(result, 5, myp);
	SET_VECTOR_ELT(result, 6, app_err);
	SET_VECTOR_ELT(result, 7, loob_err);
	SET_VECTOR_ELT(result, 8, app_err2);
	SET_VECTOR_ELT(result, 9, loob_err2);
	UNPROTECT(4);
	Free(spez_loob);Free(sens_loob);Free(spez_app);Free(sens_app);
	Free(my_q); Free(templab);Free(err_app2);Free(err_loob2);
	return(result);
}


