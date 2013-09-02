/*
 *  rocDaim.c
 *  Daim
 *
 *  Created by Sergej Potapov on 07.06.10.
 *
 */


#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>


SEXP rocDaim(SEXP prob_oob, SEXP label, SEXP thres, SEXP labn)
{
	int i, j, N_pos, N_neg;
	double FP, TP;
	int N_thres = LENGTH(thres);
	double *th = REAL(thres);
	double *prob = REAL(prob_oob);
	double *class = REAL(label);
	
	SEXP FP_Rate, TP_Rate, cutoff, result, names_result, row_names;
	
	PROTECT(FP_Rate = allocVector(REALSXP,N_thres+1));
	PROTECT(TP_Rate = allocVector(REALSXP,N_thres+1));
	PROTECT(cutoff = allocVector(REALSXP,N_thres+1));
	
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
		REAL(cutoff)[i] = th[i];
	}
	REAL(FP_Rate)[i]=0.0, REAL(TP_Rate)[i]=0.0, REAL(cutoff)[i]=NA_REAL;
	
	PROTECT(result = allocVector(VECSXP,3));
	PROTECT(names_result = allocVector(STRSXP, 3));
	SET_STRING_ELT(names_result, 0, mkChar("FPR"));
	SET_STRING_ELT(names_result, 1, mkChar("TPR"));
	SET_STRING_ELT(names_result, 2, mkChar("cutoff"));
	setAttrib(result, R_NamesSymbol, names_result);
	UNPROTECT(2);
	
	PROTECT(row_names = allocVector(INTSXP, 2));
	INTEGER(row_names)[0] = NA_INTEGER;
	INTEGER(row_names)[1] = N_thres+1;
	setAttrib(result, R_RowNamesSymbol, row_names);
	UNPROTECT(1);
	
	SET_VECTOR_ELT(result, 0, FP_Rate);
	SET_VECTOR_ELT(result, 1, TP_Rate);
	SET_VECTOR_ELT(result, 2, cutoff);
	UNPROTECT(3);
	setAttrib(result, R_ClassSymbol, mkString("data.frame"));
	return(result);
}

