/* Adapted from "Infection Study - BS.sas", section 3 ("Frequently used macros", lines ~133-164)
   and the TYPEG format defined in section 1 (proc format, lines ~78-90).
   Original: Boris Schwartz, CESP-INSERM-U1018 Radiation Epidemiology Team - 2024.

   %supp, %svg, %tri and %QL are small utility macros the author defines once and
   then calls throughout the ~1800-line program to delete work tables, save a table
   into a named library, sort a table, and print descriptive statistics for a
   qualitative variable. They don't depend on the protected cohort data, so they are
   reproduced here unmodified and exercised end-to-end (sort -> describe -> save ->
   delete) against a small mock table built from the tumour-group codes (TYPEG
   format) used throughout the original program. */

proc format;
value TYPEG     1='Nephroblastome  '
                2='Neuroblastome'
                3='Lymphome'
                30='Hodgkin'
                31='LMNH'
                4='Tissus mous'
                5='Os '
                6='Cerveau'
                7='Gonades'
                8='Thyro�de'
                9='Retinoblastome'
                0='Autres'
				99='Overall';
run;

/*3.frequently used macros*/
/*deleting table*/
%macro supp(tab);
proc delete data=&tab;
run;
%mend;
/*saving table*/
%macro svg(tab,repe,nom);
data &repe..&nom;
set &tab;
run;
%mend;
/*sort table*/
%macro tri(tab,var);
proc sort data=&tab;
by &var;
run;
%mend;
/*descriptive statistics for QUALITATIVE variable*/
%macro QL(tab,var,nom,g);
proc freq data=&tab;
table &var;
run;
%mend;

/* Mock cohort slice: one row per patient, typeg = tumour-group code (see
   TYPEG format above), sex = 1/2. Stands in for the protected patient-level
   extract the original macros are called against ("t_individ1y_inf_grave"). */
data cohort;
	length id 8 typeg 8 sex 8;
	input id typeg sex;
	format typeg typeg.;
	datalines;
1 1 1
2 2 2
3 2 1
4 4 1
5 5 2
6 1 2
7 8 1
8 2 1
9 6 2
10 1 1
;
run;

/* Exercise the utility macros the way the original program does: sort, then
   describe the qualitative variable, then save the sorted table to a
   library, then delete the intermediate work copy. */
%tri(cohort, typeg);
%QL(cohort, typeg, "Tumour group", 0);

libname outlib "./output";
%svg(cohort, outlib, cohort_sorted);

proc print data=outlib.cohort_sorted noobs;
	var id typeg sex;
run;

%supp(cohort);
