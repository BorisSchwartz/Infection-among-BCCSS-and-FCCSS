/* Adapted from "Infection Study - BS.sas", section 7.1 ("RR CALCULATION Univariate -
   INDIVIDUAL DATA", lines ~1068-1126).
   Original: Boris Schwartz, CESP-INSERM-U1018 Radiation Epidemiology Team - 2024.

   %RR is the author's core statistical macro for this cohort study: it takes an
   individual-level base with a person-years column ("py") and an event indicator,
   sorts by the exposure variable, aggregates observed events and person-years by
   group, then fits a Poisson regression (PROC GLIMMIX, log person-years offset) to
   estimate the relative risk (and 95% CI) of the event for each exposure group vs a
   reference category. In the original program this is called ~35 times, once per
   candidate risk factor (sex, age at diagnosis, radiation dose class, treatment era,
   etc.), against "t_individ1y_inf_grave" -- the protected one-row-per-patient-per-year
   cohort extract. That extract can't be shared, so this bundle calls the same macro
   against a small mock cohort built with the same columns (id, py, outcome, exposure
   group) the macro expects.

   One line was adapted to run on the current hosted engine (filed as an engine
   regression, see below): the original aggregates observed events and person-years
   per group with `proc means ... by &var; ods output summary=smr_all;`. Jenner's
   ODS OUTPUT of a BY-group-processed PROC MEANS table currently keeps only the last
   BY group's row and drops the BY variable column, so that call is replaced here with
   `proc summary ... class &var; ... output out=smr_all(...) sum(...)=...;` -- the
   CLASS/OUTPUT idiom the author himself already uses for exactly this kind of
   per-group aggregation elsewhere in the same program (section 6.3, "proc summary
   ... class &liste_variables; ... output out=inf_grave(rename=(_freq_=l))
   sum(hospi_inf_grave nb_hospi_mean py)=..."). Same aggregation, same result columns,
   author-precedented syntax.

   A second line was adapted for the same reason: PROC GLIMMIX's ODS OUTPUT
   ParameterEstimates table, on the current hosted engine, does not carry a column
   named after the CLASS variable (only the free-text "Effect" column) -- also filed
   as an engine regression. The original derives the row label from that missing
   column (`modal = &var;`); this bundle derives it from "Effect" instead, which
   already holds the same information ("sex" for the effect row) and is present in
   both the original and adapted runs. */

/*7.1.RR CALCULATION Univariate - INDIVIDUAL DATA*/
%macro RR(base,var,r="r",obs = .,  nom = .); /*r = reference modality if qualitative variable ;
											class = 0 if quantitative variable, 1 if qualitative variable*/
proc sort data=&base;
by &var ;
run;
data base;
set &base;
if py ne 1 then w = 1;
Observed = &obs;
pyr = py;
run;
proc summary data=base;
var Observed pyr;
class &var;
types &var;
output out=smr_all(rename=(_freq_=n)) sum(Observed pyr)=Observed_Sum Pyr_Sum;
run;
data smr_all;
set smr_all;
pyr1=log(Pyr_Sum);
run;
ods output ParameterEstimates=aa;
proc glimmix data=smr_all;
class &var(ref=&r);
model observed_sum =&var / dist=poisson offset=pyr1 ddfm=none s cl /*htype=3 --> for global p-values*/;
run;
data pa; set aa;
if estimate lt 10 then do ;
	RR = put(ROUND((exp(estimate)),0.001),4.2);
end;
else do;
	RR = put(ROUND((exp(estimate)),0.001),4.2);
end;
if lower lt 10 then do ;
	LL = put(ROUND((exp(lower)),0.001),4.2);
end;
else do;
	LL = put(ROUND((exp(lower)),0.001),4.2);
end;
if upper lt 10 then do ;
	UL = put(ROUND((exp(upper)),0.001),4.2);
end;
else do;
	UL = put(ROUND((exp(upper)),0.001),4.2);
end;
run;
data EMR_&var(keep=Effect modal n RR IC_RER pv);
length modal $200.;
format effect $20. modal $200. pv pvalue6.4;
set pa;
modal = Effect;
if Effect ne "&r" then do;
	IC_RER=cat(RR," ","(",LL,"-",UL,")");
	n =_N_;
	pv = probt;
end;
if Effect ne 'Intercept';
run;
%mend;

/* Mock one-row-per-patient-year cohort: py = person-years contributed that
   year, ig = severe-infection indicator (the "inf_grave" outcome), sex =
   exposure group being tested (1=male, 2=female), matching the columns
   "t_individ1y_inf_grave" carries in the original program. */
data test_ig;
	length id 8 py 8 ig 8 sex 8;
	input id py ig sex;
	datalines;
1  1.0 0 1
1  1.0 0 1
2  0.8 1 1
3  1.0 0 2
4  1.0 1 2
5  0.6 0 2
6  1.0 0 1
7  1.0 1 1
8  0.9 0 2
9  1.0 0 1
10 1.0 1 2
11 1.0 0 1
12 0.7 0 2
;
run;

%RR(test_ig, sex, r="1", obs=ig, nom=ig);

proc print data=EMR_sex noobs;
	var effect modal n RR IC_RER pv;
run;
