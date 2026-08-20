/* Adapted from "Infection Study - BS.sas", section 6.3 ("Aggregate data") and 6.3
   ("Data management of aggregate data" / SMR-AER calculation), lines ~618-797.
   Original: Boris Schwartz, CESP-INSERM-U1018 Radiation Epidemiology Team - 2024.

   This is the paper's core risk-quantification arithmetic: for each candidate risk
   factor, aggregate observed severe-infection hospitalisations, the population-based
   expected count, and person-years, then compute:
     - SMR (standardised morbidity/hospitalisation ratio) = observed / expected
     - AER (absolute excess risk per 10,000 person-years) = (observed-expected)/py*10000
     - Exact Poisson (chi-square-based) 95% confidence intervals for both
   The original runs this over ~25 risk factors at once via a multi-way
   PROC SUMMARY CLASS/TYPES/OUTPUT and a SELECT/WHEN block that figures out, row by
   row, which single risk factor that aggregate row belongs to. This bundle reproduces
   that structure unmodified, scaled down to two risk factors (sex, typeg) against a
   small mock aggregated table standing in for "shr_aer_inf_grave" (built upstream in
   the original program from the protected FCCSS/BCCSS patient extracts). */

/*One-row-per-patient aggregate: hospi_inf_grave = observed severe-infection
  hospitalisations, nb_hospi_mean = population-expected count, py = person-years.
  sex and typeg are two of the ~25 candidate risk factors from &liste_variables. */
data shr_aer_inf_grave;
	length id 8 sex 8 typeg 8 hospi_inf_grave 8 nb_hospi_mean 8 py 8;
	input id sex typeg hospi_inf_grave nb_hospi_mean py;
	datalines;
1 1 1 1 0.31 12.4
2 1 2 0 0.22  9.8
3 2 1 0 0.18  8.1
4 2 2 1 0.29 11.0
5 1 1 0 0.20  9.5
6 2 2 0 0.15  7.4
7 1 2 1 0.27 10.6
8 2 1 0 0.19  8.9
;
run;

/*6.3.Aggregate data*/
proc summary data=shr_aer_inf_grave; var hospi_inf_grave nb_hospi_mean py;
class sex typeg;
;
types () sex typeg
;
output out=inf_grave(rename=(_freq_=l)) sum(hospi_inf_grave nb_hospi_mean py)= hospi_inf_grave nb_hospi_mean py;
run;

/*Data management of aggregate data*/
data SHR_inf_grave; set inf_grave;
format classlevel1 $50.;
length classlevel1 $50.;
length classvar1 $100.;
length effect $30.;
select;
    when (not missing(typeg)) do;
      effect = "typeg";
      classvar1 = "FPN";
	  classlevel1 = typeg;
    end;
	when (not missing(sex)) do;
       effect = "sex";
     classvar1 = "Sex";
	  classlevel1 = sex;
    end;
	otherwise do;
      effect = "Overall";
      classvar1 = "Overall";
    end;
  end;
/*Calculation*/
d = hospi_inf_grave;
Nexp = nb_hospi_mean;
AER1 = ((d-Nexp)/py)*10000 ; /* AER = (observed death - expected deaths)/py*10000*/
SMR = d / Nexp;  /* SMR = observed death / expected deaths */
lower_d = quantile('CHISQUARE', 0.025, 2*d) / 2;
upper_d = quantile('CHISQUARE', 0.975, 2*(d+1)) / 2;
SMR_lower = lower_d / Nexp;
SMR_upper = upper_d / Nexp;
AER_lower = ((lower_d-Nexp)/py)*10000 ;
AER_upper = ((upper_d-Nexp)/py)*10000 ;
LL1=round(SMR_lower,0.01);
UP1=round(SMR_upper,0.01);
LL2=round(AER_lower,0.01);
UP2=round(AER_upper,0.01);
E=round(Nexp,0.1);
SMR1=round(SMR,0.01);
AER2=round(AER1,0.01);
py=round(py,0.1);
RHR=cat(SMR1," ","(",LL1,"-",UP1,")");
AER=cat(AER2," ","(",LL2,"-",UP2,")");
OE=cat(d,"/",E);

label RHR = "RHR (95% CI)"
      AER = "AER (95% CI)"
	  OE = "O/E"
	  py = "PY";
run;

data SHR_inf_grave(keep=classvar1 effect classlevel1 py OE RhR AER ); set SHR_inf_grave;  run;
data SHR_inf_grave;
retain effect classvar1 classlevel1 py OE RHR AER; set SHR_inf_grave;
run;

proc print data=SHR_inf_grave noobs label;
run;
