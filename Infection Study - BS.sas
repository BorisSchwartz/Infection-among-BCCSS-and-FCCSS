
/************************************************************************************************/
/*																							   	*/
/* SAS program created by Boris Schwartz, CESP-INSERM-U1018 Radiation Epidemiology Team - 2024 	*/
/*																							   	*/
/************************************************************************************************/
/*																							   	*/
/*			Late Infection-Related Risk among Childhood Solid Cancer Survivors 					*/
/*	A binational study from the French and British Childhood Cancer Survivor Studies			*/
/*																								*/
/************************************************************************************************/

/*1.Set up*/
/*Options to speed up processing time*/ 
OPTIONS nofmterr nonotes nosource error=0;
/*Processing time*/
%put Starting compilation of DATA step: %QSYSFUNC(DATETIME(), DATETIME20.3);
%let startTime=%QSYSFUNC(DATETIME());

/*date*/
%let year=%sysfunc(year(%sysfunc(date())));
%let month=%sysfunc(putn(%sysfunc(month(%sysfunc(date()))), Z2.0));
%let day=%sysfunc(putn(%sysfunc(day(%sysfunc(date()))), Z2.0));
%let date=&year.&month.&day;

/*path*/	/*CHANGE PATHS*/
%let path=C:\Users\b_schwartz\Desktop\Infection_study;

/*databases*/
%let path=C:\Users\b_schwartz\Desktop\Divers\Charlotte\29_01_2024\Immuno K2 Sein Thyr\new codes;
/*outputs*/
%let path_out=C:\Users\b_schwartz\Desktop\Divers\Charlotte\29_01_2024\sorties;
/*scripts*/
%let path_script=C:\Users\b_schwartz\Desktop\Divers\Charlotte\29_01_2024\new codes pgm;

/*Libraries*/ 	/*CHANGE LIBRARIES*/
/*Main lib.*/
libname inf "C:\Users\b_schwartz\Desktop\Divers\Charlotte\29_01_2024\Immuno K2 Sein Thyr\new codes";
/*Library for statistical results*/
libname RR "C:\Users\b_schwartz\Desktop\Divers\Charlotte\29_01_2024\Immuno K2 Sein Thyr\new codes\RR";
libname multi "C:\Users\b_schwartz\Desktop\Divers\Charlotte\29_01_2024\Immuno K2 Sein Thyr\new codes\multi";
libname SHR "C:\Users\b_schwartz\Desktop\Divers\Charlotte\29_01_2024\Immuno K2 Sein Thyr\new codes\RR";

/*Load script for Sankey plot*/
%include "&path_script\rawtosankey.sas";
%include "&path_script\sankey.sas";
%include "&path_script\sankeybarchart.sas";

/*Format*/
proc format  ;
value TYPE      10='Néphroblastome  '
                11='Autres tumeurs du rein '
                12='Tumeurs du rein sans précision'
                20='Neuroblastome'
                21='Ganglioneuroblastome'
                22='Autre tumeur du syst nerv sympat'
                30='Hodgkin'
                31='LMNH'
                40='Rhabdomyosarcome'
                41='Autres Tum embryo parties molles'
                42='Autres sarcome mous type adulte '
                50='Ewing '
                51='Ostéosarcome '
                52='Autres OS'
                60='Medulloblastome'
                61='Autres T cérébrales'
                62='T cérébrales sans histo'
                70='Séminome'
                71='Tératome malin'
                72='Autres Tum gonades'
                80='Thyroïde'
                90='Retinoblastome Unilatéral'
                91='Retinoblastome Bilatéral'
                00='Autres'
                                01='Autres : cortico-surrénalome'
                                02='Autres : UCNT'
                                05='Autres : foie'
                06='Autres : Autres carcinomes' ;

value TYPEG     1='Nephroblastome  '
                2='Neuroblastome'
                3='Lymphome'
                30='Hodgkin'
                31='LMNH'
                4='Tissus mous'
                5='Os '
                6='Cerveau'
                7='Gonades'
                8='Thyroïde'
                9='Retinoblastome'
                0='Autres' 
				99='Overall';
run;
/*2.Code list for infections (ICD-9-10)*/
%global liste_toute_bact4;
%global liste_toute_bact3;
%let liste_toute_bact4 =  	("A010","A021","A022","A170","A321","A327","A420","A427","A430","A481","A483","A488","A514","A544","A548",
							"G042","I330","J201","J202",
							"M010","M013","M902","N118","N119","N390","R572","R578",
							"0020","0031","0032","0130","3207","0270","0391","0399","4828","0408","0918","0985","0988","4660","5750","5751","7110","7114","7300","7301","7302","7307","5838","5845",
							"5900","5908","5990","7855")	;
%let liste_toute_bact3 =	("A39","A40","A41","A49","B95","B96","G00","G01","G02","J13","J14","J15","J16","J17","J18","J85","K65","K81","M00","M86","N10","O85","R65","036","038","041",
							"320","321","323","421","480","481","482","483","484","485","486","513","567","670");

%global liste_sepsis4;
%global liste_sepsis3;
%let liste_sepsis4 =  	("A021","A327","A427","0031","0270","0391","0399");
%let liste_sepsis3 =	("A40","A41","O85","038","670");

%global liste_choc_sris4;
%global liste_choc_sris3;
%let liste_choc_sris4 =  	("A483","R572","R578","0408","7855");
%let liste_choc_sris3 =		("R65");

%global liste_encaps4;
%global liste_encaps3;
%let liste_encaps4 =  	("A390","A395","A399","1394","A391","G001","A403","M001","A413","J201","G000","B963","A492");
%let liste_encaps3 =	("J13","J14");

%global liste_severe4;
%global liste_severe3;
%let liste_severe4 =  	("A483","R572","R578","0408","7855","A390","A395","A399","1394","A391","G001","A403","M001","A413","J201","G000","B963","A492");
%let liste_severe3 =	("J13","J14","R65");

%global liste_low4;
%global liste_low3;
%let liste_low4 =  	("A010","A022","A170","A321","A420","A430","A481","A488","A514","A544","A548",
					"G042","I330","J201","J202",
					"M010","M013","M902","N118","N119","N390",
					"0020","0032","0130","3207","4828","0918","0985","0988","4660","5750","5751","7110","7114","7300","7301","7302","7307","5838","5845",
					"5900","5908","5990")	;
%let liste_low3 =	("A39","A49","B95","B96","G00","G01","G02","J13","J14","J15","J16","J17","J18","J85","K65","K81","M00","M86","N10","036","041",
					"320","321","323","421","480","481","482","483","484","485","486","513","567");

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
%if &g=1 %then %do;
	axis3 label=("Effectif");
	axis4 label=(&nom);
	proc gchart data=&tab;
	vbar &var / discrete percent raxis=axis3 maxis=axis4;
	run;
	quit;
%end;
%mend;

/*4.Sankey diagram */
/*Transpose individual data to align infection dates*/
%tri(inf.t_all_inf, id); /*toute_bact: individual base with one row per patient for all infections, with 187 columns of potential infection records (date, type, code)*/
proc transpose data = inf.t_all_inf out=tran;/*prefix t_ for test database*/
var dte_toute_bact_1-dte_toute_bact_187;
by id;
run;
proc transpose data = inf.t_all_inf out=tran3 prefix = type_code;
var type_code_1-type_code_187;
by id;
run;
proc transpose data = inf.t_all_inf out=tran5 prefix = code;
var code_1-code_187;
by id;
run;
data tranF;
merge tran tran3 tran5;
format col1 ddmmyy10.;
if col1 = . then delete;
run;
data tranF;
set tranF;
by id;
format COL1 ddmmyy10.;
lag_date = lag(col1);
/*Outcome dates are identified and accounted if at least 14 days have elapsed since the previous infection.*/
if (col1-lag(col1) le 14) and first.id ne 1 then supp = 1;
rename col1 = dte_toute_bact;
run;
data tranF2;
set tranF;
by id dte_toute_bact;
retain cpt 0;
if first.id  then cpt = 1;
else do;
	if supp ne 1 then do;
		cpt = cpt+1;
	end;
end;
run;
data tranF3;
set tranF2;
format inf $20. inf_num 8.;
if substr(code1, 1, 4) in (&liste_sepsis4) or substr(code1, 1, 3) in (&liste_sepsis3) 		then inf = "sepsis";
if substr(code1, 1, 4) in (&liste_choc_sris4) or substr(code1, 1, 3) in (&liste_choc_sris3) then inf = "choc_sris";
if substr(code1, 1, 4) in (&liste_low4) or substr(code1, 1, 3) in (&liste_low3) 			then inf = "low";
if inf = "low" 			then inf_num = 1;
if inf = "sepsis" 		then inf_num = 2;
if inf = "choc_sris" 	then inf_num = 3;
run;
%tri(tranF3, id cpt);
/*Diagnosis according to type of diagnosis*/
proc freq data = tranF3;
table inf*type_code1/missing;
run;
/*Priority is given to the most serious infection in the event of a tie (i.e. less than 14 days): choc_sris > sepsis > low*/
proc sql;
create table test as
select id, cpt, max(inf_num) as inf_ok, min(dte_toute_bact) as dte_inf_ok format ddmmyy10.
from tranF3
group by id, cpt;
quit;
data data_sankey;
set test;
format infection $15.;
if inf_ok = 1 then infection = "Low grade";
if inf_ok = 2 then infection = "Sepsis";
if inf_ok = 3 then infection = "Severe";
run;

/*Example Sankey diagram in the cohort (can be declined in subgroup if needed by adding some selection lines in the program)*/
%macro sankey_inf(n);
/*Sankey formatting : one line per individual per infection, so here 15 lines per patient (max inf. = 15)*/
data sankey_inf;
set data_sankey;
if supp = 1 then delete;
keep id infection cpt;
run;
/*Patient identification and reference database*/
data sankey_pat;
set sankey_inf;
by id;
if first.id;
run;
/*Artificial data (no infection)*/
%macro artif();
%do i = 2 %to 15;
	data artif_&i;
	set sankey_pat;
	cpt = &i;
	infection = "No infection";
	run;
%end;
data artif_all;
set artif_2-artif_15;
run;
%tri(artif_all, id cpt);
%mend;
%artif();
/*Artificial data are selected for missing real lines*/
proc sql noprint;
create table artif as
select b.*
from sankey_inf as a right join artif_all as b
on a.id = b.id and a.cpt = b.cpt
where a.cpt is null;
quit;
/*Merging with real data*/
data base_sankey;
set sankey_inf artif;
run;
%tri(base_sankey nodupkey, id cpt);
/*keeping &n first infections*/
data base_sankey_&n;
set base_sankey;
if cpt gt &n then delete;
run;
/*Patients from analysis database only*/
proc sql;
create table base_sankey_&n as
select *
from base_sankey_&n
where id in(select id from inf.test_analyse_toute_bact);/*Add selection criteria here if needed*/
quit;
/*Analysis of &n first infections*/
title "Sankey diagram - Infection";
footnote ;
/****************************/
/*Sankey programs to download from: SAS Sankey macro created by Shane Rosanbalm of Rho, Inc. 2015 https://github.com/RhoInc/sas-sankeybarchart*/
/****************************/
%sankeybarchart(data=base_sankey_&n
   ,subject=id
   ,yvar=infection
   ,xvar=cpt
   ,yvarord=%quote(Low grade, Sepsis, Severe, No infection)
   ,colorlist=VLIGB BIOY BIYG LIGGR
   );
%mend;
ods rtf file ="&path_out.\&date._Sankey_infections.rtf";
/*3, 5, 10 and 15th first infections, resp.*/
%sankey_inf(3);
%sankey_inf(5);
%sankey_inf(10);
%sankey_inf(15);
ods rtf close;

/*5.Calculation of PY : from a base with one line per patient to a base with one line per patient per year*/
/*Calculation of PY for each year. Event stopping follow-up = death or infection related hospitalization */ 
/*Example for inf_grave endpoint*/
%macro calcul_py(base,ev);/*base: individual base; ev: indicator of event of interest*/
%macro py_card(an);/*an: year*/
/*calcul des PY*/
data py_&an;
set inf.test_analyse_&base;
format entree ddmmyy10. sortie ddmmyy10.;
annee_deces = year(date_deces);
if &ev ne 1 then dte_&ev = .;
annee_&ev = year(dte_&ev);
/*Traitement des patients ayant eu une hospitalisation &ev entre 2006 et 2018*/
if &ev = 1 then do;
	if annee_&ev lt &an then delete; 
	if annee_&ev = &an then do;
		entree = mdy(01,01,&an);
		/*Entry at 04/01/1997 for UK so PY recalculated for this year (remove 3 months i.e. 3/12)*/
		if &an = 1997 then do;
			entree = mdy(04,01,&an);
		end;
		sortie = dte_&ev;
		py_&an = (sortie - entree)/365.25;
	end;
	if annee_&ev gt &an then do;
		entree = mdy(01,01,&an);
		/*Entry at 04/01/1997 for UK so PY recalculated for this year*/
		if &an = 1997 then do;
			entree = mdy(04,01,&an);
		end;
		sortie = mdy(12,31,&an);
		py_&an = 1;
	end;
end;
if deces = 1 then do;
	if annee_deces lt &an then delete;
	if annee_deces = &an then do;
		/*If the patient died during the year, he is only counted if he was not hospitalized during the year. If this is the case, the patient has already been counted in the previous step*/
		if &ev ne 1 then do; 
			entree = mdy(01,01,&an);
			/*Entry at 04/01/1997 for UK so PY recalculated for this year*/
			if &an = 1997 then do;
				entree = mdy(04,01,&an);
			end;
			sortie = date_deces;
			py_&an = (sortie - entree)/365.25;
		end;
	end;
	if annee_deces gt &an then do;
		if annee_&ev gt &an then do;
			entree = mdy(01,01,&an);
			/*Entry at 04/01/1997 for UK so PY recalculated for this year*/
			if &an = 1997 then do;
				entree = mdy(04,01,&an);
			end;
			sortie = mdy(12,31,&an);
			py_&an = 1;
		end;
		if annee_&ev = . then do;
			entree = mdy(01,01,&an);
			/*Entry at 04/01/1997 for UK so PY recalculated for this year*/
			if &an = 1997 then do;
				entree = mdy(04,01,&an);
			end;
			sortie = mdy(12,31,&an);
			py_&an = 1;
		end;
	end;
end;
if deces = 0 and &ev = 0 then do;
	entree = mdy(01,01,&an);
	/*Entry at 04/01/1997 for UK so PY recalculated for this year*/
	if &an = 1997 then do;
		entree = mdy(04,01,&an);
	end;
	sortie = mdy(12,31,&an);
	py_&an = 1;
end;
/*age_pat : age in the middle of the year*/
age_pat = floor((mdy(06,30,&an) - ddn)/365.25);
annee = &an;
if &ev = . then &ev = 0;
run;
%mend;
%py_card(1997);%py_card(1998);%py_card(1999);%py_card(2000);%py_card(2001);%py_card(2002);%py_card(2003);%py_card(2004);%py_card(2004);%py_card(2005);
%py_card(2006);%py_card(2007);%py_card(2008);%py_card(2009);%py_card(2010);%py_card(2011);%py_card(2012);
%py_card(2013);%py_card(2014);%py_card(2015);%py_card(2016);%py_card(2017);%py_card(2018);
/*data fusion: one line per individual per year. Age evolves with year*/
data t_individ1y_&ev;/*prefix t_ for test database*/
retain py py_1997 py_1998 py_1999 py_2000 py_2001 py_2002 py_2003 py_2004 py_2005
	py_2006 py_2007 py_2008 py_2009 py_2010 py_2011 py_2012 py_2013 py_2014 py_2015 py_2016 py_2017 py_2018;
set py_1997 py_1998 py_1999 py_2000 py_2001 py_2002 py_2003 py_2004 py_2005 
	py_2006 py_2007 py_2008 py_2009 py_2010 py_2011 py_2012 py_2013 py_2014 py_2015 py_2016 py_2017 py_2018;
py = sum(py_1997,py_1998,py_1999,py_2000,py_2001,py_2002,py_2003,py_2004,py_2005,py_2006,py_2007,py_2008,py_2009,py_2010,py_2011,py_2012,py_2013,py_2014,py_2015,py_2016,py_2017,py_2018);
/*Drop temporary variables*/
drop py_1997-py_2018;
run;
/*deleting temporary tables*/
%supp(py_1997-py_2018);
%mend;
/*Macro running for each endpoint*/
/*severe (called inf_grave or choc_sris)*/
%calcul_py(choc_sris,inf_grave);

/*6.SHR/AER Calculation: EXAMPLE FOR INF_GRAVE (To adapt and reproduce for each endpoint)*/
%global liste_variables;
%let liste_variables = 	sex c_age_diag cl_rt_rate cl_rt_thymus bipulm_bs cl_rt_pituitary cl_rt_foie splenectomy H R T autre pb_k pulm_choc_sris coeur_choc_sris hepat_choc_sris greffe_choc_sris Insuff_immun;

/*6.1.FCCSS*/
/*Import rates*/
 data TAUX_fr    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile "&path.\taux_esnd.txt" delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat POP_ESND best32. ;
       informat SEX best32. ;
       informat AGE best32. ;
       informat ANNEE best32. ;
       informat NB_TOUTE_BACT best32. ;
       informat NB_LOW best32. ;
       informat NB_SEPSIS best32. ;
       informat NB_INF_GRAVE best32. ;
       informat TX_TOUTE_BACT nlnum32. ;
       informat TX_LOW nlnum32. ;
       informat TX_SEPSIS nlnum32. ;
       informat TX_INF_GRAVE nlnum32. ;
       format POP_ESND best12. ;
       format SEX best12. ;
       format AGE best12. ;
       format ANNEE best12. ;
       format NB_TOUTE_BACT best12. ;
       format NB_LOW best12. ;
       format NB_SEPSIS best12. ;
       format NB_INF_GRAVE best12. ;
       format TX_TOUTE_BACT nlnum12. ;
       format TX_LOW nlnum12. ;
       format TX_SEPSIS nlnum12. ;
       format TX_INF_GRAVE nlnum12. ;
    input
                POP_ESND
                SEX
                AGE
                ANNEE
                NB_TOUTE_BACT
                NB_LOW
                NB_SEPSIS
                NB_INF_GRAVE
                TX_TOUTE_BACT
                TX_LOW
                TX_SEPSIS
                TX_INF_GRAVE
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;

/*The corresponding reference rate is added to each line according to year, sex and age*/
data fr;
set t_individ1y_inf_grave;
if country = "FRA";
/*Age in the middle of the year*/
age_pat = round((mdy(06,30,annee)-date_nais)/365.25);
keep id sex age_pat annee py choc_sris;
run;
%tri(fr, sex age_pat annee);
data taux;
set taux_fr;
age_pat = age;
/*Inf grave*/
keep sex age_pat annee tx_inf_grave;
run;
%tri(taux, sex age_pat annee);
/*Merging with reference rates and calculation of expected numbers*/
data fr_rate;
merge fr taux;
by sex age_pat annee;
/*expected (called "nb_hospi_mean")*/
nb_hospi_mean = py*tx_inf_grave;
/*observed*/
hospi_inf_grave = choc_sris;
/*Removal of surplus lines*/
if py = . then delete;
run;
/*Aggregated database for patient follow-up: one line per patient*/
proc sql;
create table fr_ag as
select distinct id, sum(hospi_inf_grave) as hospi_inf_grave, sum(nb_hospi_mean) as nb_hospi_mean, sum(py) as py
from fr_rate
group by id;
quit;
/*Definition of formats*/
data fr_ag;
set fr_ag;
format c_age_fin $20. c_fup_fin $20. c_age_2006 $20. /*c_age_diag $20.*/ ttt_era $20. trt $20.;
run;
%tri(t_individ1y_inf_grave, id);
data last_ligne;
set t_individ1y_inf_grave;
by id;
if last.id;
if country = "FRA";
/*drop py in this base to avoid overwriting data*/
drop py;
run;
/* ADD TREATMENT DATA*/
proc sql;
create table fr_ag as
select *
from fr_ag as a left join last_ligne as b
on a.id = b.id;
quit;

/*6.2. BCCSS*/
/*Import rates*/
data WORK.Taux_uk;
%let _EFIERR_ = 0; /* set the ERROR detection macro variable */
infile "&path.\index_infections_2.csv" delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
   informat sex best32. ;
   informat age best32. ;
   informat calyr best32. ;
   informat risk_score best32. ;
   informat popn best32. ;
   informat cases best32. ;
   informat rate best32. ;
   format sex best12. ;
   format age best12. ;
   format calyr best12. ;
   format risk_score best12. ;
   format popn best12. ;
   format cases best12. ;
   format rate best12. ;
input
            sex
            age
            calyr
            risk_score
            popn
            cases
            rate
;
if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
run;
/*The corresponding reference rate is added to each line according to year, sex and age.*/
data uk;
set t_individ1y_inf_grave;
if country = "ENG";
/*Age in the middle of the year*/
age_pat = round((mdy(06,30,annee)-date_nais)/365.25);
keep id sex age_pat annee py choc_sris;
run;
%tri(uk, sex age_pat annee);
data taux;
set taux_uk;
age_pat = age;
annee = calyr;
/*Inf grave: score=3 (low-grade=1, sepsis=2)*/
where risk_score = 3;
keep sex age_pat annee rate;
run;
%tri(taux, sex age_pat annee);
/*Merging with reference rates and calculation of expected numbers*/
data uk_rate;
merge uk taux;
by sex age_pat annee;
/*Calculation of py since april 1997 (i.e. on remove 3 months)*/
if annee = 2017 then py = py - (3/12);
/*expected (called "nb_hospi_mean" to harmonise with FCCSS data)*/
nb_hospi_mean = py*rate/100000;
/*observed*/
hospi_inf_grave = choc_sris;
/*Removal of surplus lines*/
if py = . then delete;
run;

/*Aggregated database for patient follow-up: one line per patient*/
proc sql;
create table uk_ag as
select distinct id, sum(hospi_inf_grave) as hospi_inf_grave, sum(nb_hospi_mean) as nb_hospi_mean, sum(py) as py
from uk_rate
group by id;
quit;
/*Definition of formats*/
data uk_ag;
set uk_ag;
format c_age_fin $20. c_fup_fin $20. c_age_2006 $20. /*c_age_diag $20.*/ ttt_era $20. trt $20.;
run;
%tri(t_individ1y_inf_grave, id);
data last_ligne;
set t_individ1y_inf_grave;
by id;
if last.id;
if country = "ENG";
/*drop py in this base to avoid overwriting data*/
drop py;
run;
/* ADD TREATMENT DATA*/
proc sql;
create table uk_ag as
select *
from uk_ag as a left join last_ligne as b
on a.id = b.id;
quit;
/*Merging FCCSS and BCCSS data*/
data shr_aer_inf_grave;
set fr_ag uk_ag;
keep &liste_variables deces debut sortie hospi_inf_grave /*=obs*/ nb_hospi_mean/*=exp*/ py id;
run;

/*6.3.Aggregate data*/
proc summary data=shr_aer_inf_grave; var hospi_inf_grave nb_hospi_mean py;/*******************************/
class 					&liste_variables;
 ;
types () 				&liste_variables
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
	when (not missing(chimiotherapie)) do;
      effect = "chimiotherapie";
      classvar1 = "Chemotherapy";
	  classlevel1 = chimiotherapie;
    end;
    when (not missing(typeg)) do;
      effect = "typeg";
      classvar1 = "FPN";
	  classlevel1 =typeg;
    end;
	when (not missing(ttt_era)) do;
       effect = "ttt_era";
     classvar1 = "Treatment era";
	  classlevel1 = ttt_era;
    end;
    when (not missing(c_age_diag)) do;
      effect = "c_age_diag";
      classvar1 = "Age at diagnosis";
      classlevel1 =  c_age_diag;
    end;
	when (not missing(sex)) do; 
       effect = "sex";
     classvar1 = "Sex";
	  classlevel1 = sex;
    end;
	when (not missing(trt)) do;
      effect = "trt";
      classvar1 = "FPN treatment";
	  classlevel1 = trt;
    end;
	when (not missing(splenectomie)) do;
	  effect = "splenectomie";
      classvar1 = "Splenectomy";
	  classlevel1 = splenectomie;
    end;
	when (not missing(cl_anthra)) do;
      effect = "cl_anthra";
      classvar1 = "Anthracyclines, in mg/m²";
	  classlevel1 = cl_anthra;
    end;
    when (not missing(anthra)) do;
      effect = "anthra";
      classvar1 = "Anthracyclines administration";
	  classlevel1 = anthra;
	end;
	when (not missing(cl_rt_rate)) do;
      effect = "cl_rt_rate";
      classvar1 = "Spleen radiation dose, in Gy";
	  classlevel1 = cl_rt_rate ;
    end;
	when (not missing(cl_rt_thymus)) do;
      effect = "cl_rt_thymus";
      classvar1 = "Thymus radiation dose, in Gy";
	  classlevel1 = cl_rt_thymus ;
    end;
		when (not missing(bipulm_bs )) do;
      effect = "bipulm_bs";
      classvar1 = "Bipulmonary field";
	  classlevel1 = bipulm_bs ;
    end;
		when (not missing(cl_rt_pituitary)) do;
      effect = "cl_rt_pituitary";
      classvar1 = "Pituitary radiation dose, in Gy";
	  classlevel1 = cl_rt_pituitary ;
    end;
		when (not missing(cl_rt_foie)) do;
      effect = "cl_rt_foie";
      classvar1 = "Liver radiation dose, in Gy";
	  classlevel1 = cl_rt_foie ;
    end;
		when (not missing(splenectomy)) do;
      effect = "splenectomy";
      classvar1 = "Splenectomy";
	  classlevel1 = splenectomy  ;
    end;
		when (not missing(H)) do;
      effect = "H";
      classvar1 = "H = rth with D80% pituitary  > 35 Gy";
	  classlevel1 = H ;
    end;
		when (not missing(R)) do;
      effect = "R";
      classvar1 = "S = rth with V5Gy spleen > 80%";
	  classlevel1 = R;
    end;
		when (not missing(T)) do;
      effect = "T";
      classvar1 = "Ty = rth with D60% thymus > 20 Gy and age < 3 years";
	  classlevel1 = T;
    end;
		when (not missing(autre)) do;
      effect = "autre";
      classvar1 = "Other treatment of rth";
	  classlevel1 = autre;
    end;
		when (not missing(pb_k)) do;
      effect = "pb_k";
      classvar1 = "Second cancer";
	  classlevel1 = pb_k;
    end;
		when (not missing(pulm_choc_sris)) do;
      effect = "pulm_choc_sris";
      classvar1 = "Pulmonary sequelae";
	  classlevel1 = pulm_choc_sris;
    end;
		when (not missing(coeur_choc_sris)) do;
      effect = "coeur_choc_sris";
      classvar1 = "Cardiac sequelae";
	  classlevel1 = coeur_choc_sris ;
    end;
		when (not missing(hepat_choc_sris)) do;
      effect = "hepat_choc_sris";
      classvar1 = "Liver failure";
	  classlevel1 = hepat_choc_sris ;
    end;
		when (not missing(greffe_choc_sris)) do;
      effect = "greffe_choc_sris";
      classvar1 = "Solid organ transplant ";
	  classlevel1 = greffe_choc_sris ;
    end;
		when (not missing(Insuff_immun)) do;
      effect = "Insuff_immun";
      classvar1 = "Immunodeficiency";
	  classlevel1 = Insuff_immun;
    end;
	when (not missing(vaccin)) do;
      effect = "vaccin";
      classvar1 = "Vaccination";
      classlevel1 =  vaccin;
    end;
	/*Add variables if needed*/
	/******/
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
/*deleting temporary base*/
%supp(inf_grave);

/*6.4.Pvalues of RHR */
%macro RHR(var); 
data Temp; set shr_aer_inf_grave; 
pyr = py;
if pyr ne 1 then w = 1;
ln_y = log(pyr);
cases = hospi_inf_grave;
Nexp =  nb_hospi_mean;
run;
proc sort data=Temp; by &var ;run;
proc means data=temp sum; 
var cases Nexp pyr; 
ods output summary=RHR_all; 
by &var ;
run;
/*Creation of indicators*/
data RHR_all; 
set RHR_all; 
excess=(cases_Sum-Nexp_Sum);
pyr1=log(Pyr_Sum); 
/*adding a negligeable value to avoid division by zero*/
Expected1=log(Nexp_Sum+0.000000001 );
SMR = (cases_Sum/Nexp_Sum);
RHR = ((cases_Sum-Nexp_Sum)/pyr_sum)*10000;
run;
ods output Tests3=RHR_&var;
proc glimmix data=RHR_all; 
class &var;
/*Calculation of global p values for SHR*/
model cases_Sum =&var / dist=poisson offset=Expected1 ddfm=none s cl htype=3; run;
data RHR_&var;
format effect $100. ;
set RHR_&var;
effect = "&var";
run;
%supp(rhr_all);
%mend;
/*Base SHR*/
%macro base_RHR();
/*List of variables for output table (in correct format)*/
%macro RHR_liste(liste); /*Macro for automatically selecting variables in outputs*/
	%global liste_var_RHR;
	/*Nb of variables*/
	%let nb_elem = %sysevalf(1+%sysfunc(count((&liste),%str( ))));		/*Number of elements in the list (if values are added to/removed from lists, the program still runs)*/
	%let liste_var_RHR = ;
	%do i = 1 %to &nb_elem;					/*Loop over the number of elements in the list*/
		%let variab = %scan(&liste,&i,’ ’); 
/**/
		%RHR(&variab);
/**/
		%let variab_var = RHR_%scan(&liste,&i,’ ’);  /*We select the ith variable in the list*/
		%let liste_var_RHR = &liste_var_RHR &variab_var;
	%end;
%mend;
%RHR_liste(&liste_variables);
data pv_RHR (keep= effect  probF); format effect $100.  probF pvalue6.4;
set &liste_var_RHR
; 
run;
%supp(&liste_var_RHR);
data pv_RHR (rename=probF=p_RHR); set pv_RHR; run;
%tri(pv_RHR, effect);
%svg(pv_RHR,ShR,pval_RHR_mean_inf_grave);
%supp(pv_RHR);
%mend;
%base_RHR();

/*6.5.Pvalues of AER */
%macro AER(var); 
data Temp; set shr_aer_inf_grave; 
pyr = py;
if pyr ne 1 then w = 1;
ln_y = log(pyr);
cases = hospi_inf_grave;
Nexp =  nb_hospi_mean;
run;
proc sort data=Temp; by &var ;run;
proc means data=temp sum; 
var cases Nexp pyr; 
ods output summary=aer_all; 
by &var ;
run;
/*Creation of indicators*/
data aer_all; set aer_all; 
excess=(cases_Sum-Nexp_Sum);
pyr1=log(Pyr_Sum); 
/*adding a negligeable value to avoid division by zero*/
Expected1=log(Nexp_Sum+0.000000001 );
SMR = (cases_Sum/Nexp_Sum);
AER = ((cases_Sum-Nexp_Sum)/pyr_sum)*1000;
run;
ods output Tests3=aer_&var;
proc glimmix data=aer_all; 
class &var;
/*Calculation of global p values for AER*/
model excess =&var / dist=poisson offset=pyr1 ddfm=none s cl htype=3; run;
data aer_&var;
format effect $100. ;
set aer_&var;
effect = "&var";
run;
%supp(aer_all);
%mend;
/*Base AER*/
%macro base_aer();
/*List of variables for output table (in correct format)*/
%macro aer_liste(liste); /*Macro for automatically selecting variables in outputs*/
	%global liste_var_aer;
	/*Nb of variables*/
	%let nb_elem = %sysevalf(1+%sysfunc(count((&liste),%str( ))));		/*Number of elements in the list (if values are added to/removed from lists, the program still runs)*/
	%let liste_var_aer = ;
	%do i = 1 %to &nb_elem;					/*Loop over the number of elements in the list*/
		%let variab = %scan(&liste,&i,’ ’); 
/**/
		%AER(&variab);
/**/
		%let variab_var = aer_%scan(&liste,&i,’ ’);  /*We select the ith variable in the list*/
		%let liste_var_aer = &liste_var_aer &variab_var;
	%end;
%mend;
%aer_liste(&liste_variables);
data pv_AER (keep= effect  probF); format effect $100.  probF pvalue6.4;
set &liste_var_aer
; 
run;
%supp(&liste_var_aer);
data pv_AER (rename=probF=p_AER); set pv_AER; run;
%tri(pv_AER, effect);
%svg(pv_AER,ShR,pval_AER_mean_inf_grave);
%supp(pv_AER);
%mend;
%base_AER();
%supp(temp);

/*6.6.final table formatting*/
%macro tab_uni();
/*Correspondence table effect/classvar1*/
data corresp;
set shr_inf_grave;
keep classvar1 effect;
run;
%tri(corresp, effect);
%tri(shr_inf_grave, effect);
%tri(shr.pval_aer_mean_inf_grave, effect);
%tri(shr.pval_rhr_mean_inf_grave, effect);
/*Merge tables to add classvar1*/
data pval_aer2;
merge shr.pval_aer_mean_inf_grave corresp;
by effect;
if p_AER ne .;
run;
%tri(pval_aer2 nodupkey, effect p_AER);
data pval_RhR2;
merge shr.pval_RhR_mean_inf_grave corresp;
by effect;
if p_Rhr ne .;
run;
%tri(pval_RhR2 nodupkey, effect p_RhR);

data All_uni; 
format paer PVALUE6.4;
format pRhR PVALUE6.4;
merge pval_aer2 pval_RhR2 shr_inf_grave;
by effect;
if first.effect then paer = p_AER;
if first.effect then pRhR = p_RhR;
run;
%supp(pval_aer2 pval_RhR2);
/*******************Grouping results and pval in different tables********************/
data calcul (keep=Effect classvar1 classlevel1 tne pyr OE rhr AER); set All_uni; run;
data p_hetero(keep=effect prhr paer); set All_uni;run;
/* p hetero in line**********************************************************************************/
data pval_hetero; set p_hetero;if prhr=. then delete; 
classlevel1 = "p-heterogeneity";
prhr = round(prhr, 0.001);
paer = round(paer, 0.001);
/*num to char*/
rhr = input(prhr, $20.); 
AER = input(paer, $20.); 
if rhr = "0" 	then rhr = "<.001";
if AER = "0" 	then AER = "<.001";
drop prhr paer;
run;
data all; set calcul pval_hetero ; if classvar1 = "" then classvar1 = effect; run;
%tri(all, effect classlevel1);
data all2; set all;
	if effect = 'Overall' then ordre = 1;
	if effect = 'sex' then ordre = 2;
	if effect = 'typeg' then ordre = 3;
	if effect = 'c_age_diag' then ordre = 4;
	if effect = 'ttt_era' then ordre = 5;
	if effect = 'trt' then ordre = 6;
	if effect = 'chimiotherapie' then ordre = 7;
	if effect = 'anthra' then ordre = 8;
	if effect = 'cl_anthra' then ordre = 9;
	if effect = 'cl_rt_rate' then ordre = 10;
	if effect = 'cl_rt_thymus' then ordre = 11;
	if effect = 'bipulm_bs' then ordre = 12;
	if effect = 'cl_rt_pituitary' then ordre = 13;
	if effect = 'cl_rt_foie' then ordre = 14;
	if effect = 'splenectomy' then ordre = 15;
	if effect = 'H' then ordre = 16;
	if effect = 'R' then ordre = 17;
	if effect = 'T' then ordre = 18;
	if effect = 'autre' then ordre = 19;
	if effect = 'pb_k' then ordre = 20;
	if effect = 'pulm_choc_sris' then ordre = 21;
	if effect = 'coeur_choc_sris' then ordre = 22;
	if effect = 'hepat_choc_sris' then ordre = 23;
	if effect = 'greffe_choc_sris' then ordre = 24;
	if effect = 'Insuff_immun' then ordre = 25;
	if effect = 'Vaccin' then ordre = 26;
	/*Add variables if needed*/
run;
data all3; set all2; format level $100.; 
level = classlevel1;
run;
%tri(all3, classvar1);
data all4; set all3; by classvar1; if first.classvar1 ne 1 then classvar1 = ""; run;
proc sort data=all4 out= all5; by ordre classlevel1; run; 
data x;  set all5 (drop= classlevel1 ordre ); run;
data all6;  set x; format factor $100.; 
factor = classvar1;
run;
data x2;  set all6 (drop=effect classvar1); run;
data Table_1; retain factor level; set x2; run;
/**************************formatting and saving************/
options mprint mlogic;
ods rtf file= "&path_out.\table_SHR_AER_inf_grave_&date..rtf" 
style=statistical fontscale=75 ;
proc report data=Table_1
split='#' headline headskip
   style(report)=[rules       = groups
                  background  = white
                  bordercolor = white
				   borderwidth = .2cm
				  frame = hsides]
   style(header)=[background  = white
                  font_size   = 7pt
                  font_face   = 'Arial'
                  borderbottomcolor=black
				  bordertopcolor=black
				  bordertopwidth=.5pt
                  foreground=black
                  bordertopwidth=.5pt
                  just = c]
   style(column)=[font_size   = 7pt
                  font_face   = 'Arial'
                  cellwidth   = 2 cm
                  just        = l];
define rhr / 'shr (95% IC)'   center;
define AER / 'AER (95% IC)'   center;
define OE / 'O/E' center;
define level / 'Level' center;
define factor / 'Factor' center;
title1 "Table-XX - Observed and expected numbers of severe infection related hospitalization, Relative Hospitalization Ratios and Absolute Excess Risks per 10000 Person-Years";
run;
options orientation=portrait;
ods rtf close;
%mend;
%tab_uni();


/*7.Macro for relative risks from Poisson regression*/
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
proc means data=base sum; 
var Observed pyr; 
ods output summary=smr_all; 
by &var ;
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
data EMR_&var(keep=Effect &var modal n RR IC_RER pv); 
format effect $20. modal $200. pv pvalue6.4;
set pa; 
modal = &var;
if &var ne &r then do; 
	IC_RER=cat(RR," ","(",LL,"-",UL,")");
	n =_N_;
	pv = probt;
end; 
if Effect ne 'Intercept'; 
run;
%svg(EMR_&var,RR,RR_&nom._&var);
/********* RR CALCULATION Univariate - Global Pvalues*/
%macro RR_glob(base,var,r="r", nom = ., obs = .); /*r = reference modality if qualitative variable ;
											class = 0 if quantitative variable, 1 if qualitative variable*/
proc means data=base sum; 
var Observed pyr; 
ods output summary=smr_all; 
by &var ;
run;
data smr_all; 
set smr_all; 
pyr1=log(Pyr_Sum); 
run;
ods output Tests3=global_&var;
proc glimmix data=smr_all; 
class &var(ref=&r);
model observed_sum = &var / dist=poisson offset=PYR1 ddfm=none s cl htype=3  /*for global p-values*/;
run;
%svg(global_&var,RR,RR_glob_&nom._&var);
%mend;
%RR_glob(base= base, var = &var, r = &r, nom = &nom, obs = &obs);
%mend;

/*7.2.RR MULTIVARIATE*/
%macro RR_multi(base,nom, obs, model);/*base=base to use; nom=name of output base; obs = event; model=number of the model (for not erasing outputs)*/
proc sort data=&base; 
by  &&liste_variables_&model;
run;
data base1; 
set &base; 
Observed = &obs;
pyr = py;
run;
/***********************/
proc means data=base1 sum; 
var Observed pyr; 
ods output summary=smr_all; 
by &&liste_variables_&model;
run;
%supp(base1);
data smr_all; 
set smr_all; 
pyr1=log(Pyr_Sum); 
run;
ods output Tests3=global_multi;
proc glimmix data=smr_all; 
class &&liste_variables_ref_&model ;
model Observed_Sum = &&liste_variables_&model /  dist=poisson offset=pyr1 ddfm=none s cl htype=3  /*htype=3 --> for global pvalues*/  ;
run;
/************************/
ods output ParameterEstimates=aa;
proc glimmix data=smr_all; 
class &&liste_variables_ref_&model ;
model Observed_Sum = &&liste_variables_&model / dist=poisson offset=pyr1 ddfm=none s cl /*htype=3 --> for global pvalues*/  ;
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
data EMR_multi(keep=Effect &&liste_variables_&model modal n RR IC_RER pv); 
format effect $20. modal $200.;
set pa; 
/*Recovery of modalities*/
%let nb_elem = %sysevalf(1+%sysfunc(count((&&liste_variables_&model),%str( ))));
%do i = 1 %to &nb_elem;					/*Loop over the number of elements in the list*/
	%let v_indiv = %scan(&&liste_variables_&model,&i,’ ’);  /* select the ith variable in the list*/ 
	var_&i = put(&v_indiv, 20.);
	var_&i = STRIP(var_&i);
	if var_&i ne "_" and var_&i ne "" then modal = &v_indiv;
%end;
	IC_RER=cat(RR," ","(",LL,"-",UL,")");
	n =_N_;
	pv = probt;
if Effect ne 'Intercept'; 
modal = STRIP(modal);
run;
%svg(EMR_multi,RR, RR_&nom._&model._multi);
%svg(global_multi,RR, RR_glob_&nom._&model._multi);
%mend;

/*7.3.OUTPUT TABLES RR (pv. indiv + global) */
%macro fusion_RR(nom, model, obs);
/*Grouping of results*/
data EMR_all;
format effect $20. pv pvalue6.4;
set &liste_indiv;
run;
data EMR_pv_global (keep= effect probF);
format effect $20. ;
set &liste_glob;
run;
%tri(EMR_pv_global, effect);
%tri(EMR_all, effect IC_RER modal);
/*table global univariate*/
data EMR_all;
format pv_globale PVALUE6.4;
format modalite $200.;
merge EMR_all EMR_pv_global;
by effect;
if first.effect then pv_globale = probF; 
if n = . then IC_RER = " Ref.";
if IC_RER = "1.0 ( . - . )" 	then IC_RER = " Ref."; 
if IC_RER = "1.00 ( .  - .  )" 	then IC_RER = " Ref.";
modal = STRIP(modal);
modalite = modal;
run;
/*Add results from multivariate*/
data EMR_all_multi;
format effect $20. pv pvalue6.4;
set RR.RR_&nom._&model._multi;
run;
/*Grouping of results*/
data EMR_pv_global_multi (keep= effect probF);
format effect $20. ;
set RR.RR_glob_&nom._&model._multi;
run;
%tri(EMR_pv_global_multi, effect);
%tri(EMR_all_multi, effect IC_RER modal);
/*overall table multivariate*/
data EMR_all_multi;
format pv_globale PVALUE6.4;
format pv PVALUE6.4;
/*format treatment treatment.;*/
merge EMR_all_multi EMR_pv_global_multi;
by effect;
if first.effect then pv_globale = probF;
if n = . then IC_RER = "Ref.";
if IC_RER = "1.00 ( .  - .  )" then IC_RER = " Ref."; 
run;
%ql(emr_all_multi, effect);
data emr_all2;
format modalite $200.;
format pv_globale PVALUE6.4;
format pv PVALUE6.4;
set emr_all_multi;
if n = . then IC_RER = " Ref.";
if IC_RER = "1.00 ( .  - .  )" 	then IC_RER = " Ref."; 
if IC_RER = "1.0 ( . - . )" 	then IC_RER = " Ref.";
modalite = modal;
/*Change variable names to merge with univariate results*/
format pv_globale PVALUE6.4 pv_globale_multi PVALUE6.4;
/**/format pv PVALUE6.4 pv_multi PVALUE6.4;
IC_RER_multi = IC_RER;
pv_multi = pv;
pv_globale_multi = pv_globale;
drop IC_RER pv pv_globale;
run;
/*Merge univariate and multivariate results*/
%tri(emr_all, effect modalite);
%tri(emr_all2, effect modalite);
data EMR_fusion;
format pv_globale PVALUE6.4 pv_globale_multi PVALUE6.4;
format pv PVALUE6.4 pv_multi PVALUE6.4;
merge emr_all emr_all2;
by effect modalite;
run;
/*Summary table*/
options mprint mlogic;
options orientation=portrait;
ods rtf file="&path_out.\&date RR_&nom._&model..rtf"
              style=statistical fontscale=85
              nogfootnote nogtitle;
title; footnote;
title1 "RR for &NOM. model &model.";
proc report data = EMR_fusion
	style(report)=[borderrightcolor=white borderleftcolor=white background=white ]
	style(summary)=[frame=void background=white ]
	style(report column header summary)=[background=white borderrightcolor=white borderleftcolor=white];
column 	effect
		("Modalités" modalite)
		("RR (95% CI)" IC_RER)
      	("p-value" pv)
		("p-value globale" pv_globale)
		("RR Multi.(95% CI)" IC_RER_multi)
      	("p-value Multi." pv_multi)
		("p-value globale Multi." pv_globale_multi)
       ;
run;
footnote;
title;
footnote;
options orientation=portrait;
ods rtf close;
%mend;
/*Macro for automatically selecting variables in outputs*/
%macro sorties_RR(liste,n, obs, model); /*liste: list of variables, n:name, obs:event, model:number of the model (for not erasing outputs)*/
%put &liste;
/*Nb de variables*/
%let nb_elem = %sysevalf(1+%sysfunc(count((&liste),%str( ))));		/*Number of elements in the list (if values are added to/removed from lists, the program still runs)*/
%put &nb_elem;
%let liste_indiv = ;
%let liste_glob = ;
%do i = 1 %to &nb_elem;					/*Loop over the number of elements in the list*/
	%let variab_indiv = RR.rr_&n._%scan(&liste,&i,’ ’);  /*select the ith variable in the list*/
	%let liste_indiv = &liste_indiv &variab_indiv;
	%let variab_glob = RR.rr_glob_&n._%scan(&liste,&i,’ ’);  /*select the ith variable in the list*/
	%let liste_glob = &liste_glob &variab_glob;
%end;
%put &liste_indiv;
%put &liste_glob;
%fusion_RR(nom = &n, model = &model, obs = &obs);
%mend;

/*EXAMPLE for inf_grave (i.e. severe)*/
/*Univariate first*/
/*Liste de toutes les variables*/
%global liste_variables_0;
%let liste_variables_0 = 	c_rt c_rt10 anti_haemophilus Anti_meningocoque_c Pneumovax Prevenar nb_vaccins vaccin coeur_rt sex c_age_diag cl_rt_rate cl_rt_thymus bipulm_bs cl_rt_pituitary cl_rt_foie splenectomy H R T autre pb_k pulm_inf_grave coeur_inf_grave hepat_inf_grave greffe_inf_grave Insuff_immun ttt_era radiotherapie chimiotherapie alkyl mtx anthra c_anthra cl_rt_coeur cl_rt_rate10 cl_rt_rate_m;
%global liste_variables_ref_0;
%let liste_variables_ref_0 = 	coeur_rt(ref="Pas rt coeur") sex(ref="1") c_age_diag(ref="0-5 years") cl_rt_rate(ref="No RT or V5Gy=0") cl_rt_thymus(ref="No RT or D60 thymus = 0, age >3 years") 
							bipulm_bs(ref="No RT or D70% lungs = 0") cl_rt_pituitary(ref="No RT or D80 pituitary=0") cl_rt_foie(ref="No RT or mean dose=0") splenectomy(ref="0") 
							H(ref="0") R(ref="0") T(ref="0") autre(ref="0") pb_k(ref="0") pulm_inf_grave(ref="0") coeur_inf_grave(ref="0") hepat_inf_grave(ref="0") greffe_inf_grave(ref="0") Insuff_immun(ref="0") 
							ttt_era(ref="A before 1980") radiotherapie(ref="0") chimiotherapie(ref="0") alkyl(ref="0") mtx(ref="0") anthra(ref="0") 
							c_anthra(ref="No Anthracyclines") cl_rt_coeur(ref="No radiation dose to heart") cl_rt_rate10(ref="No RT or V10Gy=0") cl_rt_rate_m(ref="No radiation dose to spleen") 
							anti_haemophilus(ref="0") Anti_meningocoque_c(ref="0") Pneumovax(ref="0") Prevenar(ref="0") nb_vaccins(ref="0") vaccin(ref="0") c_rt(ref="Other") c_rt10(ref="Other");
							
/*Few data management to speed up execution and variable class creation*/
data test_ig;
set t_individ1y_inf_grave;
if chimiotherapie = . 	then chimiotherapie = 0;							
if radiotherapie = . 	then radiotherapie = 0;
if chimiotherapie = 0 and radiotherapie = 0 then trt = "A Surgery only";
if chimiotherapie = 1 and radiotherapie = 1 then trt = "D CT and RT";
if chimiotherapie = 1 and radiotherapie = 0 then trt = "C CT no RT";
if chimiotherapie = 0 and radiotherapie = 1 then trt = "B RT no CT";
/*Rectifications*/
if typeg = . then typeg = "0";
if type_iccc = "" then type_iccc = "12.Other and unspecified malignant neoplasms";
format c_age_diag $20.;
if age_diag lt 5 						then c_age_diag = "0-5 years";
if age_diag ge 5 and age_diag lt 10 	then c_age_diag = "5-10 years";
if age_diag ge 10 and age_diag lt 15 	then c_age_diag = "10-15 years";
if age_diag ge 15 						then c_age_diag = "15 years and more";
annee_diag = year(date_diag);
if annee_diag ge 2000 							then ttt_era = "E 1990 and more";/*grouping classes*/
if annee_diag ge 1970 and annee_diag lt 1980 	then ttt_era = "A before 1980";
if annee_diag ge 1980 and annee_diag lt 1990 	then ttt_era = "C 1980-1989";
if annee_diag ge 1990 and annee_diag lt 2000 	then ttt_era = "E 1990 and more";
if annee_diag lt 1970 							then ttt_era = "A before 1980";
if chimiotherapie = . then chimiotherapie = 0;
if radiotherapie = . then radiotherapie = 0;
format trt $22.;
if chimiotherapie = 0 and radiotherapie = 0 then trt = "A Surgery only";
if chimiotherapie = 1 and radiotherapie = 1 then trt = "D CT and RT";
if chimiotherapie = 1 and radiotherapie = 0 then trt = "C CT no RT";
if chimiotherapie = 0 and radiotherapie = 1 then trt = "B RT no CT";
format ct $4. rt $4.;
if chimiotherapie = 0 	then ct = "No";
if chimiotherapie = 1 	then ct = "Yes";
if radiotherapie = 0 	then rt = "No";
if radiotherapie = 1 	then rt = "Yes";
format typeg_clair $22.;
if typeg = 0 then typeg_clair = "Others";
if typeg = 1 then typeg_clair = "Nephroblastoma";
if typeg = 2 then typeg_clair = "Neuroblastoma";
if typeg = 3 then typeg_clair = "Lymphoma";
if typeg = 4 then typeg_clair = "Soft tissue sarcoma";
if typeg = 5 then typeg_clair = "Osteosarcoma";
if typeg = 6 then typeg_clair = "CNS tumor";
if typeg = 7 then typeg_clair = "Gonadal tumor";
if typeg = 8 then typeg_clair = "Thyroid tumor";
if typeg = 9 then typeg_clair = "Retinoblastoma";
format cl_rt_thymus $40.;
cl_rt_thymus="No RT or D60 thymus = 0, age >3 years"; 
if (radiotherapie=0 or D60_thymus le 0) and age_diag le 3 	then cl_rt_thymus="No RT or D60 thymus = 0, age <=3 years"; 
if D60_thymus lt 20 and D60_thymus gt 0 and age_diag le 3 	then cl_rt_thymus="0< Dose <20Gy, age <=3 years";
if D60_thymus lt 20 and D60_thymus gt 0 and age_diag gt 3 	then cl_rt_thymus="0< Dose <20Gy, age >3 years";
if D60_thymus ge 20 and age_diag le 3 						then cl_rt_thymus=">20Gy, <=3 years";
if D60_thymus ge 20 and age_diag gt 3 						then cl_rt_thymus=">20Gy, >3 years";
format bipulm_bs $40.;
bipulm_bs = "No RT or D70% lungs = 0";
if (radiotherapie=0 or D70_poumon le 0) then bipulm_bs="No RT or D70% lungs = 0";
if radiotherapie = 1 and D70_poumon lt 10 and D70_poumon gt 0 then bipulm_bs="RT but no bipulmonary field";
if D70_poumon ge 10 then bipulm_bs="RT with bipulmonary fields";
format cl_rt_pituitary $40.;
cl_rt_pituitary = "No RT or D80 pituitary=0";
if D80_hypophyse ge 35 							then cl_rt_pituitary="More than 35Gy";
if D80_hypophyse gt 0 and D80_hypophyse le 20 	then cl_rt_pituitary="No RT or D80 pituitary=0";
if D80_hypophyse gt 20 and D80_hypophyse le 35 	then cl_rt_pituitary=">20-35Gy";
if radiotherapie=0 or D80_hypophyse le 0 		then cl_rt_pituitary="No RT or D80 pituitary=0";
format cl_rt_foie $40.;
cl_rt_foie = "No RT or mean dose=0";
if radiotherapie = 0 or mean_foie le 0 		then cl_rt_foie = "No RT or mean dose=0";
if mean_foie gt 0 and mean_foie le 5 		then cl_rt_foie = "0-5 Gy";
if mean_foie gt 5 and mean_foie le 20 		then cl_rt_foie = ">5-20 Gy";
if mean_foie gt 20 and mean_foie le 35 		then cl_rt_foie = ">20 Gy";
if mean_foie gt 35 							then cl_rt_foie = ">20 Gy";
IG = inf_grave;
if T = 1 then autre = 1;
if mtx = . 		then mtx = 0;
if anthra = . 	then anthra = 0;
if alkyl = . 	then alkyl = 0;
format c_anthra $30.;
if do_anthra lt 100 						then c_anthra = "<100 mg/m²";
if do_anthra ge 100 and do_anthra le 200 	then c_anthra = "100-200 mg/m²";
if do_anthra gt 200 						then c_anthra = ">200 mg/m²";
if anthra le 0 								then c_anthra = "No Anthracyclines";
format cl_rt_coeur $30.;
if mean_coeur le 0 							then cl_rt_coeur = "No radiation dose to heart";
if mean_coeur gt 0 and mean_coeur le 5 		then cl_rt_coeur = ">0-5 Gy";
if mean_coeur gt 5 and mean_coeur le 10 	then cl_rt_coeur = ">5-10 Gy";
if mean_coeur gt 10 and mean_coeur le 15 	then cl_rt_coeur = ">10-15 Gy";
if mean_coeur gt 15 and mean_coeur le 20 	then cl_rt_coeur = ">15-20 Gy";
if mean_coeur gt 20 						then cl_rt_coeur = ">20 Gy";
if splenectomy = 1 then V20Gy_rate = 0;
if splenectomy = 1 then V15Gy_rate = 0;
if splenectomy = 1 then V5Gy_rate = 0;
if splenectomy = 1 then mean_rate = 0;
format cl_rt_rate_m $30.;
if mean_rate le 0 						then cl_rt_rate_m = "No radiation dose to spleen";
if mean_rate gt 0 and mean_rate le 5 	then cl_rt_rate_m = ">0-5 Gy";
if mean_rate gt 5 and mean_rate le 10 	then cl_rt_rate_m = ">5-10 Gy";
if mean_rate gt 10 and mean_rate le 20 	then cl_rt_rate_m = ">10-20 Gy";
if mean_rate gt 20 						then cl_rt_rate_m = ">20 Gy";
format cl_rt_rate $40.;
cl_rt_rate="No RT or V5Gy=0";
if V5Gy_rate ge 0.8 						then cl_rt_rate="More than 80%";
if V5Gy_rate gt 0 and V5Gy_rate le 0.05	 	then cl_rt_rate=">0-5%";
if V5Gy_rate gt 0.05 and V5Gy_rate le 0.5 	then cl_rt_rate=">5-50%";
if V5Gy_rate gt 0.5 and V5Gy_rate le 0.8 	then cl_rt_rate=">50-80%";
if radiotherapie=0 or V5Gy_rate = 0 		then cl_rt_rate="No RT or V5Gy=0"; 
format cl_rt_rate10 $40.;
cl_rt_rate10="No RT or V10Gy=0";
if V10Gy_rate ge 0.8 						then cl_rt_rate10="More than 80%";
if V10Gy_rate gt 0 and V10Gy_rate le 0.05 	then cl_rt_rate10=">0-5%";
if V10Gy_rate gt 0.05 and V10Gy_rate le 0.5 then cl_rt_rate10=">5-50%";
if V10Gy_rate gt 0.5 and V10Gy_rate le 0.8 	then cl_rt_rate10=">50-80%";
if radiotherapie=0 or V10Gy_rate = 0 		then cl_rt_rate10="No RT or V10Gy=0"; 
format coeur_rt $20.;
if mean_coeur le 15 then coeur_rt = "Pas rt coeur";
if mean_coeur gt 15 then coeur_rt = "RT coeur";
format c_rt $40.;
if radiotherapie = 0 then c_rt = "Other";
if radiotherapie = 1 then c_rt = "Other";
if t = 1 then c_rt = "Thymus";
if r = 1 then c_rt = "Spleen";
if h = 1 then c_rt = "Pituitary";
if V10Gy_rate>=0.8 then rate=2;  
if radiotherapie=0 then rate=0; 
if V10Gy_rate=0 then rate=0;
R10=0; if rate=2 then R10=1;
format c_rt10 $40.;
if radiotherapie = 0 then c_rt10 = "Other";
if radiotherapie = 1 then c_rt10 = "Other";
/*Priorize categories when ties (here, for inf_grave)*/
if t = 1 then c_rt10 = "Thymus";
if r10 = 1 then c_rt10 = "Spleen";
if h = 1 then c_rt10 = "Pituitary";
vaccin = 0;
if nb_vaccins gt 0 then vaccin = 1;
/*shorten the name to avoid lentgh issues later*/
ig = inf_grave;
/*Keep only usefull variables*/
keep id py ig &liste_variables_0 ;
run;
/*Univariate*/
%RR(test_ig,sex,r="1", obs = ig, nom = ig);  
%RR(test_ig,c_age_diag,r="0-5 years", obs = ig, nom = ig);  
%RR(test_ig,cl_rt_rate,r="No RT or V5Gy=0", obs = ig, nom = ig);  
%RR(test_ig,cl_rt_thymus,r="No RT or D60 thymus = 0, age >3 years", obs = ig, nom = ig);  
%RR(test_ig,bipulm_bs,r="No RT or D70% lungs = 0", obs = ig, nom = ig);  
%RR(test_ig,cl_rt_pituitary,r="No RT or D80 pituitary=0", obs = ig, nom = ig);  
%RR(test_ig,cl_rt_foie,r="No RT or mean dose=0", obs = ig, nom = ig);  
%RR(test_ig,splenectomy,r="0", obs = ig, nom = ig);  
%RR(test_ig,H,r="0", obs = ig, nom = ig);  
%RR(test_ig,R,r="0", obs = ig, nom = ig);  
%RR(test_ig,T,r="0", obs = ig, nom = ig);  
%RR(test_ig,autre,r="0", obs = ig, nom = ig);  
%RR(test_ig,pb_k,r="0", obs = ig, nom = ig);  
%RR(test_ig,pulm_inf_grave,r="0", obs = ig, nom = ig);  
%RR(test_ig,coeur_inf_grave,r="0", obs = ig, nom = ig);  
%RR(test_ig,hepat_inf_grave,r="0", obs = ig, nom = ig);  
%RR(test_ig,greffe_inf_grave,r="0", obs = ig, nom = ig);  
%RR(test_ig,Insuff_immun,r="0", obs = ig, nom = ig);  
%RR(test_ig,ttt_era,r="A before 1980", obs = ig, nom = ig);
%RR(test_ig,radiotherapie,r="0", obs = ig, nom = ig);  
%RR(test_ig,chimiotherapie,r="0", obs = ig, nom = ig);  
%RR(test_ig,alkyl,r="0", obs = ig, nom = ig);  
%RR(test_ig,mtx,r="0", obs = ig, nom = ig);  
%RR(test_ig,anthra,r="0", obs = ig, nom = ig);  
%RR(test_ig,c_anthra,r="No Anthracyclines", obs = ig, nom = ig);  
%RR(test_ig,cl_rt_coeur,r="No radiation dose to heart", obs = ig, nom = ig);  
%RR(test_ig,cl_rt_rate10,r="No RT or V10Gy=0", obs = ig, nom = ig);  
%RR(test_ig,cl_rt_rate_m,r="No radiation dose to spleen", obs = ig, nom = ig);  
%RR(test_ig,vaccin,r="0", obs = ig, nom = ig);  
%RR(test_ig,nb_vaccins,r="0", obs = ig, nom = ig);  
%RR(test_ig,anti_haemophilus,r="0", obs = ig, nom = ig);  
%RR(test_ig,Anti_meningocoque_c,r="0", obs = ig, nom = ig);  
%RR(test_ig,Pneumovax,r="0", obs = ig, nom = ig);  
%RR(test_ig,Prevenar,r="0", obs = ig, nom = ig);  
%RR(test_ig,c_rt,r="Other", obs = ig, nom = ig);  
%RR(test_ig,c_rt10,r="Other", obs = ig, nom = ig);  

/*Multivariate*/
/*List variables model 0bis (called 9)*/
%global liste_variables_9;
%let liste_variables_9 = splenectomy radiotherapie ttt_era;
%global liste_variables_ref_9;
%let liste_variables_ref_9 = splenectomy(ref="0") radiotherapie(ref="0") ttt_era(ref="A before 1980");
/*List variables model 1*/
%global liste_variables_1;
%let liste_variables_1 = c_anthra splenectomy ttt_era c_rt;
%global liste_variables_ref_1;
%let liste_variables_ref_1 = c_rt(ref="Other") splenectomy(ref="0") ttt_era(ref="A before 1980") c_anthra(ref="No Anthracyclines");
/*Liste variables model 2*/
%global liste_variables_2;
%let liste_variables_2 = pb_k pulm_inf_grave coeur_inf_grave hepat_inf_grave greffe_inf_grave Insuff_immun ;
%global liste_variables_ref_2;
%let liste_variables_ref_2 = pb_k(ref="0") pulm_inf_grave(ref="0") coeur_inf_grave(ref="0") hepat_inf_grave(ref="0") greffe_inf_grave(ref="0") Insuff_immun(ref="0");
/*Liste variables model 3*/
%global liste_variables_3;
%let liste_variables_3 = c_rt ttt_era Insuff_immun hepat_inf_grave pulm_inf_grave coeur_inf_grave pb_k splenectomy;
%global liste_variables_ref_3;
%let liste_variables_ref_3 = splenectomy(ref="0") pb_k(ref="0") c_rt(ref="Other") ttt_era(ref="A before 1980") Insuff_immun(ref="0") hepat_inf_grave(ref="0") pulm_inf_grave(ref="0") coeur_inf_grave(ref="0");

/*M0: All variables of the list (liste_variables_0) in the same model (to retrieve univariate results only. Delete consecutive multivariate result columns)*/
																											/**********************************************/
%RR_multi(test_ig,ig, ig,0);
%sorties_RR(liste = &liste_variables_0, n = ig, obs = ig, model = 0);
/*M0bis*/
%RR_multi(test_ig,ig, ig,9);
%sorties_RR(liste = &liste_variables_9, n = ig, obs = ig, model = 9);
/*M1*/
%RR_multi(test_ig,ig, ig,1);
%sorties_RR(liste = &liste_variables_1, n = ig, obs = ig, model = 1);
/*M2*/
%RR_multi(test_ig,ig, ig,2);
%sorties_RR(liste = &liste_variables_2, n = ig, obs = ig, model = 2);
/*M3*/
%RR_multi(test_ig,ig, ig,3);
%sorties_RR(liste = &liste_variables_3, n = ig, obs = ig, model = 3);
/*NOTE: In the ouput (rtf file), first columns: RR, pvalue, pvalue globale correspond to univariate analysis. 
Last 3 columns RR_multi, pvalue_multi, pvalue globale_multi correspond to multivariate analysis*/


/*8.Cumulative incidence (CIF)*/
/*EXAMPLE for inf_grave endpoint*/
/*All variables list*/
%global liste_variables;
%let liste_variables = 	sex, c_age_diag, cl_rt_rate, cl_rt_thymus, bipulm_bs, cl_rt_pituitary, cl_rt_foie, splenectomy, H, R, T, autre, pb_k, pulm_choc_sris, coeur_choc_sris, hepat_choc_sris, 
						greffe_choc_sris, Insuff_immun, ttt_era, typeg;
/*8.1.Data management*/
proc sql;
create table cif as
select distinct id, min(entree) as entry format ddmmyy10., max(sortie) as exit format ddmmyy10., max(att_age) as age_exit, max(choc_sris) as obs, 
max(date_diag) as date_diag format ddmmyy10., min(age_diag) as age_diag, max(deces) as deces, ctr, fup, &liste_variables, date_deces, dte_choc_sris, date_nais
from t_individ1y_inf_grave
group by id;
quit;
data cif;
set cif;
/*deces*/
if deces ne 1 then deces = 0;
/*competing risk: death*/
obs_dc = 0;
if deces = 1 then obs_dc = 9;
if obs = 1 then obs_dc = 1;
/*rounding ages (for accelerate calculation)*/
age_entry = (entry-date_nais)/365.25;
if age_entry = age_exit then age_exit = age_exit + 0.5;
fup = round(fup,1);
age_exit2 = round(age_exit,1);
run;
/*8.2.CIF by age*/
proc phreg data = cif plots(overlay=stratum)=cif atrisk;
model age_exit*obs_dc(0)= / rl entry=age_entry eventcode=1 ;
baseline out=courb survival=_all_ loglogs=lls cif=_all_/ nomean;
/*ods output RiskSetInfo=Atrisk;/*-->do not run with competing risks*/
run;
/*calculate the cumulative incidence and the IC95%*/
data courb;
set courb;
format stratus $15.;
cif = cif*100;
uppercif = uppercif*100;
lowercif = lowercif*100;
run;
/*8.3.Select incidence for each age and add expected incidence*/
proc sql;
create table temp as
select a.id, a.nb_hospi_Mean as nb_exp, b.age_exit2
from shr_aer_inf_grave as a left join cif as b
on a.id = b.id;
quit;
proc sql;
create table exp as 
select distinct age_exit2, sum(nb_exp) as expected
from temp
group by age_exit2;
select count(distinct id) into:nb_id
from temp;
select count(*) into:nb_lignes
from exp;
select max(age_exit2) into:max_age
from exp;
quit;
%let max_age = %sysevalf(&max_age);
%macro calc_exp();
%do i = 1 %to &max_age;
	proc sql;
	create table atrisk_&i as 
	select count(distinct id) as nombre_&i
	from temp
	where age_exit2 ge &i;
	quit;
	data atrisk_&i;
	set atrisk_&i;
	age_exit2 = &i;
	run;
%end;
data exp;
merge exp atrisk_1-atrisk_&max_age;
by age_exit2;
run;
data exp;
set exp;
%do i = 1 %to &max_age;
	if age_exit2 = &i then do;
		incidence = expected/nombre_&i*100;
	end;
%end;
if incidence = . then incidence = 0;
format stratus $20.;
stratus = "Expected";
drop nombre_1-nombre_&max_age;
run;
data exp;
set exp;
by stratus;
if first.stratus then cif = incidence;
else cif + incidence;
keep age_exit2 cif stratus;
run;
%supp(atrisk_1-atrisk_&max_age);
%mend;
%calc_exp();
/*reformatting*/
data obs;
set courb;
drop stratus;
run;
data obs;
set obs;
format stratus $30.;
stratus = "Severe infection (95%CI)";
age_exit2 = age_exit;
keep age_exit2 cif uppercif lowercif stratus;
run;
data cif_exp;
set obs exp;
run;
data cif_exp;
set cif_exp;
if stratus = "Severe infection (95%CI)" then cif_step = cif;
if stratus = "Expected" then cif_line = cif;
run;
/*8.4.Drawing graph*/
ods graphics on / reset =all height = 400px width =600px border=off imagefmt=png;
ods rtf file="&path_out.\&date._CumInc_inf_grave_age.rtf";
title ;
proc sgplot data=cif_exp noautolegend;
band x=age_exit2 upper=uppercif lower=lowercif  / legendlabel="95%CI"  transparency=.5  name="band1" group=stratus;
step x=age_exit2 y=cif_step /legendlabel="Incidence" name="step1" group=stratus;
series x=age_exit2 y=cif_line /legendlabel="Incidence" name="step1" group=stratus;
XAXIS
                          LABEL = "Attained age (years)"
                          VALUES = (15 to 70 by 5)
                          Max=50
                                     Min=5
                          OFFSETMIN= 0.05
                          OFFSETMAX= 0.05
                          valueattrs=(weight=Bold style=normal color=Black size=11pt)
                          LABELATTRS=(color=Black weight=bold size= 12pt);
YAXIS
                          LABEL  = "Cumulative severe infection incidence (%)"
                          VALUES = (0 TO 24 BY 2)
                          OFFSETMIN= 0.05
                          OFFSETMAX= 0.1
                                    Min=0
                                    Max=30
                          valueattrs=(weight=Bold style=normal color=Black size=11pt)
                          LABELATTRS=(color=Black weight=bold size= 12pt);
                                DISCRETELEGEND 
                          "Cumulative severe infection incidence (%)"
                                     /noborder
                           across=1
                           location = outside
                           position=Bottom
                           valueattrs=(color=Black size=11pt style=italic weight=bold);
						   keylegend "band1" /location = inside position = NW down = 2;;
run;
quit;
ods graphics off;
ods rtf close;

/*9.Additionnal code for sensitivity*/
/*9.1.COX ph regression*/
/*EXAMPLE for toute_bact endpoint: status = 1 if inf_grave = 1*/
proc phreg data = inf.test_analyse_toute_bact;
class splenectomy(ref="0") sex(ref="1") c_rt(ref="Other") ttt_era(ref="A before 1980") greffe_tb(ref="0") Insuff_immun(ref="0") hepat_tb(ref="0") pulm_tb(ref="0") coeur_tb(ref="0") pb_k(ref="0")
		/param=ref;
model age_sortie*status(0)= 	sex c_rt ttt_era Insuff_immun hepat_tb pulm_tb coeur_tb greffe_tb pb_k splenectomy
								/rl ties=breslow entry = age_entree;
ods output  ModelANOVA=multi.cox_inf_grave ParameterEstimates=multi.cox2_inf_grave;
title1 "Cox ph model";
run;

/*9.2.Andersen Gill method for recurrent event analyse (need previous data management, not shown. See tutorial in appendix of:
Amorim, L. D., & Cai, J. (2015). Modelling recurrent events: a tutorial for analysis in epidemiology. 
International journal of epidemiology, 44(1), 324–333. https://doi-org.proxy.insermbiblio.inist.fr/10.1093/ije/dyu222)*/
/*EXAMPLE for toute_bact endpoint: status = 1 if inf_grave = 1*/
proc phreg data = inf.t_recurrence_toute_bact covs covm;
class splenectomy(ref="0") sex(ref="1") c_rt(ref="Other") ttt_era(ref="A before 1980") greffe_tb(ref="0") Insuff_immun(ref="0") hepat_tb(ref="0") pulm_tb(ref="0") coeur_tb(ref="0") pb_k(ref="0")
		/param=ref;
model age_sortie*status(0)= 	sex c_rt ttt_era Insuff_immun hepat_tb pulm_tb coeur_tb greffe_tb pb_k splenectomy /*if we want to estimate the effect of previous events, add variable "event"*/
								/rl entry = age_entree;
ods output  ModelANOVA=multi.AG_inf_grave ParameterEstimates=multi.AG2_inf_grave;
title1 "AG model for recurrent events";
run;
/*9.3.Marginal means and rate model for recurrent event analyse (same comment than 9.2.)*/
/*EXAMPLE for toute_bact endpoint: status = 1 if inf_grave = 1*/
proc phreg data = inf.t_recurrence_toute_bact covs(aggregate);
class splenectomy(ref="0") sex(ref="1") c_rt(ref="Other") ttt_era(ref="A before 1980") greffe_tb(ref="0") Insuff_immun(ref="0") hepat_tb(ref="0") pulm_tb(ref="0") coeur_tb(ref="0") pb_k(ref="0")
		
		/param=ref;
model age_sortie*status(0)= 	sex c_rt ttt_era Insuff_immun hepat_tb pulm_tb coeur_tb greffe_tb pb_k splenectomy /*if we want to estimate the effect of previous events, add variable "event"*/
								/rl entry = age_entree;
id id;
ods output  ModelANOVA=multi.MMR_inf_grave ParameterEstimates=multi.MMR2_inf_grave;
title1 "Marginal means model for recurrent events";
run;

/*Processing time*/
%let endTime=%QSYSFUNC(DATETIME());
%put endTime %QSYSFUNC(DATETIME(), DATETIME20.3);
%let timeDiff=%sysevalf(&endTime-&startTime);
%put 'The Compile time for this program is approximately ' &timeDiff. 'seconds';  
