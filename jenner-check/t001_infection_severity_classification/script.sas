/* Adapted from "Infection Study - BS.sas", section 2 ("Code list for infections (ICD-9-10)")
   and the classification block around lines 206-215 (data tranF3).
   Original: Boris Schwartz, CESP-INSERM-U1018 Radiation Epidemiology Team - 2024.

   The original script classifies each recorded infection event into one of three
   severity categories (low grade / sepsis / severe shock-SRIS) by matching the first
   3 or 4 characters of an ICD-9/10 diagnosis code against fixed code lists. That
   classification step is self-contained logic and does not depend on the protected
   patient-level cohort database ("inf.t_all_inf") used elsewhere in the program, so
   it is reproduced here driven by a small mock event table instead of the real
   cohort extract. Two syntax adaptations were needed to run this exact logic on
   the current hosted engine (both filed as engine regressions, see the code
   comments below): the code lists are referenced as `in &listname` rather than
   `in (&listname)` (the macro variable already expands to a parenthesized list),
   and `inf` is given an explicit LENGTH statement rather than relying on the
   FORMAT statement alone to set its width. The classification rules themselves
   -- which ICD-9/10 codes count as low-grade, sepsis, or severe/SRIS infection --
   are unchanged from the original. */

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

/* Mock event table standing in for a slice of "inf.t_all_inf" (one row per
   recorded infection event per patient, already transposed to a single code column) */
data events;
	length id 8 code1 $4;
	input id code1 $;
	datalines;
1 A483
1 A021
2 A400
2 J130
3 R572
4 0270
5 M010
6 B959
;
run;

/* Same matching logic as "data tranF3" (lines 206-215 of the original program).
   Note: the original writes `in (&liste_sepsis4)` -- with an extra pair of
   parens around the macro reference, since &liste_sepsis4 already expands to
   a parenthesized list. That double-parenthesized IN list is valid SAS 9.4
   but is mishandled by the current Jenner release, so this bundle references
   the macro list without the redundant outer parens (`in &liste_sepsis4`),
   which is equivalent SAS and keeps the classification logic identical. */
data classified;
	set events;
	length inf $20.;
	format inf $20. inf_num 8.;
	if substr(code1, 1, 4) in &liste_sepsis4 or substr(code1, 1, 3) in &liste_sepsis3 		then inf = "sepsis";
	if substr(code1, 1, 4) in &liste_choc_sris4 or substr(code1, 1, 3) in &liste_choc_sris3 then inf = "choc_sris";
	if substr(code1, 1, 4) in &liste_low4 or substr(code1, 1, 3) in &liste_low3 			then inf = "low";
	if inf = "low" 			then inf_num = 1;
	if inf = "sepsis" 		then inf_num = 2;
	if inf = "choc_sris" 	then inf_num = 3;
run;

proc freq data = classified;
	table inf / missing;
run;

proc print data=classified noobs;
	var id code1 inf inf_num;
run;
