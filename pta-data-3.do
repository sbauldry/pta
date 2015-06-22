*** Purpose: To prepare data for analysis of pathways to adulthood and health
*** Author: S Bauldry 
*** Date: January 14, 2015

****** Setting data file macros and initializing tempfiles
global w1 "[path to data]"
global w4 "[path to data]"
global s16 "[path to data]"
global s18 "[path to data]"
global s19 "[path to data]"
global wgt "[path to data]"
global fst "[path to data]"

tempfile d1 d2 d3


****** calculating age of first full-time job and four-year degree
use aid imonth iday4 iyear4 h4od1m h4od1y h4lm5 h4ed1 h4ed3* h4ed4* ///
	using "$w4", replace

recode h4ed1 (6 8 = .)
recode h4ed3a h4ed3b h4ed3c (96 97 98 = .)
recode h4ed3d h4ed3e h4ed3f h4ed3g h4ed3h (7 = .)
recode h4ed4* (9997 9998 = .)
recode h4lm5 (97 98 = .)

gen birthDate = mdy(h4od1m, 15, h4od1y)
gen w4age = floor( (mdy(imonth4, iday4, iyear4) - birthDate)/364.25 )

egen highDegree = rowmax(h4ed3*)
gen highDegreeYear = .
foreach x in a b c d e f g h {
	replace highDegreeYear = h4ed4`x' if (h4ed3`x' == highDegree & highDegree > 3)
}
gen highDegreeAge = floor( (mdy(5, 15, highDegreeYear) - birthDate)/364.25 )

lab def hd 1 "no degree" 2 "voc/tech" 3 "associate's" 4 "bachelor's" ///
	       5 "post-bac" 6 "master's" 7 "PhD" 8 "professional"
lab val highDegree hd
lab var highDegree "highest degree completed"

rename h4lm5 fullTimeJobAge

forval i = 18(2)30 {
	gen edu`i' = ( highDegreeAge < `i' + 1 ) if w4age >= `i'
	gen job`i' = ( fullTimeJobAge < `i' + 1 ) if w4age >= `i'
}

keep aid birthDate w4age highDegree highDegreeAge fullTimeJobAge edu* job*
sort aid
save `d1', replace


****** calculating ages cohabited and married
use aid h4tr25 h4tr27m h4tr27y h4tr28m h4tr28y using "$s16", replace

recode h4tr27m h4tr28m (96 98 = .)
recode h4tr27y h4tr28y (9996 9998 = .)

* marriage & cohabitation
keep if h4tr25 == 1 | h4tr25 == 2

merge m:1 aid using `d1', keepusing(birthDate w4age)
drop _merge

recode h4tr27m h4tr28m (13 = 2) (14 = 5) (15 = 8) (16 = 11)

gen relStartAge = floor( (mdy(h4tr27m, 15, h4tr27y) - birthDate)/364.25 )
gen relEndAge   = floor( (mdy(h4tr28m, 15, h4tr28y) - birthDate)/364.25 )
drop h4tr27m h4tr28m h4tr27y h4tr28y

bysort aid: gen pid = _n
reshape wide h4tr25 relStartAge relEndAge, i(aid) j(pid)

forval i = 18(2)30 {
	gen mar`i' = 0
	gen coh`i' = 0
	gen rel`i' = 0
	
	forval j = 1/22 {
		replace mar`i' = 1 if h4tr25`j' == 1 & relStartAge`j' < `i' + 1 ///
		                      & relEndAge`j' > `i' & w4age >= `i'
		replace coh`i' = 1 if h4tr25`j' == 2 & relStartAge`j' < `i' + 1 ///
		                      & relEndAge`j' > `i' & w4age >= `i' 
	}

	replace coh`i' = 0 if mar`i' == 1
	replace rel`i' = 1 if coh`i' == 1
	replace rel`i' = 2 if mar`i' == 1
}

keep aid mar18-rel30
sort aid
save `d2', replace


****** calculating age of first child
use aid ptnr_id prgno lbno h4lb1 h4lb2m h4lb2y using "$s19", replace
drop if lbno > 1
merge 1:1 aid ptnr_id prgno using "$s18", keepusing(h4pg3m h4pg3y)
drop if _merge != 3
drop _merge
merge m:1 aid using `d1', keepusing(birthDate w4age)
drop _merge

recode h4lb2m h4pg3m (96 97 98 = .)
recode h4lb2y h4pg3y (9996 9997 9998 = .)
recode h4lb1 (6 7 8 = .)

foreach x in m y {
	gen `x'ChildDate = .
	replace `x'ChildDate = h4lb2`x' if h4lb1 == 0
	replace `x'ChildDate = h4pg3`x' if h4lb1 == 1
}

gen parentAge = floor( (mdy(mChildDate, 15, yChildDate) - birthDate)/364.25 )
gsort +aid +parentAge
by aid: keep if _n == 1

forval i = 18(2)30 {
	gen par`i' = ( parentAge < `i' + 1 ) if w4age >= `i'
}

keep aid parentAge par*
sort aid
save `d3', replace



****** extracting wave 1 covariates
use aid bio_sex h1gi4 h1gi5c h1gi6a-h1gi6e pa12 pb8 h1rm1 h1rf1 pa55 h1ir11 ///
    h1ir14 h1ee1 h1ed11-h1ed14 h1to3 h1to15 pa61 pa64 pc49a_2 pc49a_3 h1gh59a ///
	h1gh59b h1gh60 imonth iday iyear h1gi1m h1gi1y using "$w1", replace

recode bio_sex h1gi4 h1gi5c h1gi6a-h1gi6e (6 7 8 9 = .)
recode pa12 pb8 h1rm1 h1rf1 (11 12 96 97 98 99 = .)
recode pa55 (996 = .)
recode h1ir11 (6 8 9 = .)
recode h1ir14 (6 97 98 99 = .)
recode h1ee1 (6 8 = .)
recode h1ed11-h1ed14 (5 6 96 97 98 99 = .)
recode h1to3 (6 8 9 = .) (7 = 0)
recode h1to15 (97 = 7) (96 98 = .)
recode pc49a_2 pc49a_3 (6 8 = .)
recode pa61 (96 = .)
recode pa64 (6 = .) (7 = 0)
recode h1gh59a h1gh59b (96 98 99 = .)
recode h1gh60 (996 998 999 = .)
recode h1gi1m h1gi1y (96 = .)

gen w1age = floor( (mdy(imonth, iday, iyear + 1900) - mdy(h1gi1m, 15, h1gi1y + 1900))/364.25 )

gen female = ( bio_sex == 2 ) if !mi(bio_sex)

gen race = 3 if h1gi4 == 1
replace race = 2 if h1gi6b == 1 & mi(race)
replace race = 4 if (h1gi6c == 1 | h1gi6d == 1 | h1gi6e == 1) & mi(race)
replace race = 1 if h1gi6a == 1 & mi(race)

recode pa12 pb8 h1rm1 h1rf1 (10 = 1) 
gen paredu = max(pa12, pb8)
gen cparedu = max(h1rm1, h1rf1)
replace paredu = cparedu if mi(paredu)
lab var paredu "parent education"

gen parinc = log(pa55 + 1)
lab var parinc "parent income (log)"

gen livenv = 10 - (h1ir11 + h1ir14)
lab var livenv "living environment"

rename h1ee1 colasp
lab var colasp "college aspirations"

recode h1ed11-h1ed14 (1 = 4) (2 = 3) (3 = 2) (4 = 1)
egen gpa = rowmean(h1ed11-h1ed14)
lab var gpa "gpa"

rename h1to3 smoke
lab var smoke "ever smoked"

gen drink = 7 - h1to15
lab var drink "how often drink"

gen pobese = (pc49a_2 == 1 | pc49a_3 == 1)
replace pobese = . if mi(pc49a_2) & mi(pc49a_3)
lab var pobese "parent obese"

rename pa61 pdrink
lab var pdrink "parent how often drink"

rename pa64 psmoke
lab var psmoke "parent smoker"

gen height = 12*h1gh59a + h1gh59b
gen bmi = (h1gh60*703)/(height^2)
lab var bmi "BMI"


*** merge family structure
merge 1:1 aid using "$fst", keepusing(famst5)
drop _merge

recode famst5 (2 3 4 5 = 0), gen(biopar)
lab var biopar "two bio-parent household"


*** merge sample weights 
merge 1:1 aid using "$wgt", keepusing(gswgt4_2 psuscid region)
drop _merge

forval i = 1/3 {
	merge 1:1 aid using `d`i''
	drop _merge
}

*** calculate adolescent obesity thresholds based on age
gen obese     = ( bmi > 25.89 ) if female == 0 & w1age == 13
replace obese = ( bmi > 26.34 ) if female == 0 & w1age == 14
replace obese = ( bmi > 26.39 ) if female == 0 & w1age == 15
replace obese = ( bmi > 27.76 ) if female == 0 & w1age == 16
replace obese = ( bmi > 28.45 ) if female == 0 & w1age == 17
replace obese = ( bmi > 30.07 ) if female == 0 & w1age == 18
replace obese = ( bmi > 30.76 ) if female == 0 & w1age == 19
replace obese = ( bmi > 31.26 ) if female == 0 & w1age >= 20

replace obese = ( bmi > 28.15 ) if female == 1 & w1age == 13
replace obese = ( bmi > 28.41 ) if female == 1 & w1age == 14
replace obese = ( bmi > 27.84 ) if female == 1 & w1age == 15
replace obese = ( bmi > 29.51 ) if female == 1 & w1age == 16
replace obese = ( bmi > 30.27 ) if female == 1 & w1age == 17
replace obese = ( bmi > 30.72 ) if female == 1 & w1age == 18
replace obese = ( bmi > 30.75 ) if female == 1 & w1age == 19
replace obese = ( bmi > 31.20 ) if female == 1 & w1age >= 20

replace obese = ( bmi > 30 ) if w1age > 20
replace obese = . if mi(bmi)
lab var obese "obese"

*** select analysis sample
keep if w4age >= 30 & !mi(w4age)
keep if !mi(gswgt4_2)

*** keep analysis variables
order aid gswgt4_2 psuscid region w4age female race biopar paredu parinc ///
      livenv gpa colasp obese smoke drink pobese psmoke pdrink edu* job* ///
      mar* coh* rel* par18-par30
keep aid-par30

save "pta-data-3", replace 


*** saving data for unconditional analysis in Mplus
preserve
destring aid, replace
keep aid gswgt4_2 edu18-job30 rel18-par30
desc
outsheet using "pta-data-3.txt", replace comma nolabel noname
restore

*** running multiple imputation for conditional LCA in Mplus
mi set flong
mi reg impute race paredu-colasp obese-pdrink 
mi imp chain (mlogit) race (logit) obese smoke pobese psmoke ///
             (ologit) drink pdrink paredu livenv colasp ///
             (regress) parinc gpa = w4age biopar edu18-coh30 par18-par30, ///
             add(20) rseed(39371234) augment by(female)
save "pta-data-mi-3", replace 

*** constructing race indicators and interaction terms
qui tab race, gen(r)
rename (r1 r2 r3 r4) (white black hisp othrace)
foreach x of varlist smoke drink obese psmoke pdrink pobese {
	gen fem_`x' = female*`x'
}

*** saving data files for analysis in Mplus
forval i = 1/20 {
	preserve
	keep if _mi_m == `i'
	drop _mi* coh* mar* race white
	desc
	outsheet using "pta-data-mi-3-`i'.txt", replace comma nolabel noname
	restore
} 



