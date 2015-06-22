*** Compute predicted probabilities for pta 
*** Author: S Bauldry
*** Date: April 4, 2015

*** Store estimates from Mplus in matrix A
use "ConLCA/pta-con-est-3", replace

replace var = "_con" if var == "C#1"

gen id     = 1 if var == "W4AGE"
replace id = 2 if var == "FEMALE"
replace id = 3 if var == "BLACK"
replace id = 4 if var == "HISP"
replace id = 5 if var == "ORACE"
replace id = 6 if var == "BIOPAR"
replace id = 7 if var == "PAREDU"
replace id = 8 if var == "PARINC"
replace id = 9 if var == "LIVENV"
replace id = 10 if var == "GPA"
replace id = 11 if var == "COLASP"
replace id = 12 if var == "OBESE"
replace id = 13 if var == "SMOKE"
replace id = 14 if var == "DRINK"
replace id = 15 if var == "POBESE"
replace id = 16 if var == "PSMOKE"
replace id = 17 if var == "PDRINK"
replace id = 18 if var == "F_OBESE"
replace id = 19 if var == "F_SMOKE"
replace id = 20 if var == "F_DRINK"
replace id = 21 if var == "F_POBESE"
replace id = 22 if var == "F_PSMOKE"
replace id = 23 if var == "F_PDRINK"
replace id = 24 if var == "_con"
replace id = 25 if var == "C#2"
replace id = 26 if var == "C#3"
replace id = 27 if var == "C#4"
replace id = 28 if var == "C#5"

keep if model == 1 & pathways == 4
sort id
keep est1 est2 est3
replace est2 = est2[_n + 1] if mi(est2)
replace est3 = est3[_n + 2] if mi(est3)
drop if _n > 18
mkmat est1 est2 est3, mat(A)


*** program to compute desired predicted probabilities
capture program drop PrdPrb
program PrdPrb
	args var val i
	qui replace `var' = `val'
	mata: MPrdPrb("A")
	post pf ("`var'") (`val') (`i') (ppm[1,1]) (ppm[1,2]) (ppm[1,3]) (ppm[1,4])
end

*** mata function for calculation
mata:
void MPrdPrb(string m1) {
	A = st_matrix(m1)
	X   = st_data(., .)
	xb  = X*A
	den = 1 :+ exp( xb[.,1] ) :+ exp( xb[.,2] ) :+ exp( xb[.,3] )
	pp  = exp(xb[.,1]):/den[.,1] , exp(xb[.,2]):/den[.,1] , exp(xb[.,3]):/den[.,1]
	pp  = pp , 1 :- (pp[.,1] :+ pp[.,2] :+ pp[.,3])
	ppm = mean(pp)
	st_matrix("ppm", ppm)
}
end

*** load MI data and compute predicted probabilities
use "pta-data-mi-3", replace

qui tab race, gen(r)
rename (r1 r2 r3 r4) (white black hisp orace)
gen one = 1

order w4age female black hisp orace biopar paredu parinc livenv gpa ///
      colasp obese smoke drink pobese psmoke pdrink one _mi_m
keep w4age-_mi_m
tempfile d1
save `d1', replace


postfile pf str20 var val mi p1 p2 p3 p4 using "pta con prd prb 3", replace
forval i = 1/20 {
	use `d1', replace
	qui keep if _mi_m == `i'
	qui drop _mi_m

	PrdPrb female 1 `i'
	PrdPrb female 0 `i'

	PrdPrb paredu 4 `i'
	PrdPrb paredu 8 `i'

	PrdPrb smoke  0 `i'
	PrdPrb smoke  1 `i'

	PrdPrb psmoke 0 `i'
	PrdPrb psmoke 1 `i'

	PrdPrb drink  0 `i'
	PrdPrb drink  4 `i'

	PrdPrb pdrink 1 `i'
	PrdPrb pdrink 4 `i'
}
postclose pf

use "pta-con-prd-prb-3", replace
collapse (mean) p1 p2 p3 p4, by(val var)

tempfile g1 g2 g3 g4 g5 g6
preserve
keep if var == "female"
reshape long p, i(val) j(path)
recode path (3 = 1) (1 = 2) (4 = 3) (2 = 4)
graph bar (sum) p, over(val, relabel(1 "M" 2 "F") label(labsize(small)) ) ///
   over(path, relabel(1 "WF" 2 "W" 3 "HE" 4 "L")) ///
   ylab(0(.1).53) ytit("pred prob") tit("sex") scheme(lean2) ///
   saving(`g1', replace)
restore

preserve
keep if var == "paredu"
reshape long p, i(val) j(path)
recode path (3 = 1) (1 = 2) (4 = 3) (2 = 4)
graph bar (sum) p, over(val, relabel(1 "HS" 2 "4Y") label(labsize(small)) ) ///
   over(path, relabel(1 "WF" 2 "W" 3 "HE" 4 "L")) ///
   ylab(0(.1).53) ytit("pred prob") tit("parent education") scheme(lean2) ///
   saving(`g2', replace)
restore

preserve
keep if var == "smoke"
reshape long p, i(val) j(path)
recode path (3 = 1) (1 = 2) (4 = 3) (2 = 4)
graph bar (sum) p, over(val, relabel(1 "N" 2 "Y") label(labsize(small)) ) ///
   over(path, relabel(1 "WF" 2 "W" 3 "HE" 4 "L")) ///
   ylab(0(.1).53) ytit("pred prob") tit("adolescent smoked") scheme(lean2) ///
   saving(`g3', replace)
restore

preserve
keep if var == "psmoke"
reshape long p, i(val) j(path)
recode path (3 = 1) (1 = 2) (4 = 3) (2 = 4)
graph bar (sum) p, over(val, relabel(1 "N" 2 "Y") label(labsize(small)) ) ///
   over(path, relabel(1 "WF" 2 "W" 3 "HE" 4 "L")) ///
   ylab(0(.1).53) ytit("pred prob") tit("parent smoked") scheme(lean2) ///
   saving(`g4', replace)
restore

preserve
keep if var == "drink"
reshape long p, i(val) j(path)
recode path (3 = 1) (1 = 2) (4 = 3) (2 = 4)
graph bar (sum) p, over(val, relabel(1 "0" 2 "2/m") label(labsize(small)) ) ///
   over(path, relabel(1 "WF" 2 "W" 3 "HE" 4 "L")) ///
   ylab(0(.1).53) ytit("pred prob") tit("adolescent drinking") scheme(lean2) ///
   saving(`g5', replace)
restore

graph combine "`g1'" "`g2'" "`g3'" "`g5'" "`g4'", scheme(lean2) holes(3)
graph export "pta-fig3.pdf", replace

