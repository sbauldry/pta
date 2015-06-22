*** Analyze Mplus results for conditional LCAs
*** Author: S Bauldry
*** Date: February 11, 2015


*** Model fit statistics.
use "pta-con-model-fit-3", replace
sort model class
gen DBIC = (( BIC[_n-1] - BIC[_n] )/(BIC[_n-1]) )*100
format LL BIC* DBIC* %9.0f E* 9.2f
list model class LL BIC DBIC E, clean


*** Sample size and percentage for different classes.
use "pta-con-class-info-3", replace
sort model class
forval i = 1/7 {
	gen pn`i'     = n`i'/5094
}
format pn1-pn6 %5.2f

list model class n1-n6 if class == 4 | class == 5, clean
list model class pn1-pn6 if class == 4 | class == 5, clean
list model class p1-p6 if class == 4 | class == 5, clean


*** Figures for conditional item response probabilities
use "pta-con-cirp-3", replace

* Program for pathway figures
capture program drop PathFig
program PathFig
	args mod cls pst tit gn
	preserve
	keep if model == `mod' & class == `cls'
	
	graph twoway (connected `pst' age if series == "edu", lp(dash_dot_dot) m(o)) ///
	             (connected `pst' age if series == "job", lp(dash) m(oh))        ///
	             (connected `pst' age if series == "coh", lp(dash_dot) m(s))     ///
				 (connected `pst' age if series == "mar", lp(solid) m(d))        ///
				 (connected `pst' age if series == "par", lp(shortdash) m(dh)),  ///
				 xlab(18(2)30) ylab(0(.2)1.05) ytit("probability")               ///
				 legend( rows(1) order(1 "edu" 2 "job" 3 "coh" 4 "mar" 5 "par")  ///
				         pos(6) ) ///
				 scheme(lean2) title("`tit'") saving(`gn', replace)
end

tempfile g1 g2 g3 g4
PathFig 1 4 p1 "work (34%)" `g1'
PathFig 1 4 p2 "limited (8%)" `g2'
PathFig 1 4 p3 "work & family (35%)" `g3'
PathFig 1 4 p4 "higher ed (24%)" `g4'
grc1leg "`g3'" "`g1'" "`g4'" "`g2'", scheme(lean2) legendfrom("`g1'")
graph export "pta fig2a 2.pdf", replace

tempfile g1 g2 g3 g4 g5
PathFig 1 5 p1 "work & family (22%)" `g1'
PathFig 1 5 p2 "work & children (14%)" `g2'
PathFig 1 5 p3 "higher ed (23%)" `g3'
PathFig 1 5 p4 "work (33%)" `g4'
PathFig 1 5 p5 "limited (8%)" `g5'
grc1leg "`g1'" "`g2'" "`g4'" "`g3'" "`g5'", scheme(lean2) legendfrom("`g1'")
graph export "pta fig2b 2.pdf", replace



*** Parameter estimates for covariates
use "pta-con-est-3", replace

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

* program to compute contrasts 
* approximate se w/out full asymptotic covariance matrix
capture program drop ComCon
program ComCon
	args con p1 p2
	gen e`con' = est`p1' - est`p2'
	gen s`con' = sqrt( se`p1'^2 + se`p2'^2 )
	gen p`con' = 2*(1 - normal( abs(e`con'/s`con') ))
	gen o`con' = exp(e`con')
	gen r`con' = o`con'*s`con'
end


* 4-pathway model 1
* p1 W, p2 L, p3 WF, p4 HE [ref]
preserve
keep if model == 1 & pathways == 4
sort id
replace est4 = 0
replace se4 = 0

ComCon WF_L  3 2
ComCon W_L   1 2
ComCon HE_L  4 2
ComCon WF_HE 3 4
ComCon W_HE  1 4
ComCon WF_W  3 1

foreach x in WF_L W_L HE_L WF_HE W_HE WF_W {
	gen a`x'     = "*" if p`x' < 0.05
	replace a`x' = "**" if p`x' < 0.01
	replace a`x' = "***" if p`x' < 0.001
}

forval i = 1/17 {
	foreach x in WF_L W_L HE_L WF_HE W_HE WF_W {
		local o`x' = o`x'[`i']
		local r`x' = r`x'[`i']
		local a`x' = a`x'[`i']
	}

	dis %5.2f `oWF_L'  " (" %4.2f `rWF_L'  ") `aWF_L', " ///
	    %5.2f `oW_L'   " (" %4.2f `rW_L'   ") `aW_L', " ///
	    %5.2f `oHE_L'  " (" %4.2f `rHE_L'  ") `aHE_L', " ///
	    %5.2f `oWF_HE' " (" %4.2f `rWF_HE' ") `aWF_HE', " ///
	    %5.2f `oW_HE'  " (" %4.2f `rW_HE'  ") `aW_HE', " ///
	    %5.2f `oWF_W'  " (" %4.2f `rWF_W'  ") `aWF_W'" 
}
restore


* 5-pathway model 1
* p1 WF, p2 WC, p3 HE, p4 W, p5 L [ref]
preserve
keep if model == 1 & pathways == 5
sort id
replace est5 = 0
replace se5 = 0

ComCon WF_L  1 5
ComCon WC_L  2 5
ComCon W_L   4 5
ComCon HE_L  3 5
ComCon WF_HE 1 3
ComCon WC_HE 2 3
ComCon W_HE  4 3
ComCon WF_W  1 4
ComCon WC_W  2 4
ComCon WF_WC 1 2

foreach x in WF_L WC_L W_L HE_L WF_HE WC_HE W_HE WF_W WC_W WF_WC {
	gen a`x'     = "*" if p`x' < 0.05
	replace a`x' = "**" if p`x' < 0.01
	replace a`x' = "***" if p`x' < 0.001
}

forval i = 1/17 {
	foreach x in WF_L WC_L W_L HE_L WF_HE WC_HE W_HE WF_W WC_W WF_WC {
		local o`x' = o`x'[`i']
		local r`x' = r`x'[`i']
		local a`x' = a`x'[`i']
	}

	dis %5.2f `oWF_L'  " (" %4.2f `rWF_L'  ") `aWF_L', " ///
	    %5.2f `oWC_L'  " (" %4.2f `rWC_L'  ") `aWC_L', " ///
	    %5.2f `oW_L'   " (" %4.2f `rW_L'   ") `aW_L', " ///
	    %5.2f `oHE_L'  " (" %4.2f `rHE_L'  ") `aHE_L', " ///
	    %5.2f `oWF_HE' " (" %4.2f `rWF_HE' ") `aWF_HE', " ///
	    %5.2f `oWC_HE' " (" %4.2f `rWC_HE' ") `aWC_HE', " ///
	    %5.2f `oW_HE'  " (" %4.2f `rW_HE'  ") `aW_HE', " ///
	    %5.2f `oWF_W'  " (" %4.2f `rWF_W'  ") `aWF_W', " ///
	    %5.2f `oWC_W'  " (" %4.2f `rWC_W'  ") `aWC_W', " ///
	    %5.2f `oWF_WC' " (" %4.2f `rWF_WC' ") `aWF_WC'" 
}
restore
