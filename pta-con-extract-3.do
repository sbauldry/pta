*** Read Mplus output for conditional LCAs
*** Author: S Bauldry
*** Date: February 11, 2015


capture program drop ExtractFit
program ExtractFit
	args OutFileName PostFile model class
	
	clear
	qui insheet using "`OutFileName'"
	
	qui gen flag     = regexm(v1, "Average number of observations")
	qui replace flag = regexm(v1, "Number of Free Parameters") if flag == 0
	qui replace flag = regexm(v1, "Entropy") if flag == 0
	
	qui gen flag2     = regexm(v1, "H0 Value") 
	qui replace flag2 = regexm(v1, "Akaike") if flag2 == 0
	qui replace flag2 = regexm(v1, "Bayesian") if flag2 == 0
	qui replace flag = 1 if flag2[_n-1] == 1

	qui keep if flag

	qui replace v1 = regexr(v1, "Average number of observations", "")
	qui replace v1 = regexr(v1, "Number of Free Parameters", "")
	qui replace v1 = regexr(v1, "Mean", "")
	qui replace v1 = regexr(v1, "Entropy", "")
	qui replace v1 = regexr(v1, "P-Value", "")
	
	qui destring v1, replace
	
	forval i = 1/7 {
		local p`i' = v1[`i']
	}
	
	post `PostFile' (`model') (`class') (`p1') (`p2') (`p3') (`p4') (`p5') (`p6')
end




capture program drop ExtractClassInfo
program ExtractClassInfo
	args OutFileName PostFile model class
	
	clear
	qui insheet using "`OutFileName'"
	
	qui gen flag2 = regexm(v1, "Class Counts and Proportions")
	qui gen flag3 = regexm(v1, "Average Latent Class")
	
	qui gen flag = 0
	
	forval i = 1/`class' {
		local j = `i' + 2
		qui replace flag = 1 if flag2[_n - `j'] == 1 & flag == 0
		qui replace flag = 1 if flag3[_n - `j'] == 1 & flag == 0
	}
	
	qui keep if flag
	
	qui split v1, gen(x) parse("") destring
	
	forval i = 1/`class' {
		local n`i' = x2[`i']
		
		local j = `i' + 1
		local k = `i' + `class'
		local p`i' = x`j'[`k']
	}
	
	if (`class' < 7) {
		local m = `class' + 1
		forval i = `m'/7 {
			local n`i' = .
			local p`i' = .
		}
	}
	
	post `PostFile' (`model') (`class') (`n1') (`n2') (`n3') (`n4') (`n5') ///
	                (`n6') (`n7') (`p1') (`p2') (`p3') (`p4') (`p5') (`p6') ///
	                (`p7')
end



capture program drop ExtractConItemProb
program ExtractConItemProb
	args OutFileName PostFile model class

	clear
	qui insheet using "`OutFileName'"

	qui gen flag1     = regexm(v1, "MODEL RESULTS")
	qui replace flag1 = regexm(v1, "Categorical Latent Variables") if flag1 == 0
	qui gen flag1a    = sum(flag1)
	qui keep if flag1a == 1

	qui split v1, gen(x) parse("")
	qui drop if x1 == "MODEL" | x1 == "Two-Tailed" | x1 == "Estimate" | ///
	            x1 == "Latent" | x1 == "Thresholds"
	qui destring x2, replace
	qui gen flag2 = regexm(x1, "REL[0-9][0-9][$]1")
	qui gen p     = 1/(1 + exp(x2)) if x2 != 15 | x2 != -15
	qui replace p = 1 if x2 == 0
	qui replace p = 0 if x2 == 1
	qui replace p = 1 - (1 - 1/(1 + exp(x2))) - p[_n + 1] if flag2 == 1
	qui keep p

	qui set obs 245
	qui gen class = 0
	forval i = 1/7 {
		local j = `i'*35
		qui replace class = `i' if _n <= `j' & class == 0
	}
	
	bysort class: gen id = _n
	qui reshape wide p, i(id) j(class)
	
	qui gen series     = "edu" if _n <= 7
	qui replace series = "job" if _n <= 14 & mi(series)
	forval i = 15(2)27 {
		local j = `i' + 1
		qui replace series = "coh" if _n <= `i' & mi(series)
		qui replace series = "mar" if _n <= `j' & mi(series)
	}
	qui replace series = "par" if _n <= 35 & mi(series)

	qui sort series id
	qui bysort series: gen age = _n*2 + 16
	
	qui gen model = `model'
	qui gen class = `class'
	
	order model class series age p1-p7
	qui drop id
	
	qui save `PostFile', replace
end




capture program drop ExtractEstimates
program ExtractEstimates
	args OutFileName PostFile model class
	
	clear
	qui insheet using "`OutFileName'"

	qui gen flag1     = regexm(v1, "Categorical Latent Variables")
	qui replace flag1 = regexm(v1, "QUALITY OF NUMERICAL RESULTS") if flag1 == 0
	qui gen flag1a    = sum(flag1)
	qui keep if flag1a == 1

	qui split v1, gen(x) parse("") destring

	qui gen flag2     = regexm(x1, "Categorical")
	qui replace flag2 = regexm(x1, "Intercepts") if flag2 == 0
	drop if flag2
	drop if x2 == "ON"

	rename (x1 x2 x3 x5) (var est se pval)
	keep var est se pval

	if ( `model' == 1 ) {
		qui gen class = 0
		local rclass = `class' - 1
		forval i = 1/`rclass' {
			local j = `i'*17
			qui replace class = `i' if _n <= `j' & class == 0
			local int = _N - (`rclass' - `i')
			qui replace class = `i' in `int'
		}
	}

	if ( `model' == 2 ) {
		qui gen class = 0
		local rclass = `class' - 1
		forval i = 1/`rclass' {
			local j = `i'*23
			qui replace class = `i' if _n <= `j' & class == 0
			local int = _N - (`rclass' - `i')
			qui replace class = `i' in `int'
		}
	}

	qui reshape wide est se pval, i(var) j(class)
	qui destring, replace
	qui gen model    = `model'
	qui gen pathways = `class'
	qui save `PostFile', replace
end





*** Extract fit statistics
postutil clear
postfile PF1 model class N Np LL AIC BIC E using "pta-con-model-fit-3", replace
forval i = 1/7 {
	ExtractFit "Mplus/pta-con-`i'.out" "PF1" 1 `i'
}
postclose PF1



*** Extract class information
postutil clear
postfile PF2 model class n1 n2 n3 n4 n5 n6 n7 p1 p2 p3 p4 p5 p6 p7 ///
  using "pta-con-class-info-3", replace
forval i = 1/7 {
	ExtractClassInfo "Mplus/pta-con-`i'.out" "PF2" 1 `i'
}
postclose PF2



*** Extract conditional item response probabilities
forval i = 2/7 {
	ExtractConItemProb "Mplus/pta-con-`i'.out" "CIRP1`i'" 1 `i'
}

clear
set obs 0
forval i = 2/7 {
	append using CIRP1`i'
	append using CIRP2`i'
	erase CIRP1`i'.dta
	erase CIRP2`i'.dta
}
save "pta-con-cirp-3", replace



*** Extract estimates
forval i = 4/6 {
	ExtractEstimates "Mplus/pta-con-`i'.out" "Est1`i'" 1 `i'
}

clear
set obs 0
forval i = 4/6 {
	append using Est1`i'
	append using Est2`i'
	erase Est1`i'.dta
	erase Est2`i'.dta
}
save "pta-con-est-3", replace
