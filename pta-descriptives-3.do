*** Purpose: Descriptive statistics for PtA Health analysis
*** Author: S Bauldry
*** Date: February 4, 2015

global d1 "[path to non-mi data]"
global d2 "[path to mi data]"

use "$d1", replace

*** missing data
* rates of missing data
local cv w4age female race biopar paredu parinc livenv gpa colasp ///
         obese smoke drink pobese psmoke pdrink
foreach x of varlist `cv' {
	gen m`x' = ( mi(`x') )
}
sum m*
sum mparinc mlivenv mpobese mpsmoke mpdrink

* diagnostic plots
use "$d2", replace

capture program drop MissFig
program MissFig
	args v g
	forval i = 1/20 {
		midiagplots `v', m(`i') plottype(histogram) scheme(lean2) ///
   		   legend(pos(6) rows(1) order(1 "obs" 2 "imp" 3 "com")) ///
   		   nodraw saving(g`i', replace)
    }
    grc1leg g1_1_`v'.gph g2_2_`v'.gph g3_3_`v'.gph g4_4_`v'.gph ///
       g5_5_`v'.gph g6_6_`v'.gph g7_7_`v'.gph g8_8_`v'.gph ///
       g9_9_`v'.gph g10_10_`v'.gph g11_11_`v'.gph g12_12_`v'.gph ///
       g13_13_`v'.gph g14_14_`v'.gph g15_15_`v'.gph g16_16_`v'.gph ///
       g17_17_`v'.gph g18_18_`v'.gph g19_19_`v'.gph g20_20_`v'.gph, ///
       leg(g1_1_`v'.gph) scheme(lean2) 
    graph export `g'.pdf, replace
    forval i = 1/20 {
    	erase g`i'_`i'_`v'.gph
    }
end

MissFig parinc AFig1
MissFig livenv AFig2
MissFig pobese AFig3
MissFig psmoke AFig4
MissFig pdrink AFig5

* more detail for income
sum parinc if _mi_m == 0
local m1 = r(mean)
local s1 = r(sd)

forval i = 1/20 {
	qui sum parinc if _mi_m == `i'
	local m2 = `m1' - r(mean)
	local s2 = `s1' - r(sd)
	dis "imputation `i': dif mean = " %5.2f `m2' " dif sd = " %5.2f `s2'
}




*** Figure 1
use "$d1", replace
svyset psuscid [pw = gswgt4_2], strata(region)
tempfile Fig1
postutil clear
postfile PF fem age job edu coh mar par using `Fig1'

forval i = 18(2)30 {
	foreach x in job edu coh mar par {
		qui svy: prop `x'`i'
		mat b = e(b)
		local p`x' = b[1,2]
	}
	post PF (.) (`i') (`pjob') (`pedu') (`pcoh') (`pmar') (`ppar')
		
	forval j = 0/1 {
		foreach x in job edu coh mar par {
			qui svy: prop `x'`i' if female == `j'
			mat b = e(b)
			local p`x' = b[1,2]
		}
		post PF (`j') (`i') (`pjob') (`pedu') (`pcoh') (`pmar') (`ppar')
	}
}

postclose PF

preserve
use `Fig1', replace
replace edu = 0 if fem == 1 & age == 18
foreach x of varlist edu job coh mar par {
	replace `x' = 100*`x'
}

twoway (connected edu age if mi(fem), lp(dash_dot_dot) m(o)) ///
	   (connected job age if mi(fem), lp(dash) m(oh)) ///
	   (connected coh age if mi(fem), lp(dash_dot) m(s)) ///
	   (connected mar age if mi(fem), lp(solid) m(d)) ///
	   (connected par age if mi(fem), lp(shortdash) m(dh)), ///
	   ylab(0(20)100) xlab(18(2)30) scheme(lean2) tit("") ///
	   ytit("% respondents") legend( rows(1) pos(6) )
graph export Fig1.pdf, replace
restore




*** Table 1
use "$d2", replace

qui tab race, gen(r)
rename (r1 r2 r3 r4) (white black hisp orace)

local cv1 obese smoke drink 
local cv2 pobese psmoke pdrink
local cv3 w4age female white black hisp orace biopar paredu parinc livenv gpa colasp

foreach x of varlist `cv1' `cv2' `cv3' {
	qui mi est: mean `x' [pw = gswgt4_2]
	mat m = e(b_mi)
	local om = m[1,1]

	dis "`x' " as res %5.2f `om'
}




