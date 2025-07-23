# notes on the US bvar


# Fred variable mapping

Note the indexed variables are indexed to 2017, not 2012 as in older vintages. 

Transformations are:
log100 = log_e * 100
x100 = * 100

Real GDP: qgdp, GDPC1; log100
Real PCE: qcons, PCECC96; log100
Real Gross Private Domestic Investment: qifixnr, GPDIC1 (note: Real Gross Private Domestic Investment over nonres to get a longer time series); log100
Real Exports of Goods and Services: qc, EXPGSC1; log100
Real Imports of Goods and Services, qm, IMPGSC1; log100
GDP Deflator: pgdp, GDPDEF; log100
PCE chain price index: pc, PCEPI; log100
PCE core: pcstrip, PCEPILFE; log100
CPI all: cpiu, CPIAUCSL; log100
Potential GDP: qgdpfe, ; log100
Non-farm payroll employment: eea, PAYEMS (note: monthly data, needs to be aggregated to quarterly); log100
Civilian labour force: lc, CLF16OV (note: monthly data, needs to be aggregated to quarterly); log100
Wages (GDI: Wages & salaries (billions, SAAR)): wsd, A4102C1Q027SBEA (note: Would have preferred ECI but series too short); log100
Non-farm business sector hours: lxnfh, HOANBS; log100
Unemployment rate: ruc, UNRATE (note: As with other monthly variables, take the quarterly average); x100
Federal funds effective rate: rmfedfundns, FEDFUNDS (note: again, this is monthly so average up to quarterly); x100
3-month T bill: rmgbs3ns, TB3MS (note: monthly, so same treatment as others to make it quarterly); x100
10-year Treasury note: r10yr, GS10 (note: monthly so same thing again); x100
1-year Treasury: GS1 (note: monthly variable so will need aggregating to quarterly); x100
Moody's Seasoned Aaa Corporate Bond Yield: aaa (note: monthly variable so will need aggregating to quarterly); x100
Industrial production: INDPRO (note: Monthly so needs averaging to quarters like the others); log100










Fred API key: 3f66fc645f49a16dba8369e2539076f8



Model takes 

Conditions:
1. Federal funds rate
2. 3-month T-bill
3. 10-year T-note
4. U3 rate
5. PCE price index 
6. PCE core price index

Forecast:
1. Payroll employment
2. People in the labour force 
3. House of work 
4. Wages and salaries 

5. GDP price index 
6. CPI all
7. CPI food at home 
8. CPI medical care

9. 5-year Treasury note
10. Aaa bond 
11. Baa bond

12. Real GDPP
13. Real PCE
14. Real nonres. fixed investment
15. Real exports
16. Real imports
17. TFP
18. Real potential GDP
19. Nominal GNP
20. Nominal private nonres. fixed investment in equipment