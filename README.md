# Marginal and Conditional Multiple Inference for Linear Mixed Model Predictors

This repository contains the code which generates the tables and figures in the article 'Marginal and Conditional Multiple Inference for Linear Mixed Model Predictors' and its Supplement. 
While the case study is performed in R version 4.1.2, the simulation study is coded in Julia version 1.6.0. 

---

### Real Data Study on Covid-19 Mortality in US State Prisons

The study is replicated by running the script <tt>study/study.R</tt>. The script calls the data of state governor party affiliation <tt>study/governor-party.csv</tt> and census regions <tt>data/census-regions.csv</tt> as well as the data which the New York Times gathered [here](https://github.com/kramlinger/covid-19-data). 
It computes the estimates of the the marginal and conditional covariance matrices as well as the conditional non-centrality parameter and thereupon performs the tests in Section 4 in the main article.  

In order to replicate the bootstrap study on the reliablity of the non-centrality parameter, the script  <tt>study/bootstrap.R</tt> is called within <tt>study/study.R</tt>. If not of interest, it is advised to skip this step, as the bootstrap simulation takes time. 

The file <tt>study/study.R</tt> also generates a map of US states and produces the t-test, which assumptions are not met. Further details are provided in the main article. 

The plots of the residual analysis which are provided in the Supplement are generated in the file <tt>study/residual_analysis.R</tt>. It relies on the calculations of <tt>study/study.R</tt>, and has therefore be executed afterwards. 


---

### Performance Study

The simulation study in the main article considers a nested error regression model with a balanced panel. 
This allows for several simplifications, which are provided and derived in <tt>docs/simulations.pdf</tt>. 
Adaptations for an unbalanced panel and the Fay-Herriot model are straight-forward, but require more effort. 

All tables in the main article and supplement can be reproduced by calling the corresponding simulated scenario: 
<tt>simulation/balanced.jl</tt> for the balanced panel (Tables 1, 5, 8, 9 and 16), 
<tt>simulation/unbalanced.jl</tt> for the unbalanced panel (Table 6), 
<tt>simulation/fhm.jl</tt> for the Fay-Herriot model (Table 7), 
<tt>simulation/bootstrap.jl</tt> for the bootstrap algorithms (Tables 10-14) and 
<tt>simulation/marginal.jl</tt> for the marginal interpretation (Table 15). 

After adjusting the table parameters to the desired specification within the document, the program prints out the
resulting coverage probability together with power, relative volume and estimated non-centrality parameters in the accordingly names objects. 
