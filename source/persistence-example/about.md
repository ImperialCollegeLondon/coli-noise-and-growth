
# Example: 'persistence'

~~~~
Description outdated a bit.
~~~~

* as the sge-lnm-model, except that the protein is 'toxic'
* the toxicity of the protein model is modeled by reducing the individual growth rate as a function of the protein

## Simulation code

* here we need to update the elongation rate many times during the cell cycle 
* also absolute time matters to predict a proportion of slow growers / fast growers

## Input data

* one table with the simulation parameters: num generations, num lineages, update_period for 'instataneous' alpha
* one table with the SJ data estimated linear trend for lnm.a,sigma_1,sigma_2 and elongation rate CV with the division rate
* 'demand' tables, i.e. gene expression parameters for different 'nutrient-imposed' division rates

## What is the output

* For each demand, a table that simply lists for each imposed mu, a vector of individual cells 'effective division rate', computed from their division time, corrected for a size doubling (if Ld << 2 * Lb, 1/div_time not meaningful)

## What are the plots

* For now, just histograms of those individual cell 'division rates'
* Those histograms visually demonstrate bistability if bimodal

## To do

* my current update period for alpha is 3 seconds... can probably be relaxed a lot, will speed up simulations

* CHECKED: same results with updating alpha every 1 minute (make sense, still 20 updates in one cycle of fast growing bacteria, ~ 5X speedup)



