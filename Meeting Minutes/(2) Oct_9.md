# Meeting Minutes - Oct 9
## Research Question/Dataset Discussion
* Dataset chosen: Gut & Opioid Use Dataset
  * Note: Opioid Use is only indicated as 'Yes' in this dataset for patients who meet the  Opioid Use Disorder critera according to DSM-4.
* Tentative routes for analysis
  * Carry out basic diversity analysis to start, use the Exploratory Data Analysis route to get best ~10 variables
    * dataset might not be large enough, but still possible to get decent results
  * differential abundance analysis (need to complete additional module 19)
  * complete a predictive analysis on the gut dataset, essentially copying what was done for the skin microbiome paper, but for this  dataset), using random forest analysis (module 20) to check if any variables are correlating
  * potentially use pi-crust to look at metabolic pathways to answer a novel question 
  * Three plots to look at: one for just metadata where we plot each variable against opiod use to get correlation, one for just taxonomy, and one that's combined to see the ranking
* Tentative Research Question
  * Does combining gut microbiome data with clinical biomarkers improve the prediction accuracy of Opioid Use Disorder compared to using either alone?
## Research Proposal Questions
* What parts of the data pipeline do we need to complete for the proposal?
  * We only need to do quality control processing in qiime for the proposal, no diversity analysis yet
## To Do 
* read through the proposal, and bring a list of questions to the next meeting
* start processing metadata in qiime (after Evelyn or Avril confirm that the dataset is available on the server)

