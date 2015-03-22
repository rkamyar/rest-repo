---
title: "ReadMe"
author: "Reza Kamyar"
date: "03/22/2015"
output: html_document
---

## Files:

* `run_analysis.R`: The code for creating a tidy dataset for 'Getting and Cleaning Data' course project. The output is `data_average.txt`.

* `CodeBook.md`: Describes the variables in `run_analysis.R`.

* `averages_data.txt`: The output of `run_analysis.R`.


##  Remarks:

To start, run the script `run_analysis.R`. This performs the following five steps to create the tidy data:

* 1) By using the function `rbind()`, all the data with the same number of columns and the same type are mereged.

* 2) The coulmns which posses mean and standard deviation measures are extracted and named according to `features.txt`.

* 3) IDs and activity names (from `activity_labels.txt`) are assigned to activity data.

* 4) A 180 by 30 tidy dataset including average measures and activity types is created.  

