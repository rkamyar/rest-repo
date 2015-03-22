---
title: "CodeBook.md"
author: "Reza Kamyar"
date: "03/22/2015"
output: html_document
---

# List of Variables:

* `xTrain`, `yTrain`, `xTest`, `yTest`, `subjectTrain` and `subjectTest` are the training and testing data downloaded from:

* `xData` is the merged `xTrain` and `xTest`.  `yData` is the merged `yTrain` and `yTest`. `subjectData` is the merged `subject_train` and `subjectTest`.

* `data` is the merged `xData`, `yData` and `subjectData`.

* `features` includes the corresponding names for the `xData`. `mean_and_std_features` containts The colNames.

* A similar approach is taken with activity names through the `activities` variable.

* `averages_data` is the tidy data which includes the average measures and activity types.

