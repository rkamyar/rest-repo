# Please make sure the package 'plyr' is installed before running the script.
library(plyr)

# Step #1: Merging the training and test sets to create one data set
xTrain <- read.table("train/X_train.txt")
yTrain <- read.table("train/y_train.txt")
subject_train <- read.table("train/subject_train.txt")
xTest <- read.table("test/X_test.txt")
yTest <- read.table("test/y_test.txt")
subjectTest <- read.table("test/subject_test.txt")
xData <- rbind(xTrain, xTest)
yData <- rbind(yTrain, yTest)
subjectData <- rbind(subject_train, subjectTest)



# Step #2: Extracting only the measurements on the mean and standard deviation for each measurement
features <- read.table("features.txt")
mean_std <- grep("-(mean|std)\\(\\)", features[, 2])
xData <- xData[, mean_std]
names(xData) <- features[mean_std, 2]



# Step #3: Using descriptive activity names to name the activities in the data set
activities <- read.table("activity_labels.txt")
yData[, 1] <- activities[yData[, 1], 2]
names(yData) <- "activity"



# Step #4: labeling the data set with descriptive variable names
names(subjectData) <- "subject"
data <- cbind(xData, yData, subjectData)



# Step #5: for each activity and each subject
dataAverage <- ddply(data, .(subject, activity), function(x) colMeans(x[, 1:66]))
write.table(dataAverage, "data_average.txt", row.name=FALSE)


