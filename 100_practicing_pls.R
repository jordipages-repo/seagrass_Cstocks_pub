
library(tidyverse)
library(tidylog)
library(pls)
library(caret)


# Load the data
data("Boston", package = "MASS")


# Split the data into training and test set
# set.seed(123)
training.samples <- Boston$medv %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- Boston[training.samples, ]
test.data <- Boston[-training.samples, ]


# Computing principal component regression (PCR) ----
# Build the model on training set
# set.seed(123)
model <- train(
  medv~., data = train.data, method = "pcr",
  scale = TRUE,
  trControl = trainControl("cv", number = 10),
  tuneLength = 10
)
# Plot model RMSE vs different values of components
plot(model)
# Print the best tuning parameter ncomp that
# minimize the cross-validation error, RMSE
model$bestTune

# Summarize the final model
summary(model$finalModel)

# Make predictions
predictions <- model %>% predict(test.data)
# Model performance metrics
data.frame(
  RMSE = caret::RMSE(predictions, test.data$medv),
  Rsquare = caret::R2(predictions, test.data$medv)
)

# CONCLUSIONS
# Taken together, cross-validation identifies ncomp = 5 as the optimal number of 
# PCs that minimize the prediction error (RMSE) and explains enough variation in 
# the predictors and in the outcome. 


# Computing partial least squares ----
# Build the model on training set
set.seed(123)
model <- train(
  medv~., data = train.data, method = "pls",
  scale = TRUE,
  trControl = trainControl("cv", number = 10),
  tuneLength = 10
)
# Plot model RMSE vs different values of components
plot(model)
# Print the best tuning parameter ncomp that
# minimize the cross-validation error, RMSE
model$bestTune

# Summarize the final model
summary(model$finalModel)

# Make predictions
predictions <- model %>% predict(test.data)
# Model performance metrics
data.frame(
  RMSE = caret::RMSE(predictions, test.data$medv),
  Rsquare = caret::R2(predictions, test.data$medv)
)
