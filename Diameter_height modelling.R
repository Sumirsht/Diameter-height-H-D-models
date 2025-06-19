# Title: Diameter-Height Modelling in R
# Author: [Sumir Shrestha]
# Date: [2024-10-01]
# Description: This script performs diameter-height modelling using various non-linear models in R.

# Load necessary libraries
library(dplyr) # Data manipulation
library(tidyverse) # Data manipulation and visualization
library(QuantPsyc) # For standardized beta coefficients
library(car)  # For Durbin-Watson test
library(ggthemes) # For pre-built ggplot themes
library(lmtest) # For Breusch-Pagan test (homoscedasticity)
library(olsrr) # for Rasidual and Cook's distance plot
library(minpack.lm) #for minimum non-linear regression

# Set seed for reproducibility
set.seed(123)

# Generate random DBH values between 0 and 120
DBH <- runif(120, min = 0, max = 120)
Height <- 10 + (DBH / max(DBH)) * (25 - 10)

# Create the data frame
data <- data.frame(
  DBH = DBH,
  Height = Height
)

# View the first few rows
head(data)

#showing the value in plot
plot(data$DBH, data$Height, main = "DBH vs Height", xlab = "DBH", ylab = "Height")

#creating normal model using linear model
m1 <- lm(Height ~ DBH, data = data)
summary(m1)

#applying the model1 as per weibull, 1951 & yang et al., 1978
model1 <- nlsLM(
  Height ~ 1.3 + a * (1 - exp(-b * DBH^c)), 
  data = data, 
  start = list(a = 45, b = 0.01, c = 1)
)

summary(model1)
sjPlot::tab_model(model1)

#applying the model2 as per curtis, 1967
model2 <- nlsLM(
  Height ~ 1.3 + a * exp(b * DBH^c), 
  data = data, 
  start = list(a = 0.05653, b = 4.10361, c = 0.09888)
)

summary(model2)
sjPlot::tab_model(model2)

#applying the model3 as per wykoff et al., 1982
model3 <- nlsLM(
  Height ~ 1.3 + exp(a + (b / DBH)), 
  data = data, 
  start = list(a = 3.62474, b = -19.53154)
)
summary(model3)
sjPlot::tab_model(model3)

#applying the model4 as per Arabatzis & Burkhart, 1992
model4 <- nlsLM(
  Height ~ 1.3 + a * DBH^b, 
  data = data, 
  start = list(a = 2.40698, b = 0.58876)
)
summary(model4)
sjPlot::tab_model(model4)

#applying the model5 as per chapman-Richards cited in sharma(2009)
model5 <- nlsLM(
  Height ~ 1.3 + a * (1 - exp(-b * DBH))^c, 
  data = data, 
  start = list(a = 43.365968, b = 0.013972, c= 0.812427)
)
summary(model5)
sjPlot::tab_model(model5)

#applying the model6 as per Naslund, 1936
model6 <- nlsLM(
  Height ~ 1.3 + ( DBH / (a + b * DBH))^3 , 
  data = data, 
  start = list(a = 2.554273, b = 0.289243)
)
summary(model6)
sjPlot::tab_model(model6)

#applying the model7 as per Gompertz (winsor, 1932)
model7 <- nlsLM(
  Height ~ 1.3 + a * exp (-b * exp(-c *DBH)), 
  data = data, 
  start = list(a = 35.77751, b = 1.96648, c= 0.033792)
)
summary(model7)
sjPlot::tab_model(model7)

#applying the model8 as per Logistic model (Zeide, 1993)
model8 <- nlsLM(
  Height ~ 1.3 + a /(1 + b * exp(-c *DBH)), 
  data = data, 
  start = list(a = 34.056388, b = 4.315559, c= 0.050469)
)
summary(model8)
sjPlot::tab_model(model8)

#applying the model9 as per Meyer, 1940
model9 <- nlsLM(
  Height ~ 1.3 + a * (1 - exp(-b * DBH)), 
  data = data, 
  start = list(a = 37.662306, b = 0.022056)
)
summary(model9)
sjPlot::tab_model(model9)

##applying the model10 as per El Mamoun et al., 2013
model10 <- nlsLM(
  Height ~ 1.3 + exp (a + b / (DBH +1 )), 
  data = data, 
  start = list(a = 3.64434, b = -21.05658)
)
summary(model10)
sjPlot::tab_model(model10)

gtsummary::tbl_regression(model1)
sjPlot::tab_model(model1, model2, model3, model4, model5, model6, model7, model8, model9, model10)

# Calculate AIC and other metrics for each model
# Define metric functions
calculate_rmse <- function(observed, predicted) {
  sqrt(mean((observed - predicted)^2))
}

calculate_bias <- function(observed, predicted) {
  mean(predicted - observed)
}

calculate_mae <- function(observed, predicted) {
  mean(abs(predicted - observed))
}

# Observed values
observed <- data$Height

# Initialize an empty data frame to store results
metrics_df <- data.frame(
  Model = character(),
  AIC = numeric(),
  RMSE = numeric(),
  Mean_Bias = numeric(),
  MAE = numeric(),
  stringsAsFactors = FALSE
)

# Loop through models
for (i in 1:10) {
  # Correct model name assignment
  model_name <- paste0("model", i)
  
  # Check if model exists before using get
  if (exists(model_name)) {
    # Extract model object and predicted values
    model <- get(model_name)
    predicted <- predict(model, newdata = data)
    
    # Calculate metrics
    aic <- AIC(model)
    rmse <- calculate_rmse(observed, predicted)
    bias <- calculate_bias(observed, predicted)
    mae <- calculate_mae(observed, predicted)
    
    # Append results to the data frame
    metrics_df <- rbind(metrics_df, data.frame(
      Model = model_name,
      AIC = aic,
      RMSE = rmse,
      Mean_Bias = bias,
      MAE = mae
    ))
  } else {
    warning(paste("Model", model_name, "does not exist"))
  }
}

# Print the results
print(metrics_df)

# Optionally, save to a CSV file
write.csv(metrics_df, "model_metrics.csv", row.names = FALSE)

#copyright disclaimer
# This code is a work of fiction and is not intended to represent any real-world data or models.
# The models and metrics used are for illustrative purposes only and should not be used for actual scientific analysis.


