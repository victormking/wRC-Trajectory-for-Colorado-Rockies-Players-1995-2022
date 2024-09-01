---
title: "Colorado Rockies wRC+ Trajectory Analysis (1995-2022)"
output: html_document
date: "2024-08-31"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

# Step 1: Data Preparation

# libraries
library(tidyverse)
library(Lahman)
library(dplyr)

# We start by retrieving player data for the Colorado Rockies from a wider range of seasons (1995-2022)
# Filter the Batting dataset to only include players from the Colorado Rockies (teamID == "COL")
# within the years 1995 to 2022. We also select relevant columns such as playerID, yearID, and various 
# batting statistics (AB, H, BB, HBP, X2B, X3B, HR, SF). Also, filter to keep only players with
# more than 50 at-bats in a season to ensure a sufficient sample size for meaningful analysis.

rockies_data <- Batting %>%
  filter(teamID == "COL" & yearID %in% c(1995:2022)) %>%
  select(playerID, yearID, AB, H, BB, HBP, X2B, X3B, HR, SF) %>%
  
  filter(AB > 50)  # Filter for players with a minimum of 50 at-bats per season


# Retrieve the most played fielding positions for Rockies players from 1995 to 2022
# Filter the Fielding dataset for the same years (1995-2022). For each player (grouped by playerID), 
# calculate the total number of games played (Games) and determine their primary fielding position (POS)
# based on the position they played the most games at (which.max(Games)).
Positions <- Fielding %>%
  filter(yearID %in% c(1995:2022)) %>%
  group_by(playerID) %>%
  summarize(Games = sum(G), POS = POS[which.max(Games)], .groups = "drop")


# Join the fielding positions with the main Rockies batting data
# Use an inner join to merge the 'rockies_data' and 'Positions' dataframes on 'playerID', 
# combining batting statistics with the primary fielding position for each player.
rockies_data_joined <- rockies_data %>%
  inner_join(Positions, by = "playerID")

str(rockies_data_joined)


# Join with the People dataset to add birth information and calculate age
# Use an inner join to merge 'rockies_data_joined' with the 'People' dataframe to add birth information (birthYear, birthMonth, birthDay).
# Calculate the age of each player at the time of each season. If a player is born after July, we consider the next year for their age calculation.

rockies_data_final <- rockies_data_joined %>%
  inner_join(select(People, playerID, birthYear, birthMonth, birthDay), by = "playerID") %>%
  mutate(Age = yearID - ifelse(birthMonth >= 7, birthYear + 1, birthYear))  # Adjust age calculation for July or later births


str(rockies_data_final)

# Here we convert 'playerID' to a factor to create fixed effects for modeling
# This conversion is necessary for regression analysis to account for player-specific effects.
rockies_data_final$playerID <- as.factor(rockies_data_final$playerID)

str(rockies_data)

# This is to ensure 'rockies_data_final' has the necessary columns before calculating SLG, OBP, and wRC+
if (all(c('H', 'X2B', 'X3B', 'HR', 'AB', 'BB', 'HBP', 'SF') %in% colnames(rockies_data_final))) {
  
  # Calculate advanced statistics (SLG, OBP, and wRC+) for each player-season
  # SLG (Slugging Percentage) and OBP (On-Base Percentage) are calculated using standard formulas.
  # wRC+ (Weighted Runs Created Plus) is a simplified formula that combines SLG and OBP to give a
  # normalized measure of offensive performance, scaled to 100 (average).
  rockies_data_final <- rockies_data_final %>%
    mutate(SLG = (H - X2B - X3B - HR + 2 * X2B + 3 * X3B + 4 * HR) / AB,
           OBP = (H + BB + HBP) / (AB + BB + HBP + SF),
           wRC_plus = (OBP * 1.8 + SLG * 1.3) * 100)  # Simplified wRC+ formula
} else {
  stop("Required columns are not present")
}

# Check the structure of 'rockies_data_final' to ensure columns are present
str(rockies_data_final)

```

```{R}

# Assignment Step 2: Create a Model to Find Peak Age

# Center Age around 30 for the regression model
# This step adjusts the Age variable by subtracting 30 from each player's age to center it around 30.
# Centering helps in interpreting the coefficients of the model, particularly for polynomial terms.
rockies_data_final <- rockies_data_final %>%
  mutate(Age_centered = Age - 30)  # Center Age around 30

# Fit a quadratic regression model for wRC+ over Age with Player Fixed Effect
# The model includes a quadratic term for Age (Age_centered^2) to capture the non-linear relationship
# between age and wRC+. Player fixed effects (playerID) are included to control for individual player differences.
# This allows us to determine the typical career trajectory of wRC+ for an average Rockies player, accounting for variability among players.

wRC_plus_model <- lm(wRC_plus ~ Age_centered + I(Age_centered^2) + playerID, data = rockies_data_final)

# Display the summary of the model to inspect coefficients, R-squared value, and p-values
# The summary provides insights into the model's fit and the significance of the predictors.

summary(wRC_plus_model)

# Here I find the peak age for wRC+
# Using the coefficients of the quadratic model, I calculate the age at which wRC+ is maximized.
# The peak age is computed as the vertex of the parabola represented by the regression equation.

coef <- coef(wRC_plus_model)
peak_age <- 30 - coef['Age_centered'] / (2 * coef['I(Age_centered^2)'])
peak_age

# Calculate the maximum wRC+ value at the peak age
# Substitute the peak age back into the regression equation to calculate the maximum wRC+ value.
# This value represents the highest predicted wRC+ for the average player at their peak age.
max_wRC_plus <- coef['(Intercept)'] + coef['Age_centered'] * (peak_age - 30) + coef['I(Age_centered^2)'] * ((peak_age - 30)^2)


```

```{R}

# Assignment Step 3: Create a Trajectory Plot of the Fitted Regression Model

# Generate a range of ages to plot predictions over, with increments of 0.5 years
# This creates a data frame containing a sequence of ages covering the range found in the dataset.
# These "fake" data points will be used to project the predicted wRC+ values across different ages.

age_range <- data.frame(Age = seq(min(rockies_data_final$Age), max(rockies_data_final$Age), by = 0.5))

# Extract model coefficients for Age-centered terms
# Retrieve the estimated coefficients from the fitted model to compute predictions for wRC+.
# The coefficients include the intercept, linear term (Age_centered), and quadratic term (I(Age_centered^2)).

intercept <- coef(wRC_plus_model)["(Intercept)"]
age_centered_coef <- coef(wRC_plus_model)["Age_centered"]
age_centered_sq_coef <- coef(wRC_plus_model)["I(Age_centered^2)"]

# Predict wRC+ using the centered Age terms
# I Calculate the predicted wRC+ values using the regression equation with the age centered around 30.
# This gives us the trajectory of the predicted wRC+ values over different ages for the "average" player.

age_range$wRC_plus_pred <- intercept + age_centered_coef * (age_range$Age - 30) + age_centered_sq_coef * ((age_range$Age - 30)^2)

# Calculate peak age and max wRC+ using centered age terms
# The peak age represents the age at which a player's predicted wRC+ is at its maximum.
# Max wRC+ is the predicted value of wRC+ at the peak age, derived using the regression equation.

peak_age <- 30 - age_centered_coef / (2 * age_centered_sq_coef)
max_wRC_plus <- intercept + age_centered_coef * (peak_age - 30) + age_centered_sq_coef * ((peak_age - 30)^2)

# I plot the trajectory of wRC+ over Age with annotations
# Use ggplot2 to create a line plot showing the predicted wRC+ over age for the range of ages in the dataset.
# Annotations are added to highlight the "Max" wRC+ and "Peak age" where the predicted wRC+ is highest.

ggplot(age_range, aes(x = Age, y = wRC_plus_pred)) +
  geom_line(color = "blue", size = 1) +
  labs(title = "wRC+ Trajectory for Colorado Rockies Players (1995-2022)",
       x = "Age",
       y = "Predicted wRC+") +
  geom_vline(xintercept = peak_age, linetype = "dotted") +
  annotate("text", x = peak_age, y = max_wRC_plus, label = "Peak age", vjust = -0.5, hjust = -0.1, size = 4) +
  geom_hline(yintercept = max_wRC_plus, linetype = "dotted") +
  annotate("text", x = min(age_range$Age) + 1, y = max_wRC_plus, label = "Max", vjust = -0.5, hjust = 0, size = 4) +
  # Apply theme with a grey background
  theme_gray() +  # Use theme_gray for a basic grey background
  theme(
    panel.background = element_rect(fill = "grey90"),  # Make panel background grey
    plot.background = element_rect(fill = "grey90")    # Make plot background grey
  )

# Saving the final dataset with additional calculated values to a CSV file

write.csv(rockies_data_final, "rockies_data_final.csv", row.names = FALSE)
ggsave("wRC_plus_trajectory.pdf", width = 8, height = 6)


```

```{R}

# Peak Age: 26.68
# Max wRC+: 112.65

# The first player in my dataset, which serves as the baseline for the model's fixed effects, is Jason Bates (playerID = "batesja01"). 

```
