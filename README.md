Colorado Rockies wRC+ Trajectory Analysis (1995-2022)

Overview

This project analyzes the career trajectories of Colorado Rockies baseball players from 1995 to 2022, focusing on the wRC+ (Weighted Runs Created Plus) metric. The goal is to model and visualize the peak performance age for Rockies players and explore how their performance changes with age. The analysis uses data from the Lahman R package, applying regression techniques to understand the non-linear relationship between age and performance.

Table of Contents
Overview
Data
Methodology
Results
Usage
Dependencies
Repository Structure
Contributing
License
Data

The data used in this analysis comes from the Lahman R package, which includes comprehensive historical data for Major League Baseball. We focus specifically on the Batting, Fielding, and People tables, filtering for players from the Colorado Rockies between 1995 and 2022.

Data Preparation

Batting Data: Filtered to include only players with at least 50 at-bats in a season to ensure meaningful statistics.

Fielding Data: Retrieved the most played fielding positions for each player to account for their primary positions.

People Data: Used to calculate the age of players during each season, adjusting for mid-year births.

Methodology

The analysis follows these steps:

Data Filtering and Cleaning: Filtered relevant data from the Lahman package and merged the datasets to create a comprehensive dataset containing player statistics, primary fielding positions, and calculated ages.

Calculate Advanced Metrics: Calculated wRC+ using simplified formulas based on OBP (On-Base Percentage) and SLG (Slugging Percentage).

Regression Modeling: Applied a quadratic regression model with player fixed effects to capture the non-linear relationship between age and wRC+. The age variable was centered around 30 for better interpretability.

Visualization: Generated a trajectory plot for predicted wRC+ values over the player's age, annotated with the peak age and maximum wRC+ for better visualization and understanding.

Results

The analysis found that the average peak age for Rockies players in terms of wRC+ is approximately 26.7 years, with a maximum wRC+ of around 112.65.
The visualization represents the career trajectory of wRC+ with respect to age, with a grey background theme for clarity.
