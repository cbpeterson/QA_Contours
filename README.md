# QA_Contours

### Author: [Zachary Wooten](https://github.com/wootz101)

This repository contains code and data to reproduce most figures in the paper: ***Predictive modeling using shape statistics for interpretable and robust quality assurance of automated contours in radiation treatment planning***

 Written by Zachary Wooten, Cenji Yu, Laurence Court, and Christine Peterson in 2022.

 - **MAIN_ShapeStatistics.R** This is the main file to run that will reproduce all figures from extracted shape data. 
 
 - **shapeHist.R** This file contains the "shapeHist" function that takes in a 3-dimensional binary matrix representation of an organ contour and returns the histogram shape features

- **bad_kidney_hist_PSB.csv** : The Histogram data of shape features from 52 unaccaptable kidney contours

- **good_kidney_hist_PSB.csv** : The Histogram data of shape features from 260 accaptable kidney contours

- **bad_kidney_2D_PSB.csv** : The shape features for every axial slice per patient from 52 unaccaptable kidney contours

- **good_kidney_2D_PSB.csv** : The shape features for every axial slice per patient from 260 accaptable kidney contours

- **unlab_kidney_hist_PSB.csv** : The Histogram data of shape features from 36 unlabeled kidney contours

