# Determinants of Regional Income Inequality across Indian States: An Econometric Analysis

## Project Overview

This repository contains an econometric analysis of regional income inequality across Indian states.  
The study investigates the determinants of **Net State Domestic Product per capita (NSDP per capita)** using a **cross-sectional regression framework** for the financial year **2022–23**.

The analysis focuses on the role of **industrial output, social sector expenditure, power availability, and agricultural credit** in explaining inter-state disparities in income levels. Emphasis is placed on **rigorous model specification, diagnostic testing, and robustness checks** to ensure the validity of the empirical results.

---

## Data Sources

All data used in this study are sourced from the **Reserve Bank of India (RBI)**:

- **RBI – Handbook of Statistics on Indian States (HSIS), 2022–23**


---

## Variables Used

### Dependent Variable
- **Log of Net State Domestic Product per capita** (`ln_nsdp_pc`)

### Independent Variables
- **Log of Gross State Value Added by Industry** (`ln_gva_industry`)
- **Log of Social Sector Expenditure** (`ln_social_exp`)
- **Log of Per Capita Availability of Power** (`ln_power_pc`)
- **Square Root of Agricultural Credit** (`credit_agri_sqrt`)

All transformations were applied to address violations of **Classical Linear Regression Model (CLRM)** assumptions.

---

## Key Features of the Analysis

### 1. Descriptive Analysis
- Summary statistics for all variables  
- Examination of regional disparities in income and development indicators  

### 2. Model Specification
- Baseline OLS regression  
- Functional form corrections using logarithmic and power transformations  

### 3. Diagnostic Testing
- **Linearity**: Component-plus-residual (CR) plots  
- **Homoscedasticity**: Breusch–Pagan and Goldfeld–Quandt tests  
- **Normality**: Q–Q plots, Shapiro–Wilk and Anderson–Darling tests  
- **Multicollinearity**: Variance Inflation Factors (VIF), eigenvalues, condition indices  

### 4. Model Selection
- Subset selection using **Mallows’ Cp criterion**

### 5. Hypothesis Testing
- ANOVA tests for joint and individual significance of regressors  

### 6. Influence Analysis
- Leverage values (hat values)  
- Studentized residuals  
- Cook’s Distance, DFFITS, and DFBETAS  
- Robustness checks after removing influential observations  

---

## Results and Findings

- **Industrial output** and **power availability** exert a strong positive influence on per capita NSDP.
- **Agricultural credit** positively contributes to regional income levels.
- **Social sector expenditure** exhibits a negative and statistically significant coefficient, suggesting short-run inverse associations that may reflect redistributive targeting, inefficiencies, or lagged returns.
- The final refined model explains **over 80% of the variation** in log per capita NSDP.

---

## Instructions for Replication

1. Open the RStudio project file

2. Run the analysis:
- Execute `ma_etrx_assignment_code.R` **from top to bottom** to reproduce all results.

3. Refer to:
- `ma_etrx_assignment.pdf` for detailed methodology, diagnostics, and interpretation.

