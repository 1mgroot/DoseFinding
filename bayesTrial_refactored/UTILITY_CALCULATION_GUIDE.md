# Utility Score Calculation Guide

## Overview

The utility score quantifies the risk-benefit tradeoff for each dose level in the Bayesian dose-finding trial. It combines information about toxicity (T), efficacy (E), and immune response (I) to provide a single numerical value for decision-making.

## Utility Table Structure

The utility table defines the value of each possible outcome combination:

```
Utility Table:
I=0 (No Immune Response):
  E=0, T=0: 0   E=1, T=0: 80 
  E=0, T=1: 0   E=1, T=1: 30 
I=1 (Immune Response):
  E=0, T=0: 10   E=1, T=0: 100 
  E=0, T=1: 0   E=1, T=1: 40 
```

**Interpretation:**
- **Highest utility (100)**: Immune response + Efficacy + No toxicity
- **High utility (80)**: Efficacy + No toxicity (no immune response)
- **Medium utility (40)**: Immune response + Efficacy + Toxicity
- **Low utility (30)**: Efficacy + Toxicity (no immune response)
- **Very low utility (10)**: Immune response + No efficacy + No toxicity
- **Zero utility (0)**: No efficacy outcomes

## Calculation Process

### Step 1: Extract Posterior Probabilities
For each dose level, extract the posterior probabilities:
- `π_I`: Probability of immune response
- `π_T|I=0`: Probability of toxicity given no immune response
- `π_T|I=1`: Probability of toxicity given immune response
- `π_E|I=0`: Probability of efficacy given no immune response
- `π_E|I=1`: Probability of efficacy given immune response

### Step 2: Calculate Probability Distributions
Convert probabilities to full distributions:
- `P(T=0|I=0) = 1 - π_T|I=0`, `P(T=1|I=0) = π_T|I=0`
- `P(T=0|I=1) = 1 - π_T|I=1`, `P(T=1|I=1) = π_T|I=1`
- `P(E=0|I=0) = 1 - π_E|I=0`, `P(E=1|I=0) = π_E|I=0`
- `P(E=0|I=1) = 1 - π_E|I=1`, `P(E=1|I=1) = π_E|I=1`

### Step 3: Calculate Expected Utilities
For each immune response state (I=0, I=1):

```
Expected Utility = Σ(Utility_Table[E,T,I] × P(E) × P(T))
```

This is calculated as:
```
Utility_I0 = sum(utility_table[,,1] * (p_E_given_I0 %o% p_T_given_I0))
Utility_I1 = sum(utility_table[,,2] * (p_E_given_I1 %o% p_T_given_I1))
```

### Step 4: Calculate Total Expected Utility
Combine utilities weighted by immune response probability:
```
Total Utility = (1 - π_I) × Utility_I0 + π_I × Utility_I1
```

## Example Calculation

**Dose 3 (Stage 1):**
- `π_I = 0.761` (76.1% immune response)
- `π_T|I=0 = 0.453`, `π_T|I=1 = 0.463`
- `π_E|I=0 = 0.449`, `π_E|I=1 = 0.851`

**Probability Distributions:**
- `P(T=0|I=0) = 0.547`, `P(T=1|I=0) = 0.453`
- `P(T=0|I=1) = 0.537`, `P(T=1|I=1) = 0.463`
- `P(E=0|I=0) = 0.551`, `P(E=1|I=0) = 0.449`
- `P(E=0|I=1) = 0.149`, `P(E=1|I=1) = 0.851`

**Expected Utilities:**
- `Utility_I0 = 25.76`
- `Utility_I1 = 62.28`

**Total Utility:**
- `Total = (1-0.761) × 25.76 + 0.761 × 62.28 = 53.55`

## Interpretation

### High Utility Scores (>50)
- Indicate favorable risk-benefit profiles
- Usually associated with high immune response and efficacy
- May still have some toxicity but balanced by benefits

### Medium Utility Scores (20-50)
- Moderate risk-benefit profiles
- May have lower efficacy or higher toxicity
- Still potentially acceptable depending on thresholds

### Low Utility Scores (<20)
- Poor risk-benefit profiles
- Usually due to low efficacy, high toxicity, or both
- Generally not recommended for further development

## Decision Making

The utility scores are used for:
1. **Adaptive Randomization**: Higher utility doses get more allocation probability
2. **Final Dose Selection**: The dose with highest utility among admissible doses is selected
3. **Trial Monitoring**: Track how utility evolves across stages

## Logging Output

The system now provides detailed logging showing:
- Utility table reference
- Individual dose calculations with all intermediate steps
- Summary table comparing all doses
- Probability distributions and expected utilities

This transparency helps in:
- Debugging calculations
- Understanding dose preferences
- Validating trial decisions
- Regulatory review and documentation 