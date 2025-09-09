# Trial Design: Mathematical Details

This document outlines the mathematical framework of the Bayesian adaptive dose-finding trial design.

## 1. Trial Stages and Randomization Strategy

### 1.1 First Stage (S=1)
- **Initial Randomization**: Patients (C₁) are equally randomized to J or (J+1) dose levels
- **No Adaptive Elements**: All dose levels receive equal allocation in the first stage

### 1.2 Adaptive Stages (S=2 to S)
- **Interim Analysis**: Based on interim data from previous stages, update the admissible set A
- **Admissible Set Criteria**: Doses in set A must satisfy criteria for:
  - **Safety**: Toxicity rate below threshold
  - **Activity**: Immune response above threshold  
  - **Efficacy**: Efficacy rate above threshold
- **Adaptive Randomization**: Patients (Cₛ) are adaptively randomized to doses dⱼ within admissible set A
- **Early Termination**: If admissible set A becomes empty, trial terminates early with no Optimal Dose (OD) selected

### 1.3 Early Termination Rule
- **Condition**: Admissible set A = ∅ at any stage
- **Outcome**: No Optimal Dose (OD) selected, trial ends

## 2. Utility Function: Risk-Benefit Tradeoff

The core of the decision-making process is the utility function, which quantifies the risk-benefit tradeoff for each dose. The utility of a dose dⱼ is the expected utility over all possible outcomes for toxicity (T), efficacy (E), and immune response (I).

The utility score w(yₜ, yₑ, yᵢ) for a given outcome combination is defined in the `utility_table` in `config.R`. The expected utility U(dⱼ) is calculated as:

U(dⱼ) = Σₓₜ₌₀¹ Σₓₑ₌₀¹ Σₓᵢ₌₀¹ w(yₜ, yₑ, yᵢ) · p(Yₜ=yₜ, Yₑ=yₑ, Yᵢ=yᵢ | dⱼ)

Assuming conditional independence of toxicity and efficacy given the immune response, this simplifies to:

U(dⱼ) = Σₓᵢ₌₀¹ p(Yᵢ=yᵢ | dⱼ) (Σₓₜ₌₀¹ Σₓₑ₌₀¹ w(yₜ, yₑ, yᵢ) · p(Yₜ=yₜ | Yᵢ=yᵢ, dⱼ) · p(Yₑ=yₑ | Yᵢ=yᵢ, dⱼ))

## 3. Admissible Set Definition

A dose dⱼ is considered admissible if it satisfies the following criteria based on the posterior probabilities of toxicity, efficacy, and immune response:

- **Safety**: The posterior probability that the toxicity rate (πₜⱼ) is below a certain threshold (Φₜ) must be greater than a certainty cutoff (cₜ).
  P(πₜⱼ ≤ Φₜ | Data) > cₜ

- **Efficacy**: The posterior probability that the efficacy rate (πₑⱼ) is above a certain threshold (Φₑ) must be greater than a certainty cutoff (cₑ).
  P(πₑⱼ > Φₑ | Data) > cₑ

- **Immune Response**: The posterior probability that the immune response rate (πᵢⱼ) is above a certain threshold (Φᵢ) must be greater than a certainty cutoff (cᵢ).
  P(πᵢⱼ > Φᵢ | Data) > cᵢ

## 4. Adaptive Randomization Mechanism

### 4.1 Posterior Probability Calculation
For each admissible dose dⱼ, calculate the posterior probability γⱼ that it is the Optimal Dose:

γⱼ = Pr(OD = dⱼ | Dₙ), for j ∈ A

where Dₙ represents the observed data up to stage n.

### 4.2 Randomization Probabilities
The allocation probability for each admissible dose dⱼ is proportional to its expected utility U(dⱼ):

rⱼ = U(dⱼ) / Σₖ∈ₐ U(dₖ)

### 4.3 Control Arm Management
When a control arm is present, patient assignment follows:

Control Assignment = rₒ · α · min(γₘₐₓ, 1/Tₐᵢ)

where:
- rₒ: Baseline randomization ratio for control
- α: Calibration parameter
- γₘₐₓ = max{γⱼ, j∈A}: Maximum posterior probability among admissible doses
- Tₐᵢ: Trial arm index parameter

## 5. Final Optimal Dose (OD) Selection

### 5.1 Selection Condition
The final Optimal Dose (OD) is selected only if the trial completes without early termination.

### 5.2 OD Definition
The Optimal Dose is the dose dⱼ within the admissible set A that maximizes the estimated utility Û(dⱼ):

OD = argmaxⱼ∈ₐ [Û(dⱼ)]

### 5.3 Probability of Correct Selection (PoC)
Before final selection, a Proof of Concept (PoC) check is performed:

Pr(Πᵢ < δ Πᵢⱼ | Dₙ) > Cₚₒc

where:
- Πᵢ and Πᵢⱼ: Parameters related to efficacy/toxicity
- δ: Predefined threshold
- Dₙ: Observed data up to stage n
- Cₚₒc: Pre-specified probability threshold for correct selection

## 6. Statistical Design Considerations

### 6.1 Calibration
- **Calibrated Cₚₒc**: The PoC threshold is determined through calibration process
- **Parameter Tuning**: All thresholds and cutoffs are calibrated for optimal performance

### 6.2 Error Rate Control
- **Familywise Type I Rate**: Controlled at 0.05
- **Power Target**: 80% - 90% power, or lower depending on trial constraints

### 6.3 Dose-Response Curves
The simulation accommodates different underlying dose-response curves:
- Flat immune response curve
- Toxicity curve  
- Efficacy curve

## 7. Implementation Notes

### 7.1 Stage-wise Execution
1. **Stage 1**: Equal randomization to all dose levels
2. **Interim Analysis**: Update admissible set based on posterior probabilities
3. **Adaptive Randomization**: Allocate patients based on utility scores
4. **Early Termination Check**: Terminate if admissible set is empty
5. **Final Selection**: Choose OD with highest utility from admissible set

### 7.2 Computational Requirements
- Bayesian posterior calculations using conjugate priors
- Isotonic regression for monotonicity constraints
- Utility calculation with detailed logging
- Adaptive randomization with control arm management

