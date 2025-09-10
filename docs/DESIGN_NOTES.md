# Design Notes for Bayesian Adaptive Trial Simulation

This document serves as a planning space for complex coding tasks related to the Bayesian adaptive dose-finding trial simulation.

---

## COMPREHENSIVE TRIAL DESIGN COMPLIANCE ANALYSIS

### Task: Review Current Implementation Against TRIAL_DESIGN.md Specifications

#### 1. TRIAL_DESIGN.md vs Current Implementation Comparison

**✅ FULLY IMPLEMENTED:**

**Section 1.1 - First Stage (S=1):**
- ✅ Equal randomization to all dose levels
- ✅ No adaptive elements in first stage
- ✅ Implementation: `alloc_probs <- rep(1/length(trial_config$dose_levels), length(trial_config$dose_levels))`

**Section 1.2 - Adaptive Stages (S=2 to S):**
- ✅ Interim analysis based on posterior probabilities
- ✅ Admissible set criteria (Safety, Efficacy, Immune Response)
- ✅ Adaptive randomization within admissible set
- ✅ Early termination when admissible set is empty

**Section 2 - Utility Function:**
- ✅ Expected utility calculation: `U(d_j) = Σ_yᵢ p(Yᵢ=yᵢ | d_j) (Σ_yₜ Σ_yₑ w(yₜ, yₑ, yᵢ) · p(Yₜ=yₜ | Yᵢ=yᵢ, d_j) · p(Yₑ=yₑ | Yᵢ=yᵢ, d_j))`
- ✅ Utility table implementation in config.R
- ✅ Conditional independence assumption implemented

**Section 3 - Admissible Set Definition:**
- ✅ Safety: `P(πₜⱼ ≤ Φₜ | Data) > cₜ`
- ✅ Efficacy: `P(πₑⱼ > Φₑ | Data) > cₑ`
- ✅ Immune Response: `P(πᵢⱼ > Φᵢ | Data) > cᵢ`
- ✅ Implementation: `get_admissible_set()` function

**Section 4.2 - Randomization Probabilities:**
- ✅ `rⱼ = U(dⱼ) / Σₖ∈ₐ U(dₖ)`
- ✅ Implementation: `adaptive_randomization()` function

**Section 5.1 - Selection Condition:**
- ✅ Final OD selected only if trial completes without early termination

**Section 5.2 - OD Definition:**
- ✅ `OD = argmaxⱼ∈ₐ [Û(dⱼ)]`
- ✅ Implementation: `select_final_od_with_poc()` function

**Section 7.1 - Stage-wise Execution:**
- ✅ Step 1: Equal randomization (Stage 1)
- ✅ Step 2: Interim analysis
- ✅ Step 3: Adaptive randomization
- ✅ Step 4: Early termination check
- ✅ Step 5: Final selection

**❌ MISSING OR DIFFERENT:**

**Section 4.1 - Posterior Probability γⱼ Calculation:**
- ❌ **Missing**: `γⱼ = Pr(OD = dⱼ | Dₙ)` calculation
- ❌ **Missing**: Function to calculate posterior probability that each admissible dose is the Optimal Dose
- ❌ **Note**: This is only needed for control arm management (Section 4.3)

**Section 4.3 - Control Arm Management:**
- ❌ **Missing**: Control arm assignment formula: `Control Assignment = rₒ · α · min(γₘₐₓ, 1/Tₐᵢ)`
- ❌ **Missing**: Control arm parameters: `rₒ`, `α`, `Tₐᵢ`
- ❌ **Missing**: `γₘₐₓ = max{γⱼ, j∈A}` calculation
- ❌ **Note**: Marked as future task

**Section 5.3 - Probability of Correct Selection (PoC):**
- ⚠️ **PARTIALLY IMPLEMENTED**: Current PoC implementation differs from TRIAL_DESIGN.md specification
- ❌ **Different**: TRIAL_DESIGN.md specifies: `Pr(Πᵢ < δ Πᵢⱼ | Dₙ) > Cₚₒc`
- ⚠️ **Current Implementation**: Uses normal approximation and combined efficacy measure
- ❌ **Missing**: Proper Πᵢ and Πᵢⱼ parameter calculation as specified
- ❌ **Missing**: Calibrated Cₚₒc threshold determination process

**Section 6.1 - Calibration:**
- ❌ **Missing**: Calibrated Cₚₒc threshold determination
- ❌ **Missing**: Parameter tuning process for optimal performance
- ❌ **Missing**: Calibration framework for thresholds and cutoffs

**Section 6.2 - Error Rate Control:**
- ❌ **Missing**: Familywise Type I Rate control at 0.05
- ❌ **Missing**: Power target specification (80% - 90%)
- ❌ **Missing**: Error rate control mechanisms

**Section 6.3 - Dose-Response Curves:**
- ⚠️ **PARTIALLY IMPLEMENTED**: Current simulation uses fixed probability matrices
- ❌ **Missing**: Accommodation for different underlying dose-response curves
- ❌ **Missing**: Flat immune response curve option
- ❌ **Missing**: Configurable dose-response curve specifications

**Section 7.2 - Computational Requirements:**
- ✅ **Implemented**: Bayesian posterior calculations using conjugate priors
- ✅ **Implemented**: Isotonic regression for monotonicity constraints
- ✅ **Implemented**: Utility calculation with detailed logging
- ❌ **Missing**: Adaptive randomization with control arm management

#### 2. Configuration Analysis

**Current Config Parameters:**
- ✅ `phi_T`, `c_T` - Toxicity thresholds
- ✅ `phi_E`, `c_E` - Efficacy thresholds  
- ✅ `phi_I`, `c_I` - Immune response thresholds
- ✅ `c_poc`, `delta_poc` - PoC parameters
- ✅ `enable_early_termination`, `log_early_termination` - Early termination
- ❌ **Missing**: Control arm parameters (`r_o`, `alpha`, `T_AI`)
- ❌ **Missing**: Calibration parameters
- ❌ **Missing**: Error rate control parameters

#### 3. Priority Implementation Tasks

**High Priority (Core Functionality):**
1. **PoC Implementation Correction** - Align with TRIAL_DESIGN.md Section 5.3 specification
2. **γⱼ Calculation** - For future control arm implementation
3. **Configuration Parameter Addition** - Add missing parameters

**Medium Priority (Advanced Features):**
4. **Control Arm Management** - Complete Section 4.3 implementation
5. **Calibration Framework** - Implement Section 6.1
6. **Error Rate Control** - Implement Section 6.2

**Low Priority (Optimization):**
7. **Dose-Response Curve Flexibility** - Implement Section 6.3
8. **Advanced Logging Features** - Enhanced debugging capabilities

#### 4. Critical Issues to Address

**1. PoC Implementation Mismatch:**
- **Issue**: Current PoC uses normal approximation instead of proper Bayesian calculation
- **Impact**: May not provide correct probability of correct selection
- **Solution**: Implement proper `Pr(Πᵢ < δ Πᵢⱼ | Dₙ)` calculation

**2. Missing Calibration:**
- **Issue**: No calibration process for thresholds
- **Impact**: Trial performance may be suboptimal
- **Solution**: Implement calibration framework

**3. Control Arm Preparation:**
- **Issue**: No γⱼ calculation for control arm management
- **Impact**: Cannot implement control arm in future
- **Solution**: Implement γⱼ calculation function

---

## CURRENT TASK: Enhanced Adaptive Randomization Implementation

### Task: Implement Enhanced Adaptive Randomization According to TRIAL_DESIGN.md Section 4

#### 1. Current State Analysis

**What Exists:**
- Basic utility-based adaptive randomization in `adaptive_randomization()` function
- Utility calculation using `get_expected_utility()` 
- Admissible set identification with `get_admissible_set()`
- Early termination logic (COMPLETED)
- PoC implementation (COMPLETED)

**What's Missing (TRIAL_DESIGN.md Section 4):**
- **Posterior probability γ_j calculation**: `γ_j = Pr(OD = d_j | D_n)` for each admissible dose (for control arm management)
- **Control arm management**: Implementation of control arm assignment formula
- **γ_max calculation**: Required for control arm management

**TRIAL_DESIGN.md Requirements:**
1. **Section 4.1**: Calculate posterior probability γ_j that each admissible dose is the Optimal Dose (for control arm)
2. **Section 4.2**: Randomization probabilities are already correct: `r_j = U(d_j) / Σ_k∈A U(d_k)` ✅ IMPLEMENTED
3. **Section 4.3**: Control arm management (future task)

#### 2. Required Changes

**Files to Modify:**
- `dose_decision.R`: Enhance `adaptive_randomization()` function
- `config.R`: Add parameters for γ-based randomization
- `main.R`: Update workflow to use enhanced randomization

**New Functions to Implement:**
- `calculate_gamma_probabilities()`: Calculate γ_j for each admissible dose
- `enhanced_adaptive_randomization()`: Combine utility and γ-based randomization

#### 3. Implementation Steps

**Step 1: Implement γ_j Calculation**
- Files: `dose_decision.R`
- Functions: `calculate_gamma_probabilities(admissible_set, posterior_summaries, config)`
- Mathematical Basis: `γ_j = Pr(OD = d_j | D_n)` where OD is the dose with highest utility
- Expected outcome: Posterior probability that each admissible dose is the optimal dose

**Step 2: Implement γ_j Calculation for Control Arm**
- Files: `dose_decision.R`
- Functions: `calculate_gamma_probabilities(admissible_set, posterior_summaries, config)`
- Mathematical Basis: `γ_j = Pr(OD = d_j | D_n)` for control arm management
- Expected outcome: Posterior probability that each admissible dose is the optimal dose

**Step 3: Prepare for Control Arm Management**
- Files: `config.R`
- Parameters: Add control arm parameters (r_o, α, T_AI) for future implementation
- Expected outcome: Configuration ready for control arm implementation

**Step 4: Update Main Workflow**
- Files: `main.R`
- Functions: Ensure enhanced randomization is used in trial simulation
- Expected outcome: Trial uses enhanced adaptive randomization

**Step 5: Add Logging and Validation**
- Files: `dose_decision.R`
- Functions: Add detailed logging of γ calculations and randomization decisions
- Expected outcome: Clear visibility into randomization logic

#### 4. Mathematical Implementation Details

**γ_j Calculation for Control Arm:**
```
γ_j = Pr(OD = d_j | D_n) = Pr(U(d_j) > U(d_k) for all k ≠ j | D_n)
```

**Control Arm Assignment (Future):**
```
Control Assignment = r_o * α * min(γ_max, 1/T_AI)
```
where γ_max = max{γ_j, j∈A}

#### 5. Testing Strategy

**Unit Tests:**
- Test γ_j calculation with known posterior distributions
- Test control arm parameter calculations
- Validate mathematical correctness of probability calculations

**Integration Tests:**
- Test complete trial simulation with γ_j calculation
- Verify that current utility-based randomization works correctly
- Prepare for future control arm integration

**Validation Criteria:**
- γ_j probabilities sum to 1 across admissible doses
- Randomization probabilities are non-negative and sum to 1
- Enhanced randomization produces reasonable allocation patterns

#### 6. Documentation Updates

**Files to Update:**
- `CODE_MAP.md`: Document new functions
- `simulation_notebook.qmd`: Add examples of enhanced randomization
- `TRIAL_DESIGN.md`: Verify implementation matches design specifications

#### 7. Risk Assessment

**Potential Issues:**
- Computational complexity of γ_j calculation
- Numerical stability with small probabilities
- Integration with control arm management

**Fallback Plans:**
- Maintain current utility-based randomization as primary method
- γ_j calculation only for control arm when implemented
- Implement robust numerical methods for probability calculations

**Validation Steps:**
- Compare with theoretical expectations
- Test with edge cases (single admissible dose, equal utilities)
- Verify convergence properties over multiple stages

---

## Implementation Plan Template

When you provide complex tasks, I will use this template to create detailed plans:

### Task: [Task Name]

#### 1. Current State Analysis
- What exists in the current codebase
- What functionality is missing
- Dependencies and constraints

#### 2. Required Changes
- Files that need modification
- New functions to implement
- Configuration updates needed

#### 3. Implementation Steps
1. **Step 1**: [Specific action]
   - Files: [list of files]
   - Functions: [list of functions]
   - Expected outcome: [description]

2. **Step 2**: [Specific action]
   - Files: [list of files]
   - Functions: [list of functions]
   - Expected outcome: [description]

#### 4. Testing Strategy
- Unit tests to write
- Integration tests needed
- Validation criteria

#### 5. Documentation Updates
- Files to update
- New documentation needed

#### 6. Risk Assessment
- Potential issues
- Fallback plans
- Validation steps

---

## Priority Implementation Tasks

Based on the design analysis, here are the recommended implementation priorities:

### High Priority
1. **γ_j Calculation Implementation** - For control arm management
2. **Control Arm Parameter Preparation** - Configuration setup
3. **Validation of Current Randomization** - Ensure TRIAL_DESIGN.md compliance

### Medium Priority
4. **Control Arm Management** - Future task
5. **Additional Calibration Features** - Optimization improvements

### Low Priority
6. **Advanced Logging Features** - Debugging and analysis

---

## Code Structure Guidelines

When implementing new features, I will follow these guidelines:

### Configuration Management
- Add new parameters to `config.R`
- Use descriptive parameter names
- Include parameter validation

### Function Design
- Single responsibility principle
- Clear input/output specifications
- Comprehensive error handling
- Detailed logging for debugging

### Testing Strategy
- Unit tests for mathematical functions
- Integration tests for trial workflows
- Validation against known results

### Documentation
- Update relevant README files
- Add function documentation
- Include usage examples

---

This document will be updated with specific implementation plans when you provide complex tasks. Each plan will include detailed steps, file modifications, and testing strategies for your review before implementation.
