# Discrete-functions-in-cryptography
A repository dedicated to the course discrete functions in cryptography. It implements a set of functions to search for characteristics and properties of Boolean functions and mappings.

# Implemented functions
## procedures.py
> The main file with functions for analyzing Boolean functions and mappings

#### **Algebraic Property Analysis**:
1. **deg** — Computes the degree of a function.  
2. **algebraic_immunity_1**, **algebraic_immunity_2** — Determine algebraic immunity.  
3. **Ann** — Returns the set of annihilators for a function.  

---

#### **Correlation and Spectral Properties**:
1. **correlative_immunity_K2** — Evaluates the order of correlation immunity.  
2. **cross_correlation**, **cross_correlation_def** — Compute cross-correlation between functions.  
3. **cross_correlation_k** — Finds the order of uncorrelatedness.  
4. **wht**, **inv_wht** — Walsh-Hadamard Transform (WHT) and inverse WHT.
5. **nonlinearity** — Calculates the nonlinearity of a Boolean function
   
---

#### **Linear and Differential Cryptanalysis**:
1. **DDT** — Generates a Difference Distribution Table.  
2. **LAT** — Generates a Linear Approximation Table.  

---

#### **Statistical Properties**:
1. **hamming_weight** — Computes Hamming weight (or distance between functions).  
2. **is_balanced** — Checks if a function is balanced.  
3. **m_balanced** — Evaluates the m-balanced property (strong equidistribution).  

---

#### **Utility Functions**:
1. **next_ham** — Generates the next number with the same Hamming weight (Gosper’s algorithm).  
2. **tensor_product**, **mat_mul** — Matrix operations (tensor product, multiplication).  
3. **scalar_product** — Scalar product of vectors.  

## format.py
> Additional functions of working with representations of function, substituting values in functions

#### **Representation Conversions**:
1. **a_to_x** — Converts an argument into a PDNF(Principle Disjunctive Normal Form) term.  
2. **a_to_v** — Converts a number into a binary vector.  
3. **vv_to_pdnf** — Vector of Values (VoV) → PDNF.  
4. **pdnf_to_vv** — PDNF → Vector of Values.  
5. **vv_to_anf** — VoV → ANF (Algebraic Normal Form).  
6. **anf_to_vv** — ANF → VoV.
7. **sub_a** — substitutes the value of a in ANF

## tests.py
> Some example use cases 
- Representation conversion tests (**format_test**).  
- Correlation property validation (**cross_correlation_test**).  
- S-box analysis (**DDT_and_LAT_test**).  
- Generating functions with specific properties (**correlative_immunity_test**, **algebraic_immunity_test**).  

## gf2.py
> Solves system of linear equations in GF(2) for algebraic_immunity_2 function
