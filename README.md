Title: Verification of Binegativity

Description: Python programs and reference results for testing two-qubit binegativity counterexamples under local amplitude damping and one-way LOCC, together with methods for reproducing the figures in the paper.

Paper: Local amplitude damping can increase two-qubit binegativity
Authors: Xiao-Lan Zong, Wei Song, and Ming Yang
Preprint: https://arxiv.org/abs/2610.34322
Program package version: 1.0
Release date: 2026-11-07

1.Purpose

This archive provides the verification programs and reference results used in the paper. The programs verify three specific counterexamples, perform finite-grid verification on a continuous-parameter counterexample family, and generate the curves in the paper for binegativity, negativity, and concurrence as functions of local amplitude damping. All input matrices and parameters have already been written into the programs. After installing the required Python dependencies, the programs can be run; no additional data files, network connection, or optimization solver are required.

The programs use exact rational arithmetic and floating-point arithmetic for verification. The content verified by each is described in Section 5.

2.File Description

README.txt
This description file, containing the title and introduction of the supplementary material.
verify_binegativity.py
The program used to verify the results and generate the figures.
verification_results.json
The reference results supplied with the package; the program does not read this file when running.
requirements.txt
Lists the NumPy and Matplotlib versions used in this verification.
verification_run.log
The run log of this verification, recording the execution command, start time, exit status, elapsed time, and the complete standard output and standard error output.
validation_report.json
Records the test environment and verification results.

No image files are included in the archive. They are generated after running the program; see below for details.

3.Runtime Environment and Usage

Use the dependency versions listed in requirements.txt. Python 3.11 or higher is required. The programs were tested on Linux with Python 3.12.14, NumPy 2.3.5, and Matplotlib 3.10.8. Detailed environment information is given in validation_report.json.

First extract the entire ZIP archive to a writable folder, then enter that folder and execute the following commands. An isolated virtual environment is recommended.

Linux or macOS:

python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python verify_binegativity.py

Windows, using the Python launcher:

py -3 -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe verify_binegativity.py

If the required dependencies are already installed, only the last command needs to be executed. The program does not depend on LaTeX or a graphical desktop environment; plotting uses Matplotlib's non-interactive Agg backend.

The program performs verification through Python assertions, so do not use -O or -OO when running, and do not set PYTHONOPTIMIZE to a nonzero value. If assertions have been disabled, the program will report an error and stop.

After a successful run, the program outputs a JSON object and exits with status code 0. The following files are saved in the folder where the script resides, independent of the current working directory:

verification_results.json
binegativity_damping.pdf
binegativity_damping.png

Among these, verification_results.json will overwrite the reference results supplied with the package. To preserve the original reference results, save a copy before running, or re-extract a fresh copy of the package for running. If the two image files already exist, they will also be overwritten.

The supplied run log and validation report record the test situation at release; the program does not regenerate these two files.

This test run took about 3 seconds; the actual time depends on the runtime environment. If an assertion fails or another exception occurs, the program exits with a nonzero status code. In this case, check the cause of the error rather than skipping the relevant verification.

4.Expected Results and Their Correspondence to the Paper

The following three data blocks in the JSON output correspond respectively to the three examples in the paper. Section titles are used here to indicate the corresponding locations, so that changes in section numbering after paper re-typesetting do not affect the correspondence.

local_amplitude_damping
This corresponds to the section "A Counterexample under Local Amplitude Damping" of the paper.
Both the input and output states are full-rank; z = 6/17, eta = 253/289.
Input N2: 54/47125
Output N2: 1248/1047625
Gain: 618/13619125
Ratio: 2704/2601

one_way_LOCC
This corresponds to the rational-number instance in the section "A Counterexample Family" of the paper.
p = 5/17, z = 5/9; the state in this instance is not full-rank.
Input N2: 9/221
Output N2: 20/459
Gain: 1/351
Ratio: 260/243
The probabilities of the two measurement outcomes are 565/1377 and 812/1377, respectively.

full_rank_one_way_LOCC
This corresponds to the appendix "A Full-Rank Instance of the One-Way Protocol".
Both the input and output states are full-rank; z = 4/7.
Input N2: 603/28730
Output N2: 872/38675
Gain: 1567/1005550
Ratio: 22672/21105

In each data block, the exact rational numbers in input_N2, output_N2, gain, ratio, and, when present, the probability fields, are stored as strings. minimum_eigenvalues gives the floating-point calculation results for the minimum eigenvalues of the input and output states, for auxiliary checking. Zero values in non-full-rank instances are normal results.

The remaining JSON fields record the 301 grid points of the continuous-parameter counterexample family, the maximum numerical deviation of the corresponding formulas, the 429 checks performed on the aforementioned input-state family, and the normalization convention used. In addition, the program computes derivative expressions at 101 points, and takes 1001 damping probability values to plot the curves.

All curves in the figure are normalized to the initial values of the corresponding physical quantities, and the positions corresponding to eta = 253/289 and N2(output)/N2(input) = 2704/2601 are marked.

5.Exact and Numerical Verification

For the three specific instances, the program uses fractions.Fraction to perform exact rational arithmetic, verifying the channel output matrices, the trace-1 condition, matrix symmetry, the positive semidefiniteness conditions determined by the block structure, the specified partial-transpose eigenvalues and their corresponding eigenvectors, and the exact value of binegativity and the strict positivity of the gain.

For the full-rank instances, the program exactly verifies that the diagonal entries and the determinants of the nontrivial blocks are strictly positive, thereby confirming that the matrices are positive definite. For the rational-number instance with p = 5/17, the program also exactly verifies the probabilities of the measurement branches and the binegativity values of the retained branches.

Floating-point calculations are used for further cross-checking, including: verifying the completeness condition of the Kraus operators and the channel action, directly computing binegativity with an eigenvalue solver and comparing it with its spectral formula, and comparing with negativity and concurrence. The finite-grid verification of the continuous-parameter counterexample family, the finite-grid verification of the derivative expressions, the verification of the aforementioned input-state family, and the curve plotting also use floating-point arithmetic.

For the full-rank one-way LOCC instance, the program performs numerical verification on the local measurement branches; the corresponding exact argument is given in the appendix of the paper.

The absolute tolerance used for the Kraus operator verification is 2e-15, and the absolute tolerance used for comparing the above eigenvalue calculations with the formulas is 2e-14. When run in different environments, the last few digits of the floating-point results, the display effect or metadata of the images may differ slightly, but the strings representing exact rational numbers should be completely identical.

The exact proof of full rank is based on the above rational-number verification, not on the minimum eigenvalues obtained by floating-point calculation. The finite-grid verification cannot replace the analytical proof covering the entire parameter interval in the paper, nor the derivative argument therein.

Conventions:
The ordering of the computational basis is |00>, |01>, |10>, |11>.
Gamma denotes partial transposition with respect to Bob's qubit, i.e., the second qubit.
The negative part is defined as X_- = (|X| - X)/2, which is a positive semidefinite operator.
N(rho) = ||rho^Gamma||1 - 1.
N2(rho) = Tr[(rho^Gamma)-] + 2 Tr[(((rho^Gamma)-)^Gamma)-].
The operators of each branch need not be normalized.
By positive homogeneity, adding the N2 values of these unnormalized branches gives the average binegativity weighted by the branch probabilities.
