# Efficient Soft-Output Guessing for Enhanced Quantum Tanner Code Decoding

Subject to license: **[GRAND Codebase Non-Commercial Academic Research Use License 021722.pdf](GRAND%20Codebase%20Non-Commercial%20Academic%20Research%20Use%20License%20021722.pdf)**

This repository contains the official implementation of the paper 
**Efficient Soft-Output Guessing for Enhanced Quantum Tanner Code Decoding** [1]. 

**Paper:** [Physical Review A](https://doi.org/10.1103/hnds-fvry) | [arXiv](https://arxiv.org/abs/2603.18318)

It includes both MATLAB and Python environments to simulate SOGRAND and standard BP+OSD decoding of Quantum Tanner Codes, along with simulation data.

## Prerequisites

- **MATLAB** (with MEX compiler configured)
- **Python 3.10+**
- Tested on: Linux, Windows 11

## Setup
**1. Compile the MEX files in MATLAB:**
```matlab
cd src/decoder
mex SOGRAND.c
```

**2. Set up Python environment:**

Note: The Python environment is required for the Python and MATLAB simulations as MATLAB calls the ldpc package
for OSD post-processing.
```bash
python -m venv .venv
# Activate environment (Linux/macOS): source .venv/bin/activate
# Activate environment (Windows): .\.venv\Scripts\activate
pip install -r requirements.txt
```

**3. Download Quantum Tanner Code Data:**

To test the decoder on the quantum Tanner codes described in [3], obtain the corresponding data from the 
[qTanner GitHub repository](https://github.com/RebKatRad/qTanner). 
Move the repository’s `data` folder into this project’s data directory and rename it to `codes` so the final path is: 
`quantum_tanner_sogrand/data/codes/`.

This setup can also be done automatically by running:
```matlab
setup_data
```

## Run Code
### Simulate SOGRAND Decoding
To simulate SOGRAND decoding of Quantum Tanner Codes, run the following in MATLAB:
```matlab
setup_env  % Once at MATLAB startup 
simulate_quantum_tanner_code
```

`setup_env` adds required MATLAB paths and configures the Python environment for OSD post-processing.

### Simulate Standard BP and BP+OSD Decoding
To simulate the baseline BP and BP+OSD decoding, run the following in your terminal:
```bash
# Ensure your python virtual environment is activated
cd experiments
python bp_simulation.py
```
To simulate performance over maximal iteration numbers, run:
```bash
# Ensure your python virtual environment is activated
cd experiments
python bp_simulation_parallelize_iter.py
```
Note: In contrast to the results in [1], the simulation does not apply message normalization if the sum product decoder is used.
Message normalization for the sum product decoder requires modifications of the `ldpc` package which we do not provide in this repository.

### Reproduce Plots from the Paper
To reproduce the figures from [1], run the following in MATLAB:

```matlab
% (run in MATLAB):
bler_performance_curves % Figure 3 in [1]
convergence_behavior % Figure 4 in [1]
relay_bler_performance % Figure 6 in [1]
relay_different_sets % Figure 7 in [1]
hyperparameter_optimization_relay % Figure 8 in [1]
```

## Citations
- The following should be cited in association with results from this code:
  > [1] L. Rapp, M. Médard, E. Tang, and K. R. Duffy, "Efficient Soft-Output Guessing for Enhanced Quantum Tanner Code Decoding," *Phys. Rev. A* **114**, 052441 (2026) [DOI](https://doi.org/10.1103/hnds-fvry) | [arXiv](https://arxiv.org/abs/2603.18318)
  ```bibtex
  @article{rapp2026efficient,
    title = {Efficient soft-output guessing for enhanced quantum Tanner code decoding},
    author = {Rapp, Lukas and M{\'e}dard, Muriel and Tang, Eugene and Duffy, Ken R.},
    journal = {Phys. Rev. A},
    volume = {114},
    issue = {5},
    pages = {052441},
    numpages = {11},
    year = {2026},
    month = {Nov},
    publisher = {American Physical Society},
    doi = {10.1103/hnds-fvry},
    url = {https://link.aps.org/doi/10.1103/hnds-fvry}
  }
  ```

- **SOGRAND Decoder:** The decoder uses a modified version of the SOGRAND decoder to incorporate non-zero syndrome during decoding.
The original code can be found [here](https://github.com/kenrduffy/SOGRAND-C):
  > [2] P. Yuan, M. Médard, K. Galligan, and K. R. Duffy, "Soft-Output (SO) GRAND and Iterative Decoding to Outperform LDPC Codes," *IEEE Trans. Wireless Commun.*, **24**, 3386 (2025). [DOI](https://doi.org/10.1109/TWC.2025.3530880).
  ```bibtex
  @ARTICLE{10852599,
    author={Yuan, Peihong and Médard, Muriel and Galligan, Kevin and Duffy, Ken R.},
    journal={IEEE Transactions on Wireless Communications}, 
    title={Soft-Output (SO) GRAND and Iterative Decoding to Outperform LDPC Codes}, 
    year={2025},
    volume={24},
    number={4},
    pages={3386-3399},
    keywords={Codes;Parity check codes;Maximum likelihood decoding;Product codes;Noise;Accuracy;Redundancy;Standards;Complexity theory;Wireless communication;GRAND;concatenated codes;product codes;LDPC;generalized LDPC},
    doi={10.1109/TWC.2025.3530880}}
  ```

## Acknowledgments
This repository utilizes excellent open-source tools and datasets from the quantum error correction community.

- **Quantum Tanner Codes:** [1] evaluates the SOGRAND decoder on the quantum Tanner codes introduced in [3].
To reproduce these experiments, the code data can be downloaded from the [qTanner GitHub repository](https://github.com/RebKatRad/qTanner) or fetched automatically using the `setup_data.m` script. 
These quantum Tanner codes remain the property of their respective authors.
  > [3] R. K. Radebold, S. D. Bartlett, and A. C. Doherty, "Explicit instances of quantum Tanner codes," arXiv:2508.05095 (2025) [arXiv](https://arxiv.org/abs/2508.05095).
  
- **LDPC Package:** The `ldpc` Python package is used to evaluate the BP+OSD baseline and the OSD post-processing used with the SOGRAND decoder. 
Since the current PyPI release does not yet provide a standalone OSD decoder, we instead use the `standalone_osd` branch of the 
[ldpc GitHub repository](https://github.com/quantumgizmos/ldpc) (see `requirements.txt` for details).
  > [4] J. Roffe, LDPC: Python tools for low density parity check codes (2022). [PyPI](https://pypi.org/project/ldpc/).

  > [5] J. Roffe, D. R. White, S. Burton, and E. Campbell, "Decoding across the quantum low-density parity-check code landscape," *Phys. Rev. Res.* **2**, 043423 (2020).

- **RelayBP Package:** The `relay-bp` Python package is used to evaluate the BP+Relay baseline. [relay-bp GitHub repository](https://github.com/trmue/relay)
  > [6]  T. Müller, T. Alexander, M. E. Beverland, M. Bühler, B. R. Johnson, T. Maurer, and D. Vandeth, "Improved belief propagation is sufficient for real-time decoding of quantum memory," arXiv:2506.01779 (2025) [arXiv](https://arxiv.org/abs/2506.01779).


## Note
Although the SOGRAND/OBRGRAND implementation provided here is not parallelized, the algorithm is highly parallelizable. E.g.
> [7] A. Riaz, Y. Alperen, F. Ercan, W. An, J. Ngo, K. Galligan, M. Médard, K. R. Duffy, and R. Yazicigil, "A sub-0.8-pJ/bit universal soft detection decoder using ORBGRAND,” *IEEE J. Solid-State Circuits* (2025).

> [8] E. Kizilates, A. Riaz, A. Bali, J. Feng, M. Médard, K. R. Duffy, and R. T. Yazicigil, "Soft output threshold-guided CRC decoding with SOGRAND in 40nm CMOS," *IEEE CICC* (2026).

Moreover, GLDPC decoding is highly parallelizable, where all component codes can be decoded in parallel, 
but this is not exploited here. For further details on GRAND, see: [https://www.granddecoder.mit.edu/](https://www.granddecoder.mit.edu/).
