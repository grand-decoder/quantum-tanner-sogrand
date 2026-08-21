# Simulate the performance of standard BP+Relay decoding for quantum Tanner codes.

# Simulation results are stored in the folder
# data/outputs/{yy_MM_dd}/{HH_mm_ss}_{code_name}_relay_{file_prefix}.mat
# containing:
#   - stores the simulation parameters.
#   - p_list: array containing the p_physical values that were simulated.
#   - bler_deg_list: Logical block error rate (BLER).
#   - bler_list: same as bler_deg_list
#     but count a decoding as an error if not the exact error pattern is
#     recovered whereas bler_deg_list and bler_deg_per_iteration_list take
#     degeneracy into account.

import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg
import scipy
import helper

from pathlib import Path
from datetime import datetime

import os
import math

from concurrent.futures import ProcessPoolExecutor
import functools
import relay_bp

def simulate(sim_idx, Hx_global, Hz_global, p_list, num_qubits, H_global, H_dual_symplectic,
             max_simulations_per_core, max_errors_per_core, gamma0, pre_iter, num_sets, set_max_iter,
             gamma_dist_interval, stop_nconv):
    rng = np.random.default_rng()

    print(f"Process {sim_idx} started", flush=True)

    num_error_degenerated_list = np.zeros_like(p_list, dtype=int)
    num_error_list = np.zeros_like(p_list, dtype=int)
    num_decoded_list = np.zeros_like(p_list, dtype=int)

    for p_idx, p in enumerate(p_list):
        p = float(p)

        decoder_z = relay_bp.RelayDecoderF32(
            Hx_global,
            error_priors=2 * p / 3 * np.ones(Hx_global.shape[1]),
            gamma0=gamma0,
            pre_iter=pre_iter,
            num_sets=num_sets,
            set_max_iter=set_max_iter,
            gamma_dist_interval=gamma_dist_interval,
            stop_nconv=stop_nconv,
        )

        decoder_x = relay_bp.RelayDecoderF32(
            Hz_global,
            error_priors=2 * p / 3 * np.ones(Hx_global.shape[1]),
            gamma0=gamma0,
            pre_iter=pre_iter,
            num_sets=num_sets,
            set_max_iter=set_max_iter,
            gamma_dist_interval=gamma_dist_interval,
            stop_nconv=stop_nconv,
        )

        num_error_degenerated = 0
        num_error = 0
        num_decoded = 0

        while num_decoded < max_simulations_per_core and num_error_degenerated < max_errors_per_core:
            num_decoded += 1

            had_error, had_degenerate_error = helper.run_decoding_trial(
                rng, decoder_x, decoder_z, p, num_qubits,
                H_global, Hx_global.shape[0], H_dual_symplectic
            )

            if had_degenerate_error:
                num_error_degenerated += 1
            if had_error:
                num_error += 1

        num_error_degenerated_list[p_idx] = num_error_degenerated
        num_error_list[p_idx]= num_error
        num_decoded_list[p_idx] = num_decoded

    return num_decoded_list, num_error_list, num_error_degenerated_list


if __name__ == "__main__":
    #########################
    # Simulation parameters
    #########################

    # Name of the quantum Tanner code. Currently, the quantum tanner codes from
    # [3] are supported, which can be downloaded from https://github.com/RebKatRad/qTanner or via setup_data.m.
    # [3] R. K. Radebold, S. D. Bartlett, and A. C. Doherty, "Explicit instances of quantum Tanner codes," Nov. 14, 2025, arXiv: arXiv:2508.05095. doi: 10.48550/arXiv.2508.05095.
    code_name = "[[250_10_15]]" # [[250_10_15]] or [[200_10_10]]

    # Maximal number of codewords before the simulation terminates.
    max_simulations = 10000
    # Maximal number of block errors before simulation terminates.
    max_errors = 200

    # Relay parameters (see [6] and relay-bp Python package for details):

    gamma0 = 0.0700  # Uniform memory weight for the first set
    pre_iter = 80  # Maximal BP iterations for the first set
    num_sets = 60  # Maximal number of sets
    set_max_iter = 60  # Maximal BP iterations per relay set (except for the first set)
    # Sampling interval for the memory strength for each leg except for the first one.
    gamma_dist_interval_start = -0.1375  # Set the uniform distribution range for disordered memory weight selection
    gamma_dist_interval_end = 0.4975
    stop_nconv = 1  # Number of relay solutions to find before stopping (the best will be selected)

    p_start = -2
    p_end = -0.8
    p_count = 8

    file_prefix = "test"

    ##########################
    # Simulation
    ##########################

    file_name_part = f"{code_name}_relay_{file_prefix}"

    # -----------------
    gamma_dist_interval = (gamma_dist_interval_start, gamma_dist_interval_end)

    p_list = np.logspace(p_start, p_end, num=p_count)

    current_dir = Path(__file__).resolve().parent.parent
    print(current_dir)

    Hx_global, Hz_global, H_global, num_qubits, H_dual_symplectic = helper.load_code(code_name, current_dir)

    num_cores = os.cpu_count() // 2
    print(f"Detected {num_cores} cores. Spawning one task per core.")

    sim_idx_list = list(range(num_cores))

    max_simulations_per_core = math.ceil(max_simulations / num_cores)
    max_errors_per_core = math.ceil(max_errors / num_cores)

    # Use ProcessPoolExecutor to run across CPU cores
    # max_workers=None will use all available cores
    with ProcessPoolExecutor(max_workers=None) as executor:
        # Create a partial function with fixed arguments
        worker = functools.partial(simulate,
            Hx_global=Hx_global,
            Hz_global=Hz_global,
            p_list=p_list,
            num_qubits=num_qubits,
            H_global=H_global,
            H_dual_symplectic=H_dual_symplectic,
            max_simulations_per_core=max_simulations_per_core,
            max_errors_per_core=max_errors_per_core,
            gamma0=gamma0,
            pre_iter=pre_iter,
            num_sets=num_sets,
            set_max_iter=set_max_iter,
            gamma_dist_interval=gamma_dist_interval,
            stop_nconv=stop_nconv,
        )

        # Map the p_list to the workers
        results = list(executor.map(worker, sim_idx_list))

    print("finished")
    num_decoded_list_total = np.zeros((len(p_list)))
    num_error_list_total = np.zeros((len(p_list)))
    num_error_degenerated_list_total = np.zeros((len(p_list)))

    for idx, (num_decoded_list, num_error_list, num_error_degenerated_list) in enumerate(results):
        num_decoded_list_total += num_decoded_list
        num_error_list_total += num_error_list
        num_error_degenerated_list_total += num_error_degenerated_list

    bler_list = num_error_list_total / num_decoded_list_total
    bler_deg_list = num_error_degenerated_list_total / num_decoded_list_total

    file_name = helper.get_filename(file_name_part)
    save_dict = {
        "p_list": p_list,
        "bler_list": bler_list,
        "bler_deg_list": bler_deg_list,
        "number_decoded": num_decoded_list_total,
        "number_errors": num_error_list_total,
        "number_degenerated_errors": num_error_degenerated_list_total,
        "code_name": code_name,
        "max_simulations": max_simulations,
        "max_errors": max_errors,
        "p_start": p_start,
        "p_end": p_end,
        "p_count": p_count,
        "gamma0": gamma0,
        "pre_iter": pre_iter,
        "num_sets": num_sets,
        "set_max_iter": set_max_iter,
        "gamma_dist_interval": gamma_dist_interval,
        "stop_nconv": stop_nconv,
    }
    scipy.io.savemat(file_name, save_dict)
