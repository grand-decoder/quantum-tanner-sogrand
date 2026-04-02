# Simulate the performance of standard BP+OSD decoding for quantum Tanner codes.
# This script distributes the Monte Carlo simulations for one fixed maximal decoding iteration across multiple CPU cores
# whereas bp_simulation_parallelize_iter.py simulates various maximal decoding iterations in parallel.

# Simulation results are stored in the folder
# data/outputs/{yy_MM_dd}/{HH_mm_ss}_{code_name}_osd{use_osd}_iter{max_iter}_{file_prefix}.mat
# containing:
#   - stores the simulation parameters.
#   - p_list: array containing the p_physical values that were simulated.
#   - bler_deg_list: Logical block error rate (BLER).
#   - bler_list: same as bler_deg_list
#     but count a decoding as an error if not the exact error pattern is
#     recovered whereas bler_deg_list and bler_deg_per_iteration_list take
#     degeneracy into account.

import numpy as np
import scipy
import helper

from pathlib import Path

import os
import math

from concurrent.futures import ProcessPoolExecutor
import functools

def simulate(sim_idx, Hx_global, Hz_global, p_list, use_osd, max_iter, 
             osd_order, num_qubits, H_global, H_dual_symplectic, 
             max_simulations_per_core, max_errors_per_core):
    rng = np.random.default_rng()

    print(f"Process {sim_idx} started", flush=True)
    num_error_degenerated_list = np.zeros_like(p_list, dtype=int)
    num_error_list = np.zeros_like(p_list, dtype=int)
    num_decoded_list = np.zeros_like(p_list, dtype=int)

    for p_idx, p in enumerate(p_list):
        p = float(p)

        decoder_z, decoder_x = helper.create_decoders(Hx_global, Hz_global, p, use_osd, max_iter, osd_order)

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
        num_error_list[p_idx] = num_error
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
    max_simulations = 100
    # Maximal number of block errors before simulation terminates.
    max_errors = 10
    # Maximal number of decoding iterations
    max_iter = 100

    # If true, OSD with combination strategy using the LDPCv2 package is used 
    # as a post-processing step.
    use_osd = True
    # OSD order for post-processing.
    osd_order = 9

    # Physical error rates p_physical are simulated at p_count log-uniformly
    # spaced values over [10^p_start, 10^p_end].
    p_start = -3
    p_end = -0.8
    p_count = 15

    # file prefix to name simulation run.
    file_prefix = "test"


    ##########################
    # Simulation
    ##########################

    file_name_part = f"{code_name}_osd{use_osd}_iter{max_iter}_{file_prefix}"

    # -----------------
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
            use_osd=use_osd,
            max_iter=max_iter,
            osd_order=osd_order,
            num_qubits=num_qubits,
            H_global=H_global,
            H_dual_symplectic=H_dual_symplectic,
            max_simulations_per_core=max_simulations_per_core,
            max_errors_per_core=max_errors_per_core
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
        "max_iter": max_iter,
        "osd_order": osd_order,
        "use_osd": use_osd,
        "p_start": p_start,
        "p_end": p_end,
        "p_count": p_count,
    }
    scipy.io.savemat(file_name, save_dict)
