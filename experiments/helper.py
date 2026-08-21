import numpy as np
import scipy
from ldpc import BpOsdDecoder, BpDecoder
from pathlib import Path
from datetime import datetime

def gen2par(G):
    """
    Converts a systematic generator matrix G = [I | P]
    to a parity check matrix H = [P^T | I].
    """
    rows, cols = G.shape

    # Validate systematic form: left k x k block should be identity
    if cols < rows:
        raise ValueError("Generator matrix G must have cols >= rows for systematic form")

    left_block = G[:, :rows].astype(int)
    if not np.array_equal(left_block % 2, np.eye(rows, dtype=int)):
        raise ValueError("G is not in systematic form [I | P]: left block is not identity")

    # Check if G is in systematic form [I | P]
    # We assume the Gaussian elimination in g2h ensured the left block is Identity.
    # The P matrix is the right part of G.
    k = rows
    n = cols
    P = G[:, k:]

    # H = [P.T | I_(n-k)]
    # Create identity matrix for the check bits
    eye_part = np.eye(n - k, dtype=int)

    # Transpose P
    P_T = P.T

    # Concatenate to form H
    H = np.hstack((P_T, eye_part))

    return H

def g2h(G):
    """
    Converts a Generator matrix G into a Parity Check matrix H.
    
    This function handles non-systematic codes and linearly dependent rows
    by performing Gaussian elimination with full pivoting (row and column swaps),
    converting to [P^T | I], and then reversing the column swaps.
    """
    # Ensure input is a numpy array of integers
    H = G.copy().astype(int)
    rows, cols = H.shape
    col_swaps = []
    
    # Loop over the diagonal
    # effectively min(rows, cols), but we might break early if rank is exhausted
    limit = min(rows, cols)
    
    for i in range(limit):
        # --- Step 1: Find a Pivot in the entire remaining block ---
        # The MATLAB code searches the whole submatrix H(i:rows, i:cols)
        # for the first non-zero entry.
        
        sub_matrix = H[i:, i:]
        
        # specific pivot selection to match MATLAB's column-major preference isn't 
        # strictly required for correctness, but we need *some* non-zero entry.
        # np.argwhere returns (row, col) indices of non-zero elements.
        non_zeros = np.argwhere(sub_matrix)
        
        if non_zeros.size == 0:
            # If the entire sub-matrix is zeros, we are done with elimination
            continue
            
        # Pick the first non-zero found. 
        # Note: Python scans row-by-row, MATLAB scans column-by-column. 
        # This might result in a different valid H than MATLAB, but both are correct.
        r_rel, c_rel = non_zeros[0]
        
        # Convert relative indices to absolute indices
        pivot_row = i + r_rel
        pivot_col = i + c_rel
        
        # --- Step 2: Swap to bring pivot to (i,i) ---
        
        # Swap Rows if needed
        if pivot_row != i:
            H[[i, pivot_row], :] = H[[pivot_row, i], :]
            
        # Swap Columns if needed
        if pivot_col != i:
            H[:, [i, pivot_col]] = H[:, [pivot_col, i]]
            col_swaps.append((i, pivot_col))
            
        # --- Step 3: Elimination (Gaussian) ---
        # We perform elimination on ALL rows (except i) to create a unit vector at col i.
        # This creates the Identity structure required for gen2par.
        
        # Identify rows that have a 1 in column 'i', excluding row 'i' itself
        col_vec = H[:, i]
        rows_to_xor = np.where(col_vec == 1)[0]
        rows_to_xor = rows_to_xor[rows_to_xor != i]
        
        if rows_to_xor.size > 0:
            # Vectorized XOR operation
            H[rows_to_xor, :] = np.bitwise_xor(H[rows_to_xor, :], H[i, :])

    # Remove zero rows (handle dependent rows)
    # This keeps rows that are NOT all zeros
    H = H[~np.all(H == 0, axis=1)]
    
    # Convert systematic G to H
    H_parity = gen2par(H)

    # Undo the column swaps on the resulting H matrix
    # Iterate in reverse order
    for idx in range(len(col_swaps) - 1, -1, -1):
        c1, c2 = col_swaps[idx]

        # Swap columns c1 and c2
        tmp = H_parity[:, c1].copy()
        H_parity[:, c1] = H_parity[:, c2]
        H_parity[:, c2] = tmp

    assert np.all(np.mod(H_parity @ np.transpose(G), 2) == 0)

    return H_parity

def swap_z_x(v):
    """
    Swaps the first half (Z) and second half (X) of a vector.
    """
    # Ensure input is a numpy array
    v = np.array(v)

    # Check if input is a 1D vector
    if v.ndim != 1:
        raise ValueError('Input must be a 1D vector.')

    n = v.size

    # Check if the number of entries is even
    if n % 2 != 0:
        raise ValueError('Input must have an even number of entries.')

    mid = n // 2

    # Swap the first and second half
    # Python slicing v[start:end] is equivalent to MATLAB v(start+1:end)
    return np.concatenate((v[mid:], v[:mid]))


def create_decoders(Hx, Hz, p, use_osd, max_iter, osd_order, bp_method):
    """
    Creates BP (or BP+OSD) decoders for Z and X errors.
    Returns (decoder_z, decoder_x) where decoder_z decodes Z errors from X stabilizers
    and decoder_x decodes X errors from Z stabilizers.
    """
    if use_osd:
        decoder_z = BpOsdDecoder(
            Hx,
            error_rate=2 * p / 3,
            bp_method=bp_method,
            max_iter=max_iter,
            ms_scaling_factor=0.625,
            osd_method="OSD_CS",
            osd_order=osd_order
        )
        decoder_x = BpOsdDecoder(
            Hz,
            error_rate=2 * p / 3,
            bp_method=bp_method,
            max_iter=max_iter,
            ms_scaling_factor=0.625,
            osd_method="OSD_CS",
            osd_order=osd_order
        )
    else:
        decoder_z = BpDecoder(
            Hx,
            error_rate=2 * p / 3,
            bp_method=bp_method,
            max_iter=max_iter,
            ms_scaling_factor=0.625,
        )
        decoder_x = BpDecoder(
            Hz,
            error_rate=2 * p / 3,
            bp_method=bp_method,
            max_iter=max_iter,
            ms_scaling_factor=0.625,
        )
    return decoder_z, decoder_x


def load_code(code_name, base_dir):
    """
    Loads a quantum Tanner code by name and returns its matrices.
    Returns (Hx_global, Hz_global, Hx_global, num_qubits, H_dual_symplectic).
    """
    if code_name == "[[250_10_15]]":
        Hx_file = "data/codes/[[250_10_15]]/X_stabs"
        Hz_file = "data/codes/[[250_10_15]]/Z_stabs"
    elif code_name == "[[200_10_10]]":
        Hx_file = "data/codes/[[200_10_10]]/X_stabs"
        Hz_file = "data/codes/[[200_10_10]]/Z_stabs"
    else:
        raise ValueError("Unknown code name")

    Hx_global = np.loadtxt(base_dir / Hx_file).astype(int)
    Hz_global = np.loadtxt(base_dir / Hz_file).astype(int)

    H_global = scipy.linalg.block_diag(Hx_global, Hz_global)
    num_qubits = Hx_global.shape[1]

    tmp = g2h(H_global)
    H_dual_symplectic = np.hstack((tmp[:, num_qubits:], tmp[:, :num_qubits]))

    return Hx_global, Hz_global, H_global, num_qubits, H_dual_symplectic


def run_decoding_trial(rng, decoder_x, decoder_z, p, num_qubits, H_global, n_x_checks, H_dual_symplectic):
    """
    Runs a single decoding trial under the depolarizing channel at error rate p.
    Returns (had_error, had_degenerate_error).
    """
    correct_error_gf4 = rng.choice([1, 2, 3, 4], p=[1 - p, p / 3, p / 3, p / 3], size=num_qubits)
    correct_error_x = ((correct_error_gf4 == 2) | (correct_error_gf4 == 4)).astype(np.uint8)
    correct_error_z = ((correct_error_gf4 == 3) | (correct_error_gf4 == 4)).astype(np.uint8)
    correct_error = np.concatenate((correct_error_x, correct_error_z))

    syndrome = np.mod(H_global @ swap_z_x(correct_error), 2)
    syndrome_z = syndrome[:n_x_checks].astype(np.uint8)
    syndrome_x = syndrome[n_x_checks:].astype(np.uint8)

    error_hat_x = decoder_x.decode(syndrome_x)
    error_hat_z = decoder_z.decode(syndrome_z)
    error_hat = np.concatenate((error_hat_x, error_hat_z))

    error_total = np.mod(error_hat + correct_error, 2)

    had_degenerate_error = not np.all(np.mod(H_dual_symplectic @ swap_z_x(error_total), 2) == 0)
    had_error = bool(np.any(error_total))

    return had_error, had_degenerate_error


def get_filename(prefix):
    # 1. Setup Base Paths
    current_dir = Path(__file__).resolve().parent.parent
    output_dir = current_dir / "data" / "output"

    # 2. Get current time
    now = datetime.now()

    # 3. Define the date-based subfolder path
    # Using Path objects directly instead of strings
    folder_path = output_dir / now.strftime("%y_%m_%d")
    timestamp = now.strftime("%H_%M_%S")

    # 4. Create the directory (parents=True handles data/output/YY_MM_DD all at once)
    folder_path.mkdir(parents=True, exist_ok=True)

    # 5. Construct the full filename path
    filename = f"{timestamp}_{prefix}.mat"
    full_path = folder_path / filename

    return full_path