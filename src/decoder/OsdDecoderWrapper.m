classdef OsdDecoderWrapper < handle
% OSDDECODERWRAPPER  MATLAB wrapper around the OsdDecoder of the Python ldpc package.
% Initializes and calls the Python OSD decoder with the appropriate
% channel probabilities derived from a physical error rate.
%
% This class requires that the Python virtual environment .venv with the external ldpc
% package installed is set up.
%
% Note: As the pip version of the ldpc package does not provide a standalone
% OSD decoder yet, the "standalone_osd" branch of the ldpc repository must be
% used. (See README.md and requirements.txt for further details.)

    properties
        DecoderObject  % Stores the Python OsdDecoder instance
        H              % Parity check matrix
        OSDOrder       % OSD order
    end

    methods
        % --- Constructor ---
        function obj = OsdDecoderWrapper(H, p_physical, osd_order)
        % OSDDECODERWRAPPER  Constructs the wrapper and initializes the Python OSD decoder.
        %
        % Inputs:
        %   H          - Binary parity-check matrix
        %   p_physical - Physical qubit error rate
        %   osd_order  - OSD order (number of bit-flip candidates to evaluate)
            obj.H = H;
            obj.OSDOrder = int32(osd_order);
            
            % Setup Python-compatible variables
            error_rate = 2 * p_physical / 3;
            num_bits = int32(size(H, 2));
            channel_probs = py.numpy.multiply(error_rate, py.numpy.ones(num_bits));
            
            % Initialize the Python decoder
            obj.DecoderObject = py.ldpc.osd_decoder.OsdDecoder(...
                py.numpy.array(obj.H), ...
                pyargs(...
                    'osd_order', obj.OSDOrder, ...
                    'osd_method', "OSD_CS", ...
                    'channel_probabilities', py.numpy.array(channel_probs) ...
                )...
            );
        end
        
        % --- Decoding Method ---
        function result = decode(obj, syndrome, llrs)
        % DECODE  Runs the OSD decoder for a given syndrome and channel LLRs.
        %
        % Inputs:
        %   syndrome - Binary syndrome vector
        %   llrs     - Log-likelihood ratios from the channel
        %
        % Output:
        %   result - Estimated error vector as a MATLAB column vector
            % Format syndrome for Python (int64, reshaped)
            py_syndrome = py.numpy.array(syndrome, pyargs('dtype', 'int64')).reshape(int32(-1));
            
            % Format LLRs
            py_llrs = py.numpy.array(llrs).reshape(int32(-1));
            
            % Call the Python decode method
            py_result = obj.DecoderObject.decode(py_syndrome, py_llrs);
            
            % Convert back to MATLAB double and ensure column vector
            result = double(py_result);
            result = col(result); 
        end
    end
end