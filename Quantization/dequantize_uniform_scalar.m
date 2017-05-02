%-------------------------------------------------------------------------%

%Dequantize (reconstruct) an input vector that was uniformly quantized 
%using a uniform scalar quantizer.

%---- INPUTS ----

%quantized_vals: Quantization symbols (integers) produced at the quantizer.

%step_size: Uniform step size that was used for the quantizer. If step_size
%           is 0, then no quantization was applied at the encoder.

%---- OUTPUTS ----

%dequantized_vals: Values reconstructed from the input quantization
%                  symbols.


%-------------------------------------------------------------------------%

function dequantized_vals = dequantize_uniform_scalar(quantized_vals, step_size)

if step_size == 0
    %No quantization was applied at the encoder
    dequantized_vals = quantized_vals;
else
    %Compute the dequantized, or reconstructed, value corresponding to each
    %symbol produced by the quantizer. These reconstructed values should be
    %close to the original input values at the quantizer, but not exactly
    %the same due to the irreversible rounding operation at the quantizer. 
    dequantized_vals = quantized_vals.*step_size;
end

