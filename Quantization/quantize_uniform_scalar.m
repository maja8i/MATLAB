%-------------------------------------------------------------------------%

%Uniformly quantize a given input vector, using a uniform scalar quantizer.

%---- INPUTS ----

%input: Input data that is to be quantized.

%step_size: Uniform step size that will be used for the quantizer. Values 
%           for step_size should ideally be powers of 2. If step_size is 0, 
%           then apply no quantization.  

%---- OUTPUTS ----

%quantized_vals: Quantization symbols (integers) produced by the quantizer.

%-------------------------------------------------------------------------%

function quantized_vals = quantize_uniform_scalar(input, step_size)    

if step_size == 0
    %Do not quantize
    quantized_vals = input;
else
    %Compute the quantization symbols (integers)
    quantized_vals = round(input./step_size);   %Round to the nearest integer
end

    
