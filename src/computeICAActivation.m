function act2D = computeICAActivation(EEG)
% COMPUTEICAACTIVATION  Return 2-D ICA activations (nComp x nSamples).
%   act2D = COMPUTEICAACTIVATION(EEG)
%
%   Uses EEG.icaact when populated; otherwise computes the textbook
%   (icaweights * icasphere) * data formula. Falls back to the full
%   channel set when EEG.icachansind is missing or empty so the function
%   is safe to call on partially-built EEG structs.
if isfield(EEG, 'icaact') && ~isempty(EEG.icaact)
    act2D = reshape(EEG.icaact, size(EEG.icaact, 1), []);
    return
end
chanInd = 1:EEG.nbchan;
if isfield(EEG, 'icachansind') && ~isempty(EEG.icachansind)
    chanInd = EEG.icachansind;
end
data2D = reshape(EEG.data(chanInd, :, :), numel(chanInd), []);
act2D  = (EEG.icaweights * EEG.icasphere) * data2D;
end
