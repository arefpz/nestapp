function act2D = computeICAActivation(EEG)
% COMPUTEICAACTIVATION  Return 2-D ICA activations (nComp x nSamples).
%   act2D = COMPUTEICAACTIVATION(EEG)
if ~isempty(EEG.icaact)
    act2D = reshape(EEG.icaact, size(EEG.icaact,1), []);
else
    data2D = reshape(EEG.data(EEG.icachansind,:,:), numel(EEG.icachansind), []);
    act2D  = (EEG.icaweights * EEG.icasphere) * data2D;
end
end
