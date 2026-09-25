
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function TAX = stepTaxonomy()
% STEPTAXONOMY  Group the registry steps for the step picker tree.
%   TAX = STEPTAXONOMY() returns an ordered struct array of categories, each a
%   pipeline stage the user thinks in. Every category carries an ordered list
%   of operations; an operation is either a single step or several provider
%   variants of the same operation (e.g. Run ICA -> FastICA/Infomax/Picard/
%   TESA), which the picker nests so the user chooses among providers.
%
%   Fields:
%     TAX(c).name                 category (stage) name
%     TAX(c).ops(o).name          operation name (the node when >1 variant)
%     TAX(c).ops(o).variants(v).step      exact registry step name (the value)
%     TAX(c).ops(o).variants(v).provider  who supplies it (drives the in-house
%                                          flag; AARATEP/nestapp are flagged)
%
%   This is design metadata, deliberately hand-ordered - it is NOT derivable
%   from the registry. Nothing checks that it covers the listed registry steps,
%   so a step added to the registry and not placed here is dropped from the
%   picker silently, with nothing failing. Add it in both or not at all.
%
%   See also: availableSteps, stepRegistry, stepInteractiveIcon

raw = {
 'Import & Session', {
   {'Load Data',                   {{'Load Data','EEGLAB'}}}
   {'Load Channel Location',       {{'Load Channel Location','EEGLAB'}}}
   {'Choose Data Set',             {{'Choose Data Set','EEGLAB'}}}
   {'Save New Set',                {{'Save New Set','EEGLAB'}}}
 }
 'Channels', {
   {'Remove un-needed Channels',   {{'Remove un-needed Channels','EEGLAB'}}}
   {'Detect Bad Channels',         {{'Detect Bad Channels (RANSAC)','clean_rawdata'}, ...
                                    {'Detect Bad Channels (TESA)','TESA'}}}
   {'Remove Bad Channels',         {{'Remove Bad Channels','EEGLAB'}}}
   {'Remove Bad Channels (manual)',{{'Remove Bad Channels (manual)','EEGLAB'}}}
   {'Interactive Channel Reject',  {{'Interactive Channel Reject (TESA)','TESA'}}}
   {'Interpolate Channels',        {{'Interpolate Channels','EEGLAB'}}}
 }
 'TMS Artifact', {
   {'Find TMS Pulses',             {{'Find TMS Pulses (TESA)','TESA'}, ...
                                    {'Find TMS Pulses (AARATEP)','AARATEP'}}}
   {'Remove TMS Artifacts',        {{'Remove TMS Artifacts (TESA)','TESA'}}}
   {'Fix TMS Pulse',               {{'Fix TMS Pulse (TESA)','TESA'}}}
   {'Interpolate Missing Data',    {{'Interpolate Missing Data (TESA)','TESA'}, ...
                                    {'Interpolate Missing Data (AR-Blend)','AARATEP'}}}
   {'Remove Decay Artifact',       {{'Remove Decay Artifact','AARATEP'}}}
   {'Fit Artifact Model',          {{'Fit Artifact Model (TESA)','TESA'}}}
   {'Find Artifacts EDM',          {{'Find Artifacts EDM (TESA)','TESA'}}}
   {'SSP-SIR',                     {{'SSP SIR','TESA'}}}
   {'Source-Informed Sensor Cleaning', {{'Source-Informed Sensor Cleaning (SOUND)','TESA'}}}
 }
 'Reference & Resample', {
   {'Re-Reference',                {{'Re-Reference','EEGLAB'}}}
   {'Re-Sample',                   {{'Re-Sample','EEGLAB'}}}
 }
 'Detrend', {
   {'Linear Detrend',              {{'De-Trend Epoch','EEGLAB'}, ...
                                    {'TESA De-Trend','TESA'}}}
   {'Robust Detrend',              {{'Robust Detrend (TESA)','TESA'}}}
   {'Robust Demean',               {{'Robust Demean (TESA)','TESA'}}}
 }
 'Filter', {
   {'Frequency Filter',            {{'Frequency Filter','firfilt'}, ...
                                    {'Frequency Filter (TESA)','TESA'}, ...
                                    {'Frequency Filter (CleanLine)','CleanLine'}}}
   {'Modified Bandpass Filter',    {{'Modified Bandpass Filter (TESA)','TESA'}}}
   {'Median Filter 1D',            {{'Median Filter 1D','TESA'}}}
 }
 'Automatic Cleaning', {
   {'Clean Artifacts',             {{'Clean Artifacts','clean_rawdata'}}}
   {'Automatic Cleaning Data',     {{'Automatic Cleaning Data','clean_rawdata'}}}
   {'Automatic Continuous Rejection', {{'Automatic Continuous Rejection','EEGLAB'}}}
 }
 'Epochs & Trials', {
   {'Epoching',                    {{'Epoching','EEGLAB'}}}
   {'Remove Baseline',             {{'Remove Baseline','EEGLAB'}}}
   {'Remove Bad Epoch',            {{'Remove Bad Epoch','EEGLAB'}}}
   {'Remove Bad Trials',           {{'Remove Bad Trials','EEGLAB'}}}
 }
 'ICA', {
   {'Run ICA',                     {{'Run ICA (FastICA)','FastICA'}, ...
                                    {'Run ICA (Infomax)','EEGLAB'}, ...
                                    {'Run ICA (Picard)','PICARD'}, ...
                                    {'Run TESA ICA','TESA'}}}
   {'Label ICA Components',        {{'Label ICA Components','ICLabel'}}}
   {'Flag ICA Components',         {{'Flag ICA Components for Rejection','ICLabel'}}}
   {'Remove Flagged ICA Components', {{'Remove Flagged ICA Components','EEGLAB'}}}
   {'Remove ICA Components',       {{'Remove ICA Components (TESA)','TESA'}}}
 }
 'TEP Analysis', {
   {'Extract TEP',                 {{'Extract TEP (TESA)','TESA'}}}
   {'Find TEP Peaks',              {{'Find TEP Peaks (TESA)','TESA'}}}
   {'TEP Peak Output',             {{'TEP Peak Output','TESA'}}}
 }
 'Utilities & Inspection', {
   {'Quality Gate',                {{'Quality Gate','nestapp'}}}
   {'Visualize EEG Data',          {{'Visualize EEG Data','EEGLAB'}}}
   {'Manual Command',              {{'Manual Command','nestapp'}}}
 }
};

TAX = struct('name', {}, 'ops', {});
for i = 1:size(raw, 1)
    cat.name = raw{i,1};
    ops = raw{i,2};
    O = struct('name', {}, 'variants', {});
    for j = 1:numel(ops)
        O(j).name = ops{j}{1};
        vc = ops{j}{2};
        V = struct('step', {}, 'provider', {});
        for k = 1:numel(vc)
            V(k).step     = vc{k}{1};
            V(k).provider = vc{k}{2};
        end
        O(j).variants = V;
    end
    cat.ops = O;
    TAX(end+1) = cat; %#ok<AGROW>
end
end
