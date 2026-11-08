
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2023-2026 Aref Pariz and Wesley Dunne.
% Part of nestapp; see the LICENSE file for full terms.
function selected = selectDataTree(startFolder, exts)
% SELECTDATATREE  Modal folder/file browser for queuing data files.
%   selected = SELECTDATATREE(startFolder, exts) opens a modal window that
%   shows the folder tree beneath a chosen parent as a checkbox tree, with
%   a path filter so files spread across many subject folders can be picked
%   in one place. Returns a cell array of the full paths of every checked
%   file, or {} if the user cancels.
%
%   startFolder : folder to open into (the parent is browsed from here).
%   exts        : glob patterns for loadable files, default the formats
%                 nestapp understands. BrainVision .eeg/.vmrk companions
%                 are intentionally not listed - only the .vhdr header is a
%                 selectable leaf and the companions load from it.
%
%   Design notes:
%     - Folders render lazily: a folder's children are built only when it
%       is first expanded, so huge trees open instantly.
%     - The set of checked FILES is the source of truth (a containers.Map
%       keyed on full path), not the widget. The tree is just a view, so a
%       file stays checked even when filtered out of sight or never
%       materialised. Checking a folder cascades to every file beneath it,
%       including ones not yet expanded.
%     - The filter matches the path relative to the parent, so typing
%       "SPL", "_1_", or "pre" narrows by anything in folder or file names.

    if nargin < 1, startFolder = ''; end
    if nargin < 2 || isempty(exts)
        exts = {'*.set', '*.vhdr', '*.cdt', '*.cnt'};
    end

    % ---- shared state (closed over by the nested callbacks) -------------
    selected    = {};
    parentDir   = '';
    checkedSet  = containers.Map('KeyType', 'char', 'ValueType', 'logical');
    allFiles    = {};      % cached recursive full paths under parentDir
    allRelLower = {};      % their relative paths, lowercased, for filtering
    scanned     = false;   % is the allFiles cache valid for parentDir?
    prevKeys    = {};      % node keys checked at the last change event

    % ---- build the window ----------------------------------------------
    fig = appWindow('Name', 'Select Data Files', ...
        'Position', [200 140 760 600], 'Resize', 'off');
    % Keep the browser above the main window (not WindowStyle='modal',
    % whose close button is unreliable).
    try
        fig.WindowStyle = 'alwaysontop';
    catch
    end
    % Close handling: every exit path simply DELETES the figure, and the
    % bottom of this function blocks on waitfor(fig) (which returns the
    % moment the figure is destroyed). This avoids uiwait/uiresume, whose
    % close-on-X path could leave the window up and soft-lock the app. The
    % X uses a plain delete so it can never be vetoed by a callback error;
    % "selected" keeps its default {} (cancel) when the window is X-ed.
    fig.CloseRequestFcn = @(src, ~) delete(src);

    uilabel(fig, 'Text', 'Parent folder:', 'Position', [15 562 85 22]);
    parentField = uieditfield(fig, 'text', 'Editable', 'off', ...
        'Position', [100 560 520 26], 'Value', '(none selected)');
    uibutton(fig, 'Text', 'Choose...', 'Position', [630 560 115 28], ...
        'ButtonPushedFcn', @(~,~) onChoose());

    uilabel(fig, 'Text', 'Filter:', 'Position', [15 522 45 22]);
    searchField = uieditfield(fig, 'text', 'Position', [60 520 350 26], ...
        'ValueChangedFcn', @(~,~) onFilter());
    uibutton(fig, 'Text', 'Check shown', 'Position', [420 519 105 28], ...
        'ButtonPushedFcn', @(~,~) setShownChecked(true));
    uibutton(fig, 'Text', 'Uncheck shown', 'Position', [530 519 115 28], ...
        'ButtonPushedFcn', @(~,~) setShownChecked(false));
    uibutton(fig, 'Text', 'Clear filter', 'Position', [650 519 95 28], ...
        'ButtonPushedFcn', @(~,~) clearFilter());

    tree = uitree(fig, 'checkbox', 'Position', [15 60 730 450], ...
        'NodeExpandedFcn', @(~,e) onExpand(e), ...
        'CheckedNodesChangedFcn', @(~,e) onCheckChanged(e));
    % Make the check handler NON-REENTRANT. onCheckChanged briefly yields the
    % UI thread (a uifigure property round-trip), and while it is yielded the
    % next checkbox event would interrupt and re-enter it - each re-entry
    % nesting another O(nChecked) pass, so per-click cost doubled per event
    % (0.07s -> 0.14 -> ... -> 22s) and the browser appeared to freeze under
    % steady clicking. 'off' makes events queue and run one-at-a-time instead.
    tree.Interruptible = 'off';
    tree.BusyAction    = 'queue';

    countLabel = uilabel(fig, 'Text', '0 files checked', ...
        'Position', [15 18 430 22], 'FontColor', [0.2 0.2 0.2]);
    uibutton(fig, 'Text', 'Cancel', 'Position', [530 14 100 30], ...
        'ButtonPushedFcn', @(~,~) onCancel());
    uibutton(fig, 'Text', 'Use checked', 'Position', [640 14 105 30], ...
        'BackgroundColor', [0.20 0.55 0.20], 'FontColor', [1 1 1], ...
        'ButtonPushedFcn', @(~,~) onUse());

    % Seed with the start folder if it is usable, else prompt immediately.
    if ~isempty(startFolder) && isfolder(startFolder)
        setParent(startFolder);
    else
        onChoose();
    end

    waitfor(fig);   % blocks until a button or the X deletes the figure

    % Remember the browse ROOT (the data folder), not the subfolder the
    % chosen files sit in, so the next open returns here and the Preferences
    % "Data folder" setting is respected rather than clobbered. Only on a
    % confirmed selection - cancelling/X leaves the setting untouched.
    if ~isempty(selected) && ~isempty(parentDir) && isfolder(parentDir)
        setpref('nestapp', 'lastDataFolder', parentDir);
    end

    % ====================================================================
    %  Folder choice / (re)build
    % ====================================================================
    function onChoose()
        base = parentDir;
        if isempty(base), base = startFolder; end
        p = uigetdir(base, 'Select the parent folder containing your data');
        if ~isequal(p, 0)
            setParent(p);
        end
        if isvalid(fig); figure(fig); end   % uigetdir can bury the modal
    end

    function setParent(p)
        parentDir         = p;
        parentField.Value = p;
        checkedSet        = containers.Map('KeyType', 'char', 'ValueType', 'logical');
        allFiles          = {};
        allRelLower       = {};
        scanned           = false;
        searchField.Value = '';
        rebuildTree('');
        updateCount();
    end

    function rebuildTree(filterText)
        delete(tree.Children);
        prevKeys = {};
        if isempty(parentDir), return; end
        if isempty(strtrim(filterText))
            addChildren(tree, parentDir);     % lazy: top level only
        else
            buildFilteredTree(strtrim(filterText));
        end
        refreshCheckState();
        prevKeys = nodeKeys(tree.CheckedNodes);
    end

    % Add the immediate sub-folders and loadable files of FOLDERPATH under
    % PARENTNODE (a uitree or a uitreenode). Folders get a placeholder child
    % so they show an expand arrow but cost nothing until opened.
    function addChildren(parentNode, folderPath)
        d   = dir(folderPath);
        sub = d([d.isdir]);
        sub = sub(~ismember({sub.name}, {'.', '..'}));
        [~, order] = sort(lower({sub.name}));
        sub = sub(order);
        for i = 1:numel(sub)
            subPath = fullfile(folderPath, sub(i).name);
            fnode = uitreenode(parentNode, 'Text', sub(i).name, ...
                'NodeData', struct('type', 'folder', 'path', subPath, ...
                                   'populated', false));
            uitreenode(fnode, 'Text', '...', ...
                'NodeData', struct('type', 'dummy', 'path', '', ...
                                   'populated', true));
        end
        files = loadableFilesIn(folderPath, exts);
        for i = 1:numel(files)
            [~, nm, ex] = fileparts(files{i});
            uitreenode(parentNode, 'Text', [nm ex], ...
                'NodeData', struct('type', 'file', 'path', files{i}, ...
                                   'populated', true));
        end
    end

    function onExpand(e)
        node = e.Node;
        data = node.NodeData;
        if ~strcmp(data.type, 'folder') || data.populated
            return
        end
        delete(node.Children);            % drop the placeholder
        addChildren(node, data.path);
        data.populated = true;
        node.NodeData  = data;
        refreshCheckState();
        prevKeys = nodeKeys(tree.CheckedNodes);
    end

    % ====================================================================
    %  Filtering
    % ====================================================================
    function onFilter()
        rebuildTree(searchField.Value);
    end

    function clearFilter()
        searchField.Value = '';
        rebuildTree('');
    end

    function buildFilteredTree(query)
        if ~ensureScanned(); return; end
        keep    = contains(allRelLower, lower(query));
        relHits = allRelLower(keep);
        absHits = allFiles(keep);
        % Recreate the folder hierarchy for just the matching files.
        folderNodes = containers.Map('KeyType', 'char', 'ValueType', 'any');
        for i = 1:numel(absHits)
            parts  = strsplit(relHits{i}, '/');
            parent = tree;
            for k = 1:numel(parts)-1                 % ancestor folders
                keySoFar = strjoin(parts(1:k), '/');
                if isKey(folderNodes, keySoFar)
                    parent = folderNodes(keySoFar);
                else
                    fp = fullfile(parentDir, strrep(keySoFar, '/', filesep));
                    parent = uitreenode(parent, 'Text', parts{k}, ...
                        'NodeData', struct('type', 'folder', 'path', fp, ...
                                           'populated', true));
                    folderNodes(keySoFar) = parent;
                end
            end
            uitreenode(parent, 'Text', parts{end}, ...
                'NodeData', struct('type', 'file', 'path', absHits{i}, ...
                                   'populated', true));
        end
        try
            expand(tree);
        catch
        end
    end

    function ok = ensureScanned()
        ok = true;
        if scanned, return; end
        dlg = uiprogressdlg(fig, 'Title', 'Scanning folders', ...
            'Message', 'Looking for data files...', ...
            'Indeterminate', 'on', 'Cancelable', 'on');
        cl = onCleanup(@() closeValid(dlg));
        prog = @(n, folder) reportScan(dlg, n);
        [allFiles, allRelLower] = scanLoadableFiles(parentDir, exts, prog);
        allRelLower = lower(allRelLower);
        scanned = ~dlg.CancelRequested;
        ok = scanned;
    end

    % ====================================================================
    %  Bulk check / uncheck of what is currently shown
    % ====================================================================
    function setShownChecked(tf)
        if isempty(parentDir), return; end
        q = strtrim(searchField.Value);
        if isempty(q)
            % No filter: act on the whole tree (needs the full list).
            if ~ensureScanned(); return; end
            targets = allFiles;
            if tf && numel(targets) > 200
                sel = uiconfirm(fig, sprintf(['Check all %d files under ' ...
                    'this folder?'], numel(targets)), 'Check all', ...
                    'Options', {'Check all', 'Cancel'}, 'DefaultOption', 2);
                if ~strcmp(sel, 'Check all'); return; end
            end
        else
            if ~ensureScanned(); return; end
            targets = allFiles(contains(allRelLower, lower(q)));
        end
        for i = 1:numel(targets)
            if tf
                checkedSet(targets{i}) = true;
            elseif isKey(checkedSet, targets{i})
                remove(checkedSet, targets{i});
            end
        end
        refreshCheckState();
        prevKeys = nodeKeys(tree.CheckedNodes);
        updateCount();
    end

    % ====================================================================
    %  Check-state sync (model <-> widget)
    % ====================================================================
    % User toggled a checkbox: translate the change into the file-level
    % truth in checkedSet. Folder nodes cascade to every descendant file,
    % but a folder going to the *indeterminate* (partially checked) state is
    % ignored wholesale - its individual leaf deltas carry the real change.
    function onCheckChanged(e)
        curr   = nodeKeys(e.CheckedNodes);
        indet  = nodeKeys(e.IndeterminateCheckedNodes);
        added   = setdiff(curr, prevKeys);
        removed = setdiff(prevKeys, curr);
        applyDelta(added, true, indet);
        applyDelta(removed, false, indet);
        prevKeys = curr;
        updateCount();
    end

    function applyDelta(keys, makeChecked, indeterminateKeys)
        for i = 1:numel(keys)
            key  = keys{i};
            kind = key(1);                 % 'f' file, 'd' dir
            path = key(3:end);             % strip "f:" / "d:"
            if kind == 'f'
                setFileChecked(path, makeChecked);
            else
                % A folder that is now partially checked is handled by its
                % leaves; don't sweep all descendants in that case.
                if ~makeChecked && any(strcmp(indeterminateKeys, key))
                    continue
                end
                kids = collectLoadableUnder(path);
                for j = 1:numel(kids)
                    setFileChecked(kids{j}, makeChecked);
                end
            end
        end
    end

    function setFileChecked(path, tf)
        if tf
            checkedSet(path) = true;
        elseif isKey(checkedSet, path)
            remove(checkedSet, path);
        end
    end

    % Point the widget's CheckedNodes at exactly the materialised file nodes
    % whose path is in checkedSet. Parent folders show checked/indeterminate
    % automatically from their children.
    function refreshCheckState()
        fileNodes = gatherFileNodes(tree);
        keepMask  = false(1, numel(fileNodes));
        for i = 1:numel(fileNodes)
            keepMask(i) = isKey(checkedSet, fileNodes(i).NodeData.path);
        end
        keep = fileNodes(keepMask);
        if isempty(keep)
            tree.CheckedNodes = [];   % setter wants [] for "none", not typed-empty
        else
            tree.CheckedNodes = keep;
        end
    end

    % ====================================================================
    %  Finish
    % ====================================================================
    function onUse()
        if checkedSet.Count == 0
            uialert(fig, 'Check at least one file first.', 'Nothing selected');
            return
        end
        selected = sort(keys(checkedSet));
        delete(fig);            % releases waitfor and closes the window
    end

    function onCancel()
        selected = {};
        delete(fig);
    end

    function updateCount()
        n = checkedSet.Count;
        word = 'files'; if n == 1, word = 'file'; end
        countLabel.Text = sprintf('%d %s checked', n, word);
    end

    % ====================================================================
    %  Small helpers
    % ====================================================================
    % Recursively gather every loadable file under a folder on disk. Used
    % when a folder checkbox cascades to files that may not be materialised.
    function files = collectLoadableUnder(folder)
        % Files a folder-checkbox should cascade to. With a filter active,
        % cascade ONLY to the files currently visible under the folder (the
        % ones matching the filter) - not everything on disk - so checking a
        % folder queues exactly what the user sees. Unfiltered, cascade to
        % every loadable file beneath it.
        query = strtrim(searchField.Value);
        if isempty(query)
            files = scanLoadableFiles(folder, exts);
            return
        end
        if ~ensureScanned(); files = {}; return; end
        prefix      = [folder filesep];
        underFolder = startsWith(allFiles, prefix);
        matches     = contains(allRelLower, lower(query));   % allRelLower is lowercased
        files       = allFiles(underFolder & matches);
    end

    % Depth-first list of all materialised file nodes in the tree.
    function nodes = gatherFileNodes(root)
        nodes = matlab.ui.container.TreeNode.empty;
        stack = root.Children(:)';
        while ~isempty(stack)
            n = stack(1);
            stack(1) = [];
            nd = n.NodeData;
            if isstruct(nd) && isfield(nd, 'type') && strcmp(nd.type, 'file')
                nodes(end+1) = n; %#ok<AGROW>
            end
            if ~isempty(n.Children)
                stack = [n.Children(:)', stack]; %#ok<AGROW>
            end
        end
    end

    % Stable string keys for a set of nodes: "f:<path>" / "d:<path>".
    function ks = nodeKeys(nodes)
        ks = {};
        for i = 1:numel(nodes)
            nd = nodes(i).NodeData;
            if ~isstruct(nd) || ~isfield(nd, 'type'); continue; end
            switch nd.type
                case 'file',   ks{end+1} = ['f:' nd.path]; %#ok<AGROW>
                case 'folder', ks{end+1} = ['d:' nd.path]; %#ok<AGROW>
            end
        end
    end
end

% -- file-scope helpers (no shared state) ----------------------------------

function cancelled = reportScan(dlg, n)
    dlg.Message = sprintf('Found %d data file(s)...', n);
    cancelled   = dlg.CancelRequested;
end

function closeValid(dlg)
    if isvalid(dlg); close(dlg); end
end
