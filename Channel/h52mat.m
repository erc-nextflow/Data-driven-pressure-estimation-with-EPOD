% convert HDF5 file to MAT file
HDF5Prefix='Fields/';
Var2Convert = {};
Target = 'Fields';

if ~exist(Target, 'dir')
    mkdir(Target);
end

HDF5List = dir([HDF5Prefix, '*.h5']);
NFile = size(HDF5List, 1);
fprintf('begining convert\n');
for iFile = 1:NFile
    fprintf(HDF5List(iFile).name);
    subname = HDF5List(iFile).name;
    tmp = strfind(subname,'.h5');
    subname(tmp(end):end) = [];
    clear tmp
    target_name = [Target,'/',subname,'.mat'];
    if isempty(Var2Convert)
        FlagEmptyVar = 1;
        fileinfo =...
            h5info([HDF5List(iFile).folder, '/', HDF5List(iFile).name]);
        Var2Convert = repmat({''},length(fileinfo.Datasets),1);
        for iSub = 1:length(fileinfo.Datasets)
            Var2Convert{iSub} = fileinfo.Datasets(iSub).Name;
        end
    end
    save2mat([HDF5List(iFile).folder,'/',HDF5List(iFile).name],...
        Var2Convert, target_name);
    if FlagEmptyVar == 1
        Var2Convert = {};
    end
    fprintf(repmat('\b',1,length(HDF5List(iFile).name)));
end
fprintf('ending convert\n');



HDF5Prefix='Fields_Testing/';
Var2Convert = {};
Target = 'Fields_Testing';

if ~exist(Target, 'dir')
    mkdir(Target);
end

HDF5List = dir([HDF5Prefix, '*.h5']);
NFile = size(HDF5List, 1);
fprintf('begining convert\n');
for iFile = 1:NFile
    fprintf(HDF5List(iFile).name);
    subname = HDF5List(iFile).name;
    tmp = strfind(subname,'.h5');
    subname(tmp(end):end) = [];
    clear tmp
    target_name = [Target,'/',subname,'.mat'];
    if isempty(Var2Convert)
        FlagEmptyVar = 1;
        fileinfo =...
            h5info([HDF5List(iFile).folder, '/', HDF5List(iFile).name]);
        Var2Convert = repmat({''},length(fileinfo.Datasets),1);
        for iSub = 1:length(fileinfo.Datasets)
            Var2Convert{iSub} = fileinfo.Datasets(iSub).Name;
        end
    end
    save2mat([HDF5List(iFile).folder,'/',HDF5List(iFile).name],...
        Var2Convert, target_name);
    if FlagEmptyVar == 1
        Var2Convert = {};
    end
    fprintf(repmat('\b',1,length(HDF5List(iFile).name)));
end
fprintf('ending convert\n');



HDF5Prefix='Probes/';
Var2Convert = {};
Target = 'Probes';

if ~exist(Target, 'dir')
    mkdir(Target);
end

HDF5List = dir([HDF5Prefix, '*.h5']);
NFile = size(HDF5List, 1);
fprintf('begining convert\n');
for iFile = 1:NFile
    fprintf(HDF5List(iFile).name);
    subname = HDF5List(iFile).name;
    tmp = strfind(subname,'.h5');
    subname(tmp(end):end) = [];
    clear tmp
    target_name = [Target,'/',subname,'.mat'];
    if isempty(Var2Convert)
        FlagEmptyVar = 1;
        fileinfo =...
            h5info([HDF5List(iFile).folder, '/', HDF5List(iFile).name]);
        Var2Convert = repmat({''},length(fileinfo.Datasets),1);
        for iSub = 1:length(fileinfo.Datasets)
            Var2Convert{iSub} = fileinfo.Datasets(iSub).Name;
        end
    end
    save2mat([HDF5List(iFile).folder,'/',HDF5List(iFile).name],...
        Var2Convert, target_name);
    if FlagEmptyVar == 1
        Var2Convert = {};
    end
    fprintf(repmat('\b',1,length(HDF5List(iFile).name)));
end
fprintf('ending convert\n');



HDF5Prefix='Probes_Testing/';
Var2Convert = {};
Target = 'Probes_Testing';

if ~exist(Target, 'dir')
    mkdir(Target);
end

HDF5List = dir([HDF5Prefix, '*.h5']);
NFile = size(HDF5List, 1);
fprintf('begining convert\n');
for iFile = 1:NFile
    fprintf(HDF5List(iFile).name);
    subname = HDF5List(iFile).name;
    tmp = strfind(subname,'.h5');
    subname(tmp(end):end) = [];
    clear tmp
    target_name = [Target,'/',subname,'.mat'];
    if isempty(Var2Convert)
        FlagEmptyVar = 1;
        fileinfo =...
            h5info([HDF5List(iFile).folder, '/', HDF5List(iFile).name]);
        Var2Convert = repmat({''},length(fileinfo.Datasets),1);
        for iSub = 1:length(fileinfo.Datasets)
            Var2Convert{iSub} = fileinfo.Datasets(iSub).Name;
        end
    end
    save2mat([HDF5List(iFile).folder,'/',HDF5List(iFile).name],...
        Var2Convert, target_name);
    if FlagEmptyVar == 1
        Var2Convert = {};
    end
    fprintf(repmat('\b',1,length(HDF5List(iFile).name)));
end
fprintf('ending convert\n');









function save2mat(HDF5Path, Var2Convert, TargetName)
for iSub = 1:length(Var2Convert)
    tmp = h5read(HDF5Path,['/',Var2Convert{iSub}]);
    localfcn(Var2Convert, iSub, tmp);
    if iSub == 1
        save(TargetName, Var2Convert{iSub});
    else
        save(TargetName, Var2Convert{iSub}, '-append');
    end
end
end

function localfcn(Var2Convert, iSub, tmp)
    assignin( 'caller', Var2Convert{iSub}, tmp);
end




