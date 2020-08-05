%% setup workspace
wf = 'D:\job\ice'; addpath(wf); % path to iceopts.m, icenv.m, chkfile.m, tree.csv
iceopts(wf,false,false);        % warnings on (T/F), 'use' mode on (T/F)
clearvars -except a* wf;        % adapt as needed and preferred

%% find range of cell-line namespace from cell-line annotations file
assert(exist('ice000CellLineAnnotation.mat','file')==2, ...
             'Cell-line annotation MAT file required to fix column number.');
load ice000CellLineAnnotation.mat db*;
mxc = size(dbCellLineAnno,2);
clear db*;

%% load Mills proteomics data (antibody array) and extract headers
DM1 = chkfile('CCLE_RPPA_20181003.csv','Mills proteomics',',','string',true);
rnC = fields(DM1); rnC = string(rnC(2:(end-3)));
rnC(extractBefore(rnC,2)=='x') = extractAfter(rnC(extractBefore(rnC,2)=='x'),'x');
cn1 = DM1{:,DM1.Properties.VariableNames{1,1}}; cn1 = arxcellsyn(cn1); % cell lines are rows
ctm = numel(find(contains(cn1,'No CELL_SAMPLE_ID found')));
assert(ctm==0,'Achilles name not produced for %i input names.',ctm);
tm1 = (double(strrep(DM1{:,2:end},'NA','NaN')))'; % transpose
clear DM1;

%% process cell-line identifiers and shape data accordingly
[cn1s,cs] = sort(cn1); tm1 = tm1(:,cs); % initial column sort
[cn1u,cfo,cbi] = unique(cn1s); % find duplicate Achilles names
cld = unique([find(diff(cbi)==0) find(diff(cbi)==0)+1]);
if (numel(cld)>0)
    cldd = find([diff(cld)>1 1]); % identify sets of duplicates
    cldd = [1 cldd(1:(end-1))+1;cldd]; % find duplicate ranges
    for ci=1:size(cldd,2)
        tcld = cld(cldd(1,ci):cldd(2,ci));
        tm1(:,tcld(1)) = mean(tm1(:,tcld),2); % mean (for now)
    end
end
tm1 = tm1(:,cfo); % keep unique labels only
dmProteinAbArray = nan(size(tm1,1),mxc);
dmProteinAbArray(:,double(extractAfter(cn1u,'ACH-'))) = tm1;
clear c* t*;

%% process antibody identifiers to a metadata table, index and sort data
mtProteinAbArray = table('Size',[numel(rnC) 2],'VariableTypes',{'string','logical'}, ...
                          'VariableNames',{'antibody_name','caution_flag'});
rtC = contains(rnC,'_Caution'); rnC(rtC) = extractBefore(rnC(rtC),'_Caution');
mtProteinAbArray(:,:) = table(rnC,rtC);
[mtProteinAbArray,so] = sortrows(mtProteinAbArray,{'caution_flag','antibody_name'});
dmProteinAbArray = single(dmProteinAbArray(so,:));
mtProteinAbArray.row_idx = (1:numel(rnC))';
mtProteinAbArray = mtProteinAbArray(:,[end 1:(end-1)]);
clear r*C so;

%% write data and metadata variables to MAT file and CSV files
save build\mat\ice006MillsRPPAProteomics.mat *ProteinAbArray;
mat2csv(6,dmProteinAbArray,mtProteinAbArray,'antibody_array_proteomics',wf);
