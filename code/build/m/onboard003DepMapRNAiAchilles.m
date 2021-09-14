%% setup workspace
if (~exist('wf','var')) 
    wf = regexp(matlab.desktop.editor.getActiveFilename,filesep,'split');
    wf = strjoin(wf(1:(numel(wf)-4)),filesep); % ICE root folder
end
iceopts(wf,false,false);        % warnings on (T/F), 'use' mode on (T/F)
clearvars -except a* wf;        % adapt as needed and preferred

%% find range of cell-line namespace from cell-line annotations file
assert(exist('ice000CellLineAnnotation.mat','file')==2, ...
             'Cell-line annotation MAT file required to fix column number.');
load ice000CellLineAnnotation.mat db*;
mxc = size(dbCellLineAnno,2);
clear db*;

%% load Achilles dependency data (RNAi) and extract headers
DM1 = chkfile('D2_combined_gene_dep_scores.csv','RNAi dependency',',','string',true);
rnC = DM1{:,DM1.Properties.VariableNames{1,1}};
cn1 = fields(DM1); cn1 = string(cn1(2:(end-3)));
cn1 = erase(cn1,'x'); cn1 = arxcellsyn(cn1);
ctm = numel(find(contains(cn1,'No CELL_SAMPLE_ID found')));
assert(ctm==0,'Achilles name not produced for %i input names.',ctm);
tm1 = (double(strrep(DM1{:,2:end},'NA','NaN')));
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
dmRNAiAchilles = nan(size(tm1,1),mxc);
dmRNAiAchilles(:,double(extractAfter(cn1u,'ACH-'))) = tm1;
clear c* t*;

%% process gene identifiers to a metadata table, index and sort data
mtRNAiAchilles = table('Size',[numel(rnC) 2],'VariableTypes',{'string','string'}, ...
                         'VariableNames',{'gene_symbol','entrez_gene_id'});
mtRNAiAchilles(:,:) = table(extractBefore(rnC,' ('),extractBefore(extractAfter(rnC,' ('),')'));
[mtRNAiAchilles,so] = sortrows(mtRNAiAchilles,{'gene_symbol','entrez_gene_id'});
dmRNAiAchilles = single(dmRNAiAchilles(so,:));
mtRNAiAchilles.row_idx = (1:numel(rnC))';
mtRNAiAchilles = mtRNAiAchilles(:,[end 1:(end-1)]);
assert(isequal(count(mtRNAiAchilles.gene_symbol,'&'), ...
               count(mtRNAiAchilles.entrez_gene_id,'&')), ...
               'Mismatched array sizes between symbols and IDs.');
clear rnC so;

%% write data and metadata variables to MAT file and CSV files
save(['build' filesep 'mat' filesep 'ice003DepMapRNAiAchilles.mat'],'*RNAiAchilles');
mat2csv(3,dmRNAiAchilles,mtRNAiAchilles,'rnai_gene_dependency',wf);
