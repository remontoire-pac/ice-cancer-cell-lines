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

%% load CCLE metabolomics data (standard-based MS) and extract headers
DM1 = chkfile('CCLE_metabolomics_20190502.csv','Gygi proteomics',',','string',true);
rnC = DM1.Properties.VariableNames';
trnC = DM1.Properties.VariableDescriptions';
rnC(contains(trnC,'Original column heading: ')) = trnC(contains(trnC,'Original column heading: '));
rnC = strrep(rnC,"Original column heading: '","");
rnC(contains(trnC,'Original column heading: ')) = regexprep(rnC(contains(trnC,'Original column heading: ')),'''$',"");
rnC = string(rnC(3:end));
cn1 = DM1{:,DM1.Properties.VariableNames{1,2}}; cn1 = arxcellsyn(cn1); % cell lines are rows
ctm = ~(contains(cn1,'No CELL_SAMPLE_ID found')); % one not found initially
tm1 = (double(strrep(DM1{:,3:end},'NA','NaN')))'; % transpose
cn1 = cn1(ctm); tm1 = tm1(:,ctm);
clear ctm DM1 trnC;

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
dmMetaboliteMassSpec = nan(size(tm1,1),mxc);
dmMetaboliteMassSpec(:,double(extractAfter(cn1u,'ACH-'))) = tm1;
clear c* t*;

%% process metabolite identifiers to a metadata table, index and sort data
mtMetaboliteMassSpec = table('Size',[numel(rnC) 1],'VariableTypes',{'string'}, ...
                         'VariableNames',{'metabolite_name'});
mtMetaboliteMassSpec(:,:) = table(rnC);
[mtMetaboliteMassSpec,so] = sortrows(mtMetaboliteMassSpec,{'metabolite_name'});
dmMetaboliteMassSpec = single(dmMetaboliteMassSpec(so,:));
mtMetaboliteMassSpec.row_idx = (1:numel(rnC))';
mtMetaboliteMassSpec = mtMetaboliteMassSpec(:,[end 1:(end-1)]);
clear rnC so;

%% write data and metadata variables to MAT file and CSV files
save build\mat\ice008CCLEStdMetabolomics.mat *MetaboliteMassSpec;
mat2csv(8,dmMetaboliteMassSpec,mtMetaboliteMassSpec,'std_ms_metabolomics',wf);
