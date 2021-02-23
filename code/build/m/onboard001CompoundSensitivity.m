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

%% CTRPv1 small-molecule sensitivity data
% get unique cpd names and index into data
CM1 = chkfile('v10.D3.area_under_conc_curve.txt','CTRPv1 sensitivity','\t','string',true);
[rnC1,~,fui1] = unique(CM1.cpd_name);
% cross-check unique cpd names and sort metadata
mtCTRPv1Sensitivity = chkfile('v10.M1.informer_set.txt','CTRPv1 compound info','\t','string',true);
[trnC1,~,fum1] = unique(mtCTRPv1Sensitivity.cpd_name);
assert(isequal(trnC1,rnC1),'Problem with CTRPv1 compound name matching.');
mtCTRPv1Sensitivity = mtCTRPv1Sensitivity(fum1,:); clear fum* *rnC*;
% get unique ccl names, convert to ACH names, index into data
[ct1,~,cidx] = unique(CM1.ccl_name); ctx1 = arxcellsyn(ct1);
CM1.ach_name = ctx1(cidx); [cn1,~,fuj1] = unique(CM1.ach_name);
ncn1 = double(strrep(cn1,'ACH-','')); clear c*;
assert(isequal(ncn1,unique(ncn1)),'Problem with CTRPv1 cell-line sort.');
% create sparse matrix and cast into ACH namespace
tm1 = spavg(CM1.area_under_curve,fui1,fuj1); % NORMALIZE AUCs?
pmCTRPv1SensitivAUC = castcol(tm1,ncn1,mxc);
% sort data and metadata, index metadata
[mtCTRPv1Sensitivity,so1] = sortrows(mtCTRPv1Sensitivity,{'cpd_name'});
mtCTRPv1Sensitivity = tabridx(mtCTRPv1Sensitivity);
pmCTRPv1SensitivAUC = pmCTRPv1SensitivAUC(so1,:);
clear CM* f* ncn* so* t*;

%% CTRPv2 small-molecule sensitivity data
% load relative sensitivity data
tCM2 = chkfile('v20.data.curves_post_qc.txt','CTRPv2 relative sens','\t','string',true);
tCM2 = tCM2(:,{'experiment_id','master_cpd_id','area_under_curve'});
tCM2.Properties.VariableNames(3) = {'relative_auc'};
assert(isequal(tCM2{:,[1 2]},unique(tCM2{:,[1 2]},'rows')),'Extra rows in relative AUC file.');
% load absolute sensitivity data and re-capture relative AUC values
CM2 = chkfile('new-abs-auc-with-qc.txt','CTRPv2 sensitivity','\t','string',true);
CM2.relative_auc = nan(size(CM2,1),1);
for ti=1:size(CM2,1)
    try
        CM2.relative_auc(ti) = ...
            tCM2.relative_auc(tCM2.experiment_id==CM2.experiment_id(ti) & ...
                              tCM2.master_cpd_id==CM2.master_cpd_id(ti));
    catch
        continue;
    end
end
clear t*;
% get unique cpd names and index into data
[rnC2,~,fui2] = unique(CM2.master_cpd_id);
% cross-check unique cpd names and sort metadata
mtCTRPv2Sensitivity = chkfile('v20.meta.per_compound.txt','CTRPv2 compound info','\t','string',true);
[trnC2,~,fum2] = unique(mtCTRPv2Sensitivity.master_cpd_id);
assert(isequal(trnC2,rnC2),'Problem with CTRPv2 compound ID matching.');
mtCTRPv2Sensitivity = mtCTRPv2Sensitivity(fum2,:); clear fum* *rnC*;
% get unique ccl names, convert to ACH names, index into data
ctn = chkfile('v20.meta.per_cell_line.txt','CTRPv2 cell line info','\t','string',true);
[ct2,~,cidx] = unique(CM2.master_ccl_id); ctn = ctn(ismember(ctn.master_ccl_id,ct2),:);
assert(isequal(ct2,ctn.master_ccl_id),'Problem with CTRPv2 cell-line map.');
ctn.ach_name = arxcellsyn(ctn.ccl_name); creg = contains(ctn.ach_name,'No ');
CM2.ach_name = ctn.ach_name(cidx); ctr = ctn(creg,:); % NEED TO REGISTER ACH
% % writetable(ctr,[wf '\hold\data\ctrp\ctrp-v2-register-ach.csv']);
ctn = ctn(~creg,:); cki = ismember(CM2.ach_name,ctn.ach_name);
fui2 = fui2(cki); CM2 = CM2(cki,:); [cn2,~,fuj2] = unique(CM2.ach_name);
ncn2 = double(strrep(cn2,'ACH-','')); clear c*;
assert(isequal(ncn2,unique(ncn2)),'Problem with CTRPv2 cell-line sort.');
% create sparse matrices and cast into ACH namespace
tm2a = spavg(CM2.area_under_curve,fui2,fuj2);
tm2b = spavg(CM2.ec50_log2_umol,fui2,fuj2);
tm2c = spavg(100-100*CM2.pred_pv_high_conc,fui2,fuj2);
tm2d = spavg(CM2.relative_auc,fui2,fuj2);
pmCTRPv2SensitivAUCabs = castcol(tm2a,ncn2,mxc);
pmCTRPv2SensitivL2EC50 = castcol(tm2b,ncn2,mxc);
pmCTRPv2SensitivMaxKil = castcol(tm2c,ncn2,mxc);
pmCTRPv2SensitivAUCrel = castcol(tm2d,ncn2,mxc);
% sort data and metadata, index metadata
[mtCTRPv2Sensitivity,so2] = sortrows(mtCTRPv2Sensitivity,{'cpd_name'});
mtCTRPv2Sensitivity = tabridx(mtCTRPv2Sensitivity);
pmCTRPv2SensitivAUCabs = pmCTRPv2SensitivAUCabs(so2,:);
pmCTRPv2SensitivL2EC50 = pmCTRPv2SensitivL2EC50(so2,:);
pmCTRPv2SensitivMaxKil = pmCTRPv2SensitivMaxKil(so2,:);
pmCTRPv2SensitivAUCrel = pmCTRPv2SensitivAUCrel(so2,:);
clear CM* f* ncn* so* t*;

%% GDSC (v1+v2) small-molecule sensitivity data
% get unique cpd names and index into data
CM3 = chkfile('sanger-dose-response.csv','GDSC sensitivity',',','string',true);
% CM4 = CM3(CM3.DATASET=="GDSC2",:); % HOLD FOR FURTHER DATASET REVIEW
CM3 = CM3(CM3.DATASET=="GDSC1"&contains(CM3.ARXSPAN_ID,'ACH-'),:);
[mtGDSCv1Sensitivity,~,fui3] = unique(CM3(:,{'DRUG_NAME','DRUG_ID','MAX_CONC','BROAD_ID'}),'rows');
% get unique ccl names, convert to ACH names, index into data
[ct3,~,cidx] = unique(CM3.ARXSPAN_ID); ctx3 = arxcellsyn(ct3);
CM3.ach_name = ctx3(cidx); [cn3,~,fuj3] = unique(CM3.ach_name);
ncn3 = double(strrep(cn3,'ACH-','')); clear c*;
assert(isequal(ncn3,unique(ncn3)),'Problem with GDSCv1 cell-line sort.');
% create sparse matrices and cast into ACH namespace
tm3a = spavg(CM3.Z_SCORE_PUBLISHED,fui3,fuj3);
tm3b = spavg(CM3.IC50_PUBLISHED,fui3,fuj3);
tm3c = spavg(CM3.AUC_PUBLISHED,fui3,fuj3);
pmGDSCv1SensitivZScor = castcol(tm3a,ncn3,mxc);
pmGDSCv1SensitivIC50 = castcol(tm3b,ncn3,mxc);
pmGDSCv1SensitivAUC = castcol(tm3c,ncn3,mxc);
% sort data and metadata, index metadata
[mtGDSCv1Sensitivity,so3] = sortrows(mtGDSCv1Sensitivity,{'DRUG_NAME','DRUG_ID'});
mtGDSCv1Sensitivity = tabridx(mtGDSCv1Sensitivity);
pmGDSCv1SensitivZScor = pmGDSCv1SensitivZScor(so3,:);
pmGDSCv1SensitivIC50 = pmGDSCv1SensitivIC50(so3,:);
pmGDSCv1SensitivAUC = pmGDSCv1SensitivAUC(so3,:);
clear CM* f* ncn* so* t*;
 
%% Repo-PRISM small-molecule sensitivity data
CM5 = chkfile('primary-screen-replicate-collapsed-logfold-change.csv','Repo-PRISM sensitivity',',','string',true);
rnC5 = fields(CM5); rnC5 = string(rnC5(2:end-3)); rnC5 = split(rnC5,'__');
mtRePRISMSensitivity = table(rnC5(:,1),rnC5(:,2),rnC5(:,3),'VariableNames',{'broad_id','conc_umol','assay_type'});
mtRePRISMSensitivity.broad_id = strrep(mtRePRISMSensitivity.broad_id,'_','-');
mtRePRISMSensitivity.broad_cpd_id = extractBefore(mtRePRISMSensitivity.broad_id,14);
mtRePRISMSensitivity.conc_umol = round(100*double(strrep(mtRePRISMSensitivity.conc_umol,'_','.')))/100;
[mtRePRISMSensitivity,so] = sortrows(mtRePRISMSensitivity,{'broad_cpd_id','broad_id'}); clear *rnC*;
tm5 = (double(strrep(CM5{:,2:end},'NA','NaN')))'; tm5 = tm5(so,:); % transpose and sort
cn5 = CM5{:,CM5.Properties.VariableNames{1,1}}; cn5 = arxcellsyn(cn5); % cell lines are rows
ctm = contains(cn5,'No '); cn5 = cn5(~ctm); tm5 = tm5(:,~ctm); % remove FAILED_STR cell lines
ncn5 = double(strrep(cn5,'ACH-','')); clear c* so;
assert(isequal(ncn5,unique(ncn5)),'Problem with Repo-PRISM cell-line sort.');
% cast sparse matrix into ACH namespace and index metadata
pmRePRISMSensitivL2EC50 = castcol(tm5,ncn5,mxc);
mtRePRISMSensitivity = tabridx(mtRePRISMSensitivity);
clear CM* ncn* t*;

%% write data and metadata variables to MAT file and CSV files
save build\mat\ice001CompoundSensitivity.mat *Sensitiv*;
mat2csv(1,pmCTRPv1SensitivAUC,mtCTRPv1Sensitivity,'ctrpv1_cpd_sens_auc',wf);
mat2csv(1,pmCTRPv2SensitivAUCabs,mtCTRPv2Sensitivity,'ctrpv2_cpd_sens_auc_abs',wf);
mat2csv(1,pmCTRPv2SensitivL2EC50,mtCTRPv2Sensitivity,'ctrpv2_cpd_sens_l2ec50',wf);
mat2csv(1,pmCTRPv2SensitivMaxKil,mtCTRPv2Sensitivity,'ctrpv2_cpd_sens_maxkil',wf);
mat2csv(1,pmCTRPv2SensitivAUCrel,mtCTRPv2Sensitivity,'ctrpv2_cpd_sens_auc_rel',wf);
mat2csv(1,pmGDSCv1SensitivZScor,mtGDSCv1Sensitivity,'gdscv1_cpd_sens_zscor',wf);
mat2csv(1,pmGDSCv1SensitivIC50,mtGDSCv1Sensitivity,'gdscv1_cpd_sens_ic50',wf);
mat2csv(1,pmGDSCv1SensitivAUC,mtGDSCv1Sensitivity,'gdscv1_cpd_sens_auc',wf);
mat2csv(1,pmRePRISMSensitivL2EC50,mtRePRISMSensitivity,'reprism_cpd_sens_l2ec50',wf);
