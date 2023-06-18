%% Codes for the analyses in the paper 
% "The causal effects of screen uses versus reading on the brain development in early adolescents"
%% The study mainly consists of three parts. 
%  The first part is the linear mixed-effect model, which is performed in
%  the matlab with fitlme funciton.
%% 01 - Robust association between exposures and outcomes
% we prepare a dataset with 100 samples as a demonstration, the full data
% could be downloaded through ABCD study.
clear
load('sample100.mat');
% confounders
dat_age = data.interview_age;
dat_sex = nominal(data.sex);
dat_eth = nominal(data.ethn);
dat_rac = nominal(data.race);
dat_inc = data.income;
dat_edu = data.edu;
dat_mri = nominal(data.mri_info);
dat_sit = nominal(data.site_id_l);
% exposure
dat_x  = data.reading; 
% outcome
dat_y  = data.nihtbx_picvocab_agecorrected;
% lme
dat_tab = table(dat_age,dat_sex,dat_eth,dat_rac,dat_inc,dat_edu,...
    dat_mri,dat_sit,dat_x,dat_y,...
    'VariableNames',{'age','sex','eth','rac','inc','edu','mri','sit','exp','outc'});
lme = fitlme(dat_tab,'outc~age+sex+eth+rac+inc+edu+exp+(1|sit)');
% the above codes showed a example between reading time and the performance
% of picture vovabulary, the other model is similar to the above one except
% the change of exposure and outcome variables. And the site-effect would
% be changed to the mri when we use MRI-related measures as outcome.
%%  The second part is the Mendelian randomization.
% 02 - The GWAS analyses
% In order to do MR analysis, we first need to perform GWAS analyses with
% the phenotype of interst
% the codes for the preprocessing of genetic data could be found at
% https://github.com/vwarrier/ABCD_geneticQC
% the main analysis is conducted with the GCTA software
% (https://yanglab.westlake.edu.cn/software/gcta

%% %%%%%%%%%%%%%%%%%%%%% GRM %%%%%%%%%%%%%%%%%%%%%%%%%%
cmd = ['gcta64 --bfile ',genefile,' --autosome --make-grm, --out ',file_grm];
system(cmd)
% genefile = the genetic data after preprocessing
% file_grm = the output grm file for the following fastGWA
%% %%%%%%%%%%%%%%%%%%% fastGWA %%%%%%%%%%%%%%%%%%%%%%%
% make a sparse GRM from the full GRM
cmd = ['gcta64 --grm ',file_grm,' --make-bK-sparse 0.05 --out ',file_spa];
system(cmd);
% fast GWAS
cmd = ['gcta64 --bfile ',genefile,' --grm-sparse ',file_spa,...
    ' --fastGWA-mlm --pheno ',file_phe,' --qcovar ',file_qcv,...
    ' --covar ',file_sex,' --out ',file_gwa];
system(cmd);
% file_phe = phenotype
% file_qcv = quantitative covariates (e.g. age, top 20 PCs)
% file_sex = sex
% file_gwa = the results
%% 03 one sample MR analysis
% This is performed with R. see OSMR.R
%% 04 two-step MR analysis
% This is similiar to the above MR analysis except that we first calculate
% the MR effect between screen use and reading time, and also calcualte 
% the MR effect between reading time and other outcomes.

% Asumme a(1) = the effect size of first MR, 
%        a(2) = the se of the first MR.
%        b(1) = the effect size of second MR, 
%        b(2) = the se of the second MR.
% The mediation effect size equal to the product of two effect size.
m  = a*b;
% And the se:
se = sqrt((a(1)*b(2))^2 + (b(1)*a(2))^2 +a(2)^2*b(2)^2);
