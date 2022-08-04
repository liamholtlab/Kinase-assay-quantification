%Calculate Ribosome concentration based on Deff and absolute ribosome
%concentration measured from Cell paper:"mTORC1 controls phase separation
%and the Biophysical Properties of the cytoplasm by tunning crowding"
%Equation 12 in the supplementary. In the equation, 'D' refers to D/D0',
%'c_ribo' refers to c_ribo/c_ribo0'. Both are relative values.
%
% log(D/D0')=zeta*(phi_0/phi_m)/(1-phi_0/phi_m)*(1-c_ribo/c_ribo0')/(1-c_ribo*phi_0/phi_m)
% For simplicity let a=D/D0', b=zeta, c=phi_0/phi_m;
% Thus c_ribo/c_ribo0' = (1-loga(1-c)/(bc))/(1-loga(1-c)/b)

%% Fitting parameters from W303 hog1d CytoGEM Doolittle equation

% Require re-examining these parameters for specific strain based on
% hyperosmotic/hypo-osmotic shocks from Cell paper
b = 2.212; % b = 0.6; This is the parameters for BY4741 WT cytoGEM from Cell paper
db = 0.256; % db = 0.2; This is the parameters for BY4741 WT cytoGEM from Cell paper 
c = 0.4993; % c = 0.54; This is the parameters for BY4741 WT cytoGEM from Cell paper
dc = 0.0147; % dc = 0.5; This is the parameters for BY4741 WT cytoGEM from Cell paper
c_ribo0 = 23; %in the unit of uM

%% Calculate 'a' based on the input data of Deff and dDeff
D = input('Please input an array of Deff in um^2/s using []?\n');
dD = input ('Please input an array of its corresponding standard error in um^2/s using []?\n');
D0 = D(1); % hog1d Deff under DMSO
dD0 = dD(1); % hog1d Deff standard error under DMSO
% D0 = 0.192; % WT control Deff being 0.192um^2/s (Data from BY4741 WT cytoGEM 20211103)
% dD0 = 0.005246; % WT control Deff standard error (Data from BY4741 WT cytoGEM 20211103)
a = D/D0;
da = sqrt((dD./D).^2+(dD0/D0)^2).*a;

%% Calculate based on the Doolittle equation using fitted parameters and input Deff data
c_ribo_rela = (1-log(a)*(1-c)/(b*c))./(1-log(a)*(1-c)/b);

% Partial differentiation was calculated using wolfram alpha
dcribo_rela_da = -b*(c-1)^2./(a.*c.*((c-1).*log(a)+b).^2);
dcribo_rela_db = (c-1)^2*log(a)./(c*(b+(c-1)*log(a)).^2);
dcribo_rela_dc = -(c-1)*log(a).*((c-1)*log(a)+b*(c+1))./(c^2*((c-1)*log(a)+b).^2);
dc_ribo_rela = sqrt((dcribo_rela_da.*da).^2+(dcribo_rela_db.*db).^2+(dcribo_rela_dc.*dc).^2);

fprintf('The ribosome concentration in this condition is %fuM\n',c_ribo_rela*c_ribo0)
fprintf('with standard deviation of %fuM\n',dc_ribo_rela*c_ribo0)