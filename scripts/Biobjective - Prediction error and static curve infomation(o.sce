// Biobjective otimization of the model's parameters for prediction of output voltage
// of the Buck converter CC-CC using as objectives the static curve information 
// and the prediction error 
// 
// Programmer: Samir Angelo Milani Martins
// Electrical Engineering, UFSJ, Brazil
// Fifth period
// April - 2008 - LET'S WITH GOD!
//__________________________________________________________________________________

// Clearing the window and the variables
clear; clc;

// Counting the time
tic()

// Initializing variables
cont = 1; // Position of fixed points vector
vd = 24; // Voltage of supply of the buck converter
nfp = 50; // Number of fixed points that will be incorporated
k = 1; 
r = 1;
w = 1;
col_par = 1;
norm_comparison = 10^6;
max_var = 0.01;
iteration = 0;
number_of_models = 2000;  // Number of models of the pareto
max_iteration = 15; // Maximum iteration of the ELS
teta_aux(1:9) = 10^(-6); // Initializing the while loop of the ELS

// Load the datas.
loadmatfile('-ascii','buck_id.dat');    
loadmatfile('-ascii','buck_val.dat');    

// Managing data
len_data_id = size(buck_id,1)/3;
c1_id = buck_id(1:len_data_id,1);
c2_id = buck_id(len_data_id+1:2*len_data_id,1);
c3_id = buck_id(2*len_data_id+1:$,1);
buck_id = [];
buck_id = [c1_id, c2_id, c3_id];
len_data_val = size(buck_val,1)/3;
c1_val = buck_val(1:len_data_val,1);
c2_val = buck_val(len_data_val+1:2*len_data_val,1);
c3_val = buck_val(2*len_data_val+1:$,1);
buck_val = [];
buck_val = [c1_val, c2_val, c3_val];

// Backup of the datas

// Identification
y = buck_id(:,3); // Datas of output
u = buck_id(:,2); // Datas of input
n = buck_id(:,1); // Datas of samples

// Validation
yv = buck_val(:,3); // Datas of output
uv = buck_val(:,2); // Datas of input
nv = buck_val(:,1); // Datas of samples

// Static curve of Buck Converter
// v0 = 4/3vd - (vd/3)*u, where vd = 24

// Model's structure that will be utilized:
// y(k) = ay(k - 1) + by(k - 2) + c + du(k - 1) + eu(k - 1)^2 + fu(k - 2)u(k - 1) + gu(k - 2) + hu(k - 2)^2
// Maximum Lag = 2

// Fixed Points that will be utilized for the static curve incorporation 
// uf = linspace(0,4,npf)

for uf = linspace(0,4,nfp) // Define the quantity and the gap between the fixed points
  yf(cont) = (4/3)*vd - (vd/3)*uf;
  // Mounting q for each fixed point and Q of q
  // qi = [1 y u u^2]; Q = [qi], for i =1:n  
  Q(cont,:) = [1 yf(cont) uf uf^2];
  cont = cont + 1;  
end

// Matrix R that will be map the parameters

R = [0 0 1 0 0 0 0 0 0;1 1 0 0 0 0 0 0 0;0 0 0 1 0 0 1 0 0;0 0 0 0 1 1 0 1 0];  // R extended
r1 = [0 0 1 0 0 0 0 0;1 1 0 0 0 0 0 0;0 0 0 1 0 0 1 0;0 0 0 0 1 1 0 1];         // R

// Mount the matrix PSI

// Prepare for the mount of PSI matrix

//u(k-1)^2

t = u(2:(length(u)-1),:);

while k <= length(t)
e(k) = t(k)*t(k);    // u(k-1)^2
k = k+1;
end

// s = u(k-2)u(k-1)

s = u(1:(length(u)-2),:);     // u(k-2)
a = u(2:length(u)-1,:);     // u(k-1)

while r <= length(s)
d(r) = s(r)*a(r);
r = r+1;
end

//u(k-2)^2

q = u(1:(length(u)-2),:);

while w <= length(q)
p(w) = q(w)*q(w);    // u(k-2)^2
w = w+1;
end

PSI = [y(2:(length(y)-1),:),y(1:(length(y)-2),:),ones(length(e),1),u(2:(length(u)-1),:),e,d,u(1:(length(u)-2),:),p];

Y = y(3:length(y),:);
teta_ls1 = (pinv(PSI'*PSI))*(PSI'*Y);

// Extended Least Squares

for w1 = logspace(-6,0,number_of_models)// w1 = weigth of prediction's error and w2 = weigth of fixed points information
w2 = 1-w1;

// Get the parameters
teta_ls(:,col_par) = (pinv(w1*PSI'*PSI + w2*(Q*r1)'*(Q*r1)))*(w1*PSI'*Y + w2*(Q*r1)'*yf);

// Residues vector
resid = Y - PSI*teta_ls(:,col_par);
PSI_EX_2 = [PSI resid];

// Pareto-optimum set
pareto(:,col_par) = (pinv(w1*PSI_EX_2'*PSI_EX_2 + w2*(Q*R)'*(Q*R)))*(w1*PSI_EX_2'*Y + w2*(Q*R)'*yf);
iteration = 0;
  while abs((teta_aux(1) - pareto(1,col_par))/teta_aux(1)) > max_var & abs((teta_aux(2) - pareto(2,col_par))/teta_aux(2)) > max_var & abs((teta_aux(3) - pareto(3,col_par))/teta_aux(3)) > max_var & abs((teta_aux(4) - pareto(4,col_par))/teta_aux(4)) > max_var & abs((teta_aux(5) - pareto(5,col_par))/teta_aux(5)) > max_var & abs((teta_aux(6) - pareto(6,col_par))/teta_aux(6)) > max_var & abs((teta_aux(7) - pareto(7,col_par))/teta_aux(7)) > max_var & abs((teta_aux(8) - pareto(8,col_par))/teta_aux(8)) > max_var
    iteration = iteration + 1;
    teta_aux = pareto(:,col_par); 
    resid = Y - PSI_EX_2*pareto(:,col_par);
    PSI_EX_2 = [PSI resid];
    pareto(:,col_par) = (pinv(w1*PSI_EX_2'*PSI_EX_2 + w2*(Q*R)'*(Q*R)))*(w1*PSI_EX_2'*Y + w2*(Q*R)'*yf);
    if iteration > max_iteration then
      break
    end
  end
  col_par = col_par+1;
end  


// Norm of Cost functions
// cost_ls(k) = (Y - PSI*pareto(:,k))'*(Y-PSI*pareto(:,k))
// cost_sf(k) = (yf - (Q*R)*pareto(:,k))'*(yf - (Q*R)*pareto(:,k))
for k = (1:number_of_models) // 8 => number of regressors
cost_ls(k) = (Y - PSI*pareto(1:8,k))'*(Y-PSI*pareto(1:8,k));
cost_sf(k) = (yf - (Q*r1)*pareto(1:8,k))'*(yf - (Q*r1)*pareto(1:8,k));
norm_cost(k) = sqrt(cost_ls(k)^2 + cost_sf(k)^2);
end

// Positions that must be retired of the cost vector
norm2 = [cost_ls cost_sf];
coleta = 1;
e = 1>0;
for k=1:number_of_models
  s = find((norm2(k,1)-norm2(1:$,1)>0)&(norm2(k,2)-norm2(1:$,2)>0)==e);
  if length(s)>0
    posicao(coleta) = k;
    coleta = coleta +1;
  end
end

// Vector enfim contaim the plotted values
s = 1;
m = 1;
for k = 1:number_of_models
  if k>posicao(m)
    if m<length(posicao)
      m = m+1;
    else
      m = length(posicao);
    end
  end
  if k<>posicao(m)
    enfim(s,:) = norm2(k,:);
    s = s+1;
  end
end

// Comparison Between the norm of the objectives
for k = 1:length(norm_cost)
  if norm_cost(k)<norm_comparison
    norm_comparison = norm_cost(k);
    min_norm_pos = k;
  end
end

// Parameters of the pareto "min_norm_pos" minimized the norm of the cost functions
// Plot the costs functions
scf(1);
plot(enfim(:,2),enfim(:,1),'.r');
xtitle('Costs Functions','Static Curve Information','Prediction Error');

//__________________________________________________________________________________________________
//______________________________________Predictions - Identification________________________________
// Algorithm that will make the one-step-ahead prediction (Identification)

// One-step-ahead prediction
// Forecast preparation
y = buck_id(:,3); // separe the output column
yp = y; // Output's backup used in forecast
up = u; // Input's backup used in forecast
aux = 0; // Auxiliar variable
n = 1; // Vary n => Vary Forecast horizont

// Forecast
for k = 3:length(y)
   yp(k) = pareto(1,min_norm_pos)*yp(k - 1) + pareto(2,min_norm_pos)*yp(k - 2) + pareto(3,min_norm_pos) + pareto(4,min_norm_pos)*up(k - 1) + pareto(5,min_norm_pos)*up(k - 1)*up(k-1) + pareto(6,min_norm_pos)*up(k - 2)*up(k - 1) + pareto(7,min_norm_pos)*up(k - 2) + pareto(8,min_norm_pos)*up(k - 2)*up(k - 2);
   yforecastpareto1(k)=yp(k); // save the forecasts
   aux=aux+1; // help to verify if the datas can be updated
   if ((aux/n)-floor(aux/n))==0 // verify if the datas can be updated
     yp=y; // update datas
   end
end

// Algorithm that will make the free prediction (Identification)

// Free prediction
// Forecast preparation
y = buck_id(:,3); // separe the output column
yp = y; // Output backup's used in forecast
up = u; // Input's backup used in forecast
aux = 0; // Auxiliar variable
n = 10000000000; // Vary n => Vary Forecast horizont

// Forecast
for k = 3:length(y)
   yp(k) = pareto(1,min_norm_pos)*yp(k - 1) + pareto(2,min_norm_pos)*yp(k - 2) + pareto(3,min_norm_pos) + pareto(4,min_norm_pos)*up(k - 1) + pareto(5,min_norm_pos)*up(k - 1)*up(k-1) + pareto(6,min_norm_pos)*up(k - 2)*up(k - 1) + pareto(7,min_norm_pos)*up(k - 2) + pareto(8,min_norm_pos)*up(k - 2)*up(k - 2);
   yforecastparetofree(k)=yp(k); // save the forecasts
   aux=aux+1; // help to verify if the datas can be updated
   if ((aux/n)-floor(aux/n))==0 // verify if the datas can be updated
     yp=y; // update datas
   end
end

//________________________________________Model's validate__________________________________________

// Algorithm that will make the one-step-ahead prediction (Validation)

// One-step-ahead prediction
// Forecast preparation
clear y; clear yp; clear up; clear aux; clear n;
y = buck_val(:,3); // separe the output column
yp = y; // Output backup's used in forecast
up = u; // Input's backup used in forecast
aux = 0; // Auxiliar variable
n = 1; // Vary n => Vary Forecast horizont

// Forecast
for k = 3:length(y)
   yp(k) = pareto(1,min_norm_pos)*yp(k - 1) + pareto(2,min_norm_pos)*yp(k - 2) + pareto(3,min_norm_pos) + pareto(4,min_norm_pos)*up(k - 1) + pareto(5,min_norm_pos)*up(k - 1)*up(k-1) + pareto(6,min_norm_pos)*up(k - 2)*up(k - 1) + pareto(7,min_norm_pos)*up(k - 2) + pareto(8,min_norm_pos)*up(k - 2)*up(k - 2);
   yforecastpareto1_val(k)=yp(k); // save the forecasts
   aux=aux+1; // help to verify if the datas can be updated
   if ((aux/n)-floor(aux/n))==0 // verify if the datas can be updated
     yp=y; // update datas
   end
end

// Algorithm that will make the free prediction (Validation)

// Free prediction
// Forecast preparation
y = buck_val(:,3); // separe the output column
yp = y; // Output backup's used in forecast
up = u; // Input's backup used in forecast
aux = 0; // Auxiliar variable
n = 10000000000; // Vary n => Vary Forecast horizont

// Forecast
for k = 3:length(y)
   yp(k) = pareto(1,min_norm_pos)*yp(k - 1) + pareto(2,min_norm_pos)*yp(k - 2) + pareto(3,min_norm_pos) + pareto(4,min_norm_pos)*up(k - 1) + pareto(5,min_norm_pos)*up(k - 1)*up(k-1) + pareto(6,min_norm_pos)*up(k - 2)*up(k - 1) + pareto(7,min_norm_pos)*up(k - 2) + pareto(8,min_norm_pos)*up(k - 2)*up(k - 2);
   yforecastparetofree_val(k)=yp(k); // save the forecasts
   aux=aux+1; // help to verify if the datas can be updated
   if ((aux/n)-floor(aux/n))==0 // verify if the datas can be updated
     yp=y; // update datas
   end
end

// RMSE index
rmse = sqrt(sum((buck_val(2:$-1,3)-yforecastparetofree_val(3:$)).^2))/sqrt(sum((buck_val(2:$-1,:)-mean(buck_val(2:$-1,:))).^2));

//__________________________________________________________________________________________________
// Ploting the simulated datas of the pareto set with real datas

scf(2);
plot(1:998,yforecastpareto1(4:1001),'black');
plot(buck_id(3:(length(buck_id)/3),3),'red');
xtitle('One-step forecast - red(real datas), black(paretos datas) - Identification','Number of samples','Voltage');

scf(3);
plot(1:998,yforecastparetofree(4:1001),'black');
plot(buck_id(3:(length(buck_id)/3),3),'red');
xtitle('Free forecast - red(real datas), black(paretos datas) - Identification','Number of samples','Voltage');

scf(4);
plot(1:996,yforecastpareto1_val(4:999),'blue');
plot(buck_val(3:(length(buck_val)/3),3),'red');
xtitle('One-step forecast - red(real datas), black(paretos datas) - Validation','Number of samples','Voltage');

scf(5);
plot(1:996,yforecastparetofree_val(4:999),'black');
plot(buck_val(3:(length(buck_val)/3),3),'red');
xtitle('Free forecast - red(real datas), black(paretos datas) - Validation','Number of samples','Voltage');


time = toc();
