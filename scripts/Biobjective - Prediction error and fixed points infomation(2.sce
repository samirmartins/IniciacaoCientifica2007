// Biobjective otimization of the model's parameters for prediction of output voltage
// of the Buck converter CC-CC using as objectives the fixed points information 
// and the prediction error 
// 
// Programmer: Samir Angelo Milani Martins
// Electrical Engineering, UFSJ, Brazil
// Fifth period - LET'S WITH GOD!
//__________________________________________________________________________________

clear

clc
tic()
// Initializing Variables
k = 1;
r = 1;
w = 1;
col_par = 1;
norm_comparison = 10^6;
max_lag = 2;
max_var = 0.01;
iteration = 0;
number_of_models = 2000;
max_iteration = 15;
teta_aux(1:9) = 10^(-6); // Initializing the while loop
// Load the datas.
loadmatfile('-ascii','buck_id.dat');    
loadmatfile('-ascii','buck_val.dat');    

// Static curve of Buck Converter
// v0 = 4/3vd - (vd/3)*u, where vd = 24

// Fixed points obtained for u = 2.2,
vpf = (4/3)*24 - (24/3)*2.2;

// Mount the equation where the fixed points of the model can be estimated.
fixed_points = poly([vpf],'ybar');

// WARNING! VALID ONLY THIS MODEL AND THIS FIXED POINTS

sigma_capital = [0;0;0;0;-14.4];

// Model's structure that will be utilized:
// y(k) = ay(k - 1) + by(k - 2) + c + du(k - 1) + eu(k - 1)^2 + fu(k - 2)u(k - 1) + gu(k - 2) + hu(k - 2)^2
// Maximum Lag = 2

// Sigma LS:

// Backup of the datas

// Identification
y = buck_id(:,3); // Datas of output
u = buck_id(:,2); // Datas of input
n = buck_id(:,1); // Datas of samples

// Validation
yv = buck_val(:,3); // Datas of output
uv = buck_val(:,2); // Datas of input
nv = buck_val(:,1); // Datas of samples

// Mount tha matrix PSI

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

// Mounting sigma LS
// y(k) = teta_ls1(1)y(k - 1) + teta_ls1(2)y(k - 2) + teta_ls1(3) + teta_ls1(4)u(k - 1) + teta_ls1(5)u(k - 1)^2 + teta_ls1(6)u(k - 2)u(k - 1) + teta_ls1(7)u(k - 2) + teta_ls1(8)u(k - 2)^2
sigma_ls = [0 (teta_ls1(1)+teta_ls1(2)-1) (teta_ls1(5)+teta_ls1(6)+teta_ls1(8)) (teta_ls1(4)+teta_ls1(7)) teta_ls1(3)]';

// Obtaining minuscule sigma
minuscule_sigma = (sqrt(norm(sigma_ls)))/(sqrt(norm(sigma_capital)))*sigma_capital;

sigma = minuscule_sigma;
sigma(2) = sigma(2)+1;

// Matrix S that will map the parameter 
S = [0 0 0 0 0 0 0 0 0;1 1 0 0 0 0 0 0 0; 0 0 0 0 1 1 0 1 0; 0 0 0 1 0 0 1 0 0; 0 0 1 0 0 0 0 0 0]; 
s = [0 0 0 0 0 0 0 0;1 1 0 0 0 0 0 0; 0 0 0 0 1 1 0 1; 0 0 0 1 0 0 1 0; 0 0 1 0 0 0 0 0];

// Extended Least Squares

for w1 = logspace(-6,0,number_of_models)// w1 = weigth of prediction's error and w2 = weigth of fixed points information
w2 = 1-w1;

// Get the parameters
teta_ls(:,col_par) = (pinv(w1*PSI'*PSI + w2*s'*s))*(w1*PSI'*Y + w2*s'*sigma);

// Residues vector
resid = Y - PSI*teta_ls(:,col_par);
PSI_EX_2 = [PSI resid];

// Pareto-optimum set
pareto(:,col_par) = (pinv(w1*PSI_EX_2'*PSI_EX_2 + w2*S'*S))*(w1*PSI_EX_2'*Y + w2*S'*sigma);
iteration = 0;
  while abs((teta_aux(1) - pareto(1,col_par))/teta_aux(1)) > max_var & abs((teta_aux(2) - pareto(2,col_par))/teta_aux(2)) > max_var & abs((teta_aux(3) - pareto(3,col_par))/teta_aux(3)) > max_var & abs((teta_aux(4) - pareto(4,col_par))/teta_aux(4)) > max_var & abs((teta_aux(5) - pareto(5,col_par))/teta_aux(5)) > max_var & abs((teta_aux(6) - pareto(6,col_par))/teta_aux(6)) > max_var & abs((teta_aux(7) - pareto(7,col_par))/teta_aux(7)) > max_var & abs((teta_aux(8) - pareto(8,col_par))/teta_aux(8)) > max_var
    iteration = iteration + 1;
    teta_aux = pareto(:,col_par); 
    resid = Y - PSI_EX_2*pareto(:,col_par);
    PSI_EX_2 = [PSI resid];
    pareto(:,col_par) = (pinv(w1*PSI_EX_2'*PSI_EX_2 + w2*S'*S))*(w1*PSI_EX_2'*Y + w2*S'*sigma);
    if iteration > max_iteration then
      break
    end
  end
  col_par = col_par+1;
  if col_par ==992; 
    w1_used = w1;
    w2_used = w2;
  end
end  

// Norm of Cost functions
// cost_ls(k) = (Y - PSI*pareto(1:8,k))'*(Y-PSI*pareto(1:8,k))
// cost_fp(k) = (sigma - S_s*pareto(1:8,k))'*(sigma - S_s*pareto(1:8,k))
for k = (1:number_of_models) 
cost_ls(k) = (Y - PSI*pareto(1:8,k))'*(Y-PSI*pareto(1:8,k));
cost_fp(k) = (sigma - s*pareto(1:8,k))'*(sigma - s*pareto(1:8,k));
norm_cost(k) = sqrt(cost_ls(k)^2 + cost_fp(k)^2);
end

// Positions that must be retired of the cost vector
norm2 = [cost_ls cost_fp];
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
scf(1)
plot(enfim(:,2),enfim(:,1),'.');
xtitle('Costs Functions','Fixed points information','Prediction Error');

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

// Calculus of the RMSE index

rmse = sqrt(sum((buck_val(2:$-1,3)-yforecastparetofree_val(3:$)).^2))/sqrt(sum((buck_val(2:$-1,:)-mean(buck_val(2:$-1,:))).^2));

//__________________________________________________________________________________________________
// Ploting the simulated datas of the pareto set with real datas

scf(2)
plot(1:998,yforecastpareto1(4:1001),'black');
plot(buck_id(3:(length(buck_id)/3),3),'red');
xtitle('One-step forecast - red(real datas), black(paretos datas) - Identification','Number of samples','Voltage');

scf(3)
plot(1:998,yforecastparetofree(4:1001),'black');
plot(buck_id(3:(length(buck_id)/3),3),'red');
xtitle('Free forecast - red(real datas), black(paretos datas) - Identification','Number of samples','Voltage');

scf(4)
plot(1:996,yforecastpareto1_val(4:999),'black');
plot(buck_val(3:(length(buck_val)/3),3),'red');
xtitle('One-step forecast - red(real datas), black(paretos datas) - Validation','Number of samples','Voltage');

scf(5)
plot(1:996,yforecastparetofree_val(4:999),'black');
plot(buck_val(3:(length(buck_val)/3),3),'red');
xtitle('Free forecast - red(real datas), black(paretos datas) - Validation','Number of samples','Voltage');

time = toc();

