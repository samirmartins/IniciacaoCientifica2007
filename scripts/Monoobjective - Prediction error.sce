// Biobjective otimization of the model's parameters of prediction of output voltage
// of the Buck converter CC-CC using as objectives the fixed points information 
// and the prediction error 
// 
// Programmer: Samir Angelo Milani Martins
// Electrical Engineering, UFSJ, Brazil
// Fifth period - LET'S WITH GOD!
//__________________________________________________________________________________

clear

clc

// Initializing Variables
k = 1;
r = 1;
w = 1;
radician_n = 0;
radician_d = 0;
teta_aux(1:8)=10^(-6);
iteration = 0;
max_var = 0.01; // Maximum variance of the parameters
max_iteration = 15;

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

// Model that will be utilized:
// y(k) = ay(k - 1) + by(k - 2) + c + du(k - 1) + eu(k - 1)^2 + fu(k - 2)u(k - 1) + gu(k - 2) + hu(k - 2)^2
// Maximum Lag = 2

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

// Extended Least Squares

PSI = [y(2:(length(y)-1),:),y(1:(length(y)-2),:),ones(999,1),u(2:(length(u)-1),:),e,d,u(1:(length(u)-2),:),p];

Y = y(3:length(y),:);

teta_ls = (pinv(PSI'*PSI))*(PSI'*Y);

resid = Y - PSI*teta_ls;

PSI_EX_2 = [PSI resid];
teta_els = (pinv(PSI_EX_2'*PSI_EX_2)*(PSI_EX_2'*Y));

while abs((teta_aux(1) - teta_els(1))/teta_aux(1)) > max_var & abs((teta_aux(2) - teta_els(2))/teta_aux(2)) > max_var & abs((teta_aux(3) - teta_els(3))/teta_aux(3)) > max_var & abs((teta_aux(4) - teta_els(4))/teta_aux(4)) > max_var & abs((teta_aux(5) - teta_els(5))/teta_aux(5)) > max_var & abs((teta_aux(6) - teta_els(6))/teta_aux(6)) > max_var & abs((teta_aux(7) - teta_els(7))/teta_aux(7)) > max_var & abs((teta_aux(8) - teta_els(8))/teta_aux(8)) > max_var

iteration = iteration + 1;  // Number of iteration
if iteration > max_iteration then 
  break
end

teta_aux = teta_els; 
resid = Y - PSI_EX_2*teta_els;
PSI_EX_2 = [PSI resid];
teta_els = (pinv(PSI_EX_2'*PSI_EX_2)*(PSI_EX_2'*Y));

end


// Measurement of the cost function

cost_ls = (Y - PSI*teta_els(1:length(teta_els)-1))'*(Y - PSI*teta_els(1:length(teta_els)-1));
//__________________________________________________________________________________________________
//______________________________________________Predictions_________________________________________
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
   yp(k) = teta_els(1)*yp(k - 1) + teta_els(2)*yp(k - 2) + teta_els(3) + teta_els(4)*up(k - 1) + teta_els(5)*up(k - 1)*up(k-1) + teta_els(6)*up(k - 2)*up(k - 1) + teta_els(7)*up(k - 2) + teta_els(8)*up(k - 2)*up(k - 2);
   yforecast1(k)=yp(k); // save the forecasts
   aux=aux+1; // help to verify if the datas can be updated
   if ((aux/n)-floor(aux/n))==0 // verify if the datas can be updated
     yp=y; // update datas
   end
end
//__________________________________________________________________________________________________
// Algorithm that will make the six-step-ahead prediction (Identification)

// Six-step-ahead prediction
// Forecast preparation
y = buck_id(:,3); // separe the output column
yp = y; // Output backup's used in forecast
up = u; // Input's backup used in forecast
aux = 0; // Auxiliar variable
n = 6; // Vary n => Vary Forecast horizont

// Forecast
for k = 3:length(y)
   yp(k) = teta_els(1)*yp(k - 1) + teta_els(2)*yp(k - 2) + teta_els(3) + teta_els(4)*up(k - 1) + teta_els(5)*up(k - 1)*up(k-1) + teta_els(6)*up(k - 2)*up(k - 1) + teta_els(7)*up(k - 2) + teta_els(8)*up(k - 2)*up(k - 2);
   yforecast6(k)=yp(k); // save the forecasts
   aux=aux+1; // help to verify if the datas can be updated
   if ((aux/n)-floor(aux/n))==0 // verify if the datas can be updated
     yp=y; // update datas
   end
end
//__________________________________________________________________________________________________
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
   yp(k) = teta_els(1)*yp(k - 1) + teta_els(2)*yp(k - 2) + teta_els(3) + teta_els(4)*up(k - 1) + teta_els(5)*up(k - 1)*up(k-1) + teta_els(6)*up(k - 2)*up(k - 1) + teta_els(7)*up(k - 2) + teta_els(8)*up(k - 2)*up(k - 2);
   yforecastfree(k)=yp(k); // save the forecasts
   aux=aux+1; // help to verify if the datas can be updated
   if ((aux/n)-floor(aux/n))==0 // verify if the datas can be updated
     yp=y; // update datas
   end
end

//__________________________________________________________________________________________________
//________________________________________Model's validate__________________________________________

// Algorithm that will make the one-step-ahead prediction (Validation)

// One-step-ahead prediction
// Forecast preparation
y = buck_val(:,3); // separe the output column
yp = y; // Output backup's used in forecast
up = u; // Input's backup used in forecast
aux = 0; // Auxiliar variable
n = 1; // Vary n => Vary Forecast horizont

// Forecast
for k = 3:length(y)
   yp(k) = teta_els(1)*yp(k - 1) + teta_els(2)*yp(k - 2) + teta_els(3) + teta_els(4)*up(k - 1) + teta_els(5)*up(k - 1)*up(k-1) + teta_els(6)*up(k - 2)*up(k - 1) + teta_els(7)*up(k - 2) + teta_els(8)*up(k - 2)*up(k - 2);
   yforecast1_val(k)=yp(k); // save the forecasts
   aux=aux+1; // help to verify if the datas can be updated
   if ((aux/n)-floor(aux/n))==0 // verify if the datas can be updated
     yp=y; // update datas
   end
end
//__________________________________________________________________________________________________
// Algorithm that will make the six-step-ahead prediction (Validation)

// Six-step-ahead prediction
// Forecast preparation
y = buck_val(:,3); // separe the output column
yp = y; // Output backup's used in forecast
up = u; // Input's backup used in forecast
aux = 0; // Auxiliar variable
n = 6; // Vary n => Vary Forecast horizont

// Forecast
for k = 3:length(y)
   yp(k) = teta_els(1)*yp(k - 1) + teta_els(2)*yp(k - 2) + teta_els(3) + teta_els(4)*up(k - 1) + teta_els(5)*up(k - 1)*up(k-1) + teta_els(6)*up(k - 2)*up(k - 1) + teta_els(7)*up(k - 2) + teta_els(8)*up(k - 2)*up(k - 2);
   yforecast6_val(k)=yp(k); // save the forecasts
   aux=aux+1; // help to verify if the datas can be updated
   if ((aux/n)-floor(aux/n))==0 // verify if the datas can be updated
     yp=y; // update datas
   end
end
//__________________________________________________________________________________________________
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
   yp(k) = teta_els(1)*yp(k - 1) + teta_els(2)*yp(k - 2) + teta_els(3) + teta_els(4)*up(k - 1) + teta_els(5)*up(k - 1)*up(k-1) + teta_els(6)*up(k - 2)*up(k - 1) + teta_els(7)*up(k - 2) + teta_els(8)*up(k - 2)*up(k - 2);
   yforecastfree_val(k)=yp(k); // save the forecasts
   aux=aux+1; // help to verify if the datas can be updated
   if ((aux/n)-floor(aux/n))==0 // verify if the datas can be updated
     yp=y; // update datas
   end
end

//________________________________Calculus of the RMSE index________________________________________

rmse = sqrt(sum((buck_val(2:$-1,3)-yforecastfree_val(3:$)).^2))/sqrt(sum((buck_val(2:$-1,:)-mean(buck_val(2:$-1,:))).^2));

//__________________________________________________________________________________________________
//______________________________________________Graphics plot_______________________________________

// Real datas

//Output
subplot(211)
plot(buck_id(:,1),buck_id(:,3),'red');
subplot(211)
plot(buck_val(:,1),buck_val(:,3),'blue');
//xtitle('Buck converter - Output datas, Continuous line => Identification; Delineated line => Validation','SAMPLES','VOLTAGE');

//Input
subplot(212)
plot(buck_id(:,1),buck_id(:,2),'red');
subplot(212)
plot(buck_val(:,1),buck_val(:,2),'blue');
//xtitle('Buck converter - Input datas, Continuous line => Identification; Delineated line => Validation','SAMPLES','VOLTAGE');

//_________________________________________Model Identification____________________________________________

// One-step-ahead
scf(2)
plot(buck_id(:,1),buck_id(:,3),'red');
plot(buck_id(2:(length(buck_id)/3)-1,1),yforecast1(3:length(yforecast1)),'blue');
//xtitle('One-step-ahead - Identification: Continuous line => Real datas; Delineated line => Model monoobjective','SAMPLES','VOLTAGE');

// Six-step-ahead
scf(3)
plot(buck_id(:,1),buck_id(:,3),'red');
plot(buck_id(2:(length(buck_id)/3)-1,1),yforecast6(3:length(yforecast6)),'blue');
//xtitle('Six-step-ahead - Identification: Continuous line => Real datas; Delineated line => Model monoobjective','SAMPLES','VOLTAGE');

// Free prediction
scf(4)
plot(buck_id(:,1),buck_id(:,3),'red');
plot(buck_id(2:(length(buck_id)/3)-1,1),yforecastfree(3:length(yforecastfree)),'blue');
//xtitle('Free prediction - Identification: Continuous line => Real datas; Delineated line => Model monoobjective','SAMPLES','VOLTAGE');

//___________________________________________Model Validation______________________________________________
// Fault the lenght
// One-step-ahead
scf(5)
plot(buck_val(:,1),buck_val(:,3),'red');
plot(buck_val(2:(length(buck_val)/3)-1,1),yforecast1_val(3:length(yforecast1_val)),'blue');
//xtitle('One-step-ahead - Validation: Continuous line => Real datas; Delineated line => Model monoobjective','SAMPLES','VOLTAGE');

// Six-step-ahead
scf(6)
plot(buck_val(:,1),buck_val(:,3),'red');
plot(buck_val(2:(length(buck_val)/3)-1,1),yforecast6_val(3:length(yforecast6_val)),'blue');
//xtitle('Six-step-ahead - Validation: Continuous line => Real datas; Delineated line => Model monoobjective','SAMPLES','VOLTAGE');

// Free prediction
scf(7)
plot(buck_val(:,1),buck_val(:,3),'red');
plot(buck_val(2:((length(buck_val)/3)-1),1),yforecastfree_val(3:length(yforecastfree_val)),'blue');
//xtitle('Free prediction - Validation: Continuous line => Real datas; Delineated line => Model monoobjective','SAMPLES','VOLTAGE');
