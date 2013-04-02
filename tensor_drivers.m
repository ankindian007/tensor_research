file_trust = 'housing_net_guk.csv';
file_trade = 'test_item_data_50000.txt';

EQ2_tensor_edgetype = build_tensor(file_trust,file_trade);
P = parafac_als(EQ2_tensor_edgetype,50);

for i = 1:50
    %[I,T] = sort(P.U{1}(:,i)); T(1:10)
    [I1,T1] = sort(P.U{2}(:,i)); T1(1:10)
    [I2,T2] = sort(P.U{3}(:,i)); T2(1:2)
end

file_friend = 'C:\Ankit\Data\GUK\Weekly\Housing Per month (charid,friendid)_';
file_trade = 'C:\Ankit\Data\GUK\Weekly\Trade Per month (buyer,seller)_';
file_mentor = 'C:\Ankit\Data\GUK\Weekly\Mentor Per month (mentor,mentee)_';
alpha = 0;

subt = build_N_way_tensor(file_friend,file_mentor, file_trade,alpha);
output_file = ['C:\Ankit\Data\GUK\Weekly\EQ2_N_way_tensor_1.mat'];
save(output_file,'EQ2_N_way_tensor');
dim = [21540 21540 3 9];
final_tensor = sptensor(subt,ones(size(subt,1),1),dim); 
EQ2_N_way_tensor = final_tensor;
% Non-logit code commendted below
%{
P = cp_als(EQ2_N_way_tensor,50);
output_file = ['C:\Ankit\Data\GUK\Weekly\CP_ALS_output_factors=50.mat'];
save(output_file,'P');
P = cp_als(EQ2_N_way_tensor,10);
output_file = ['C:\Ankit\Data\GUK\Weekly\CP_ALS_output_factors=10.mat'];
save(output_file,'P');
%}
% Logit performed case
Z = EQ2_N_way_tensor;
T_log = elemfun(Z, @(x) log(x)+1);
P_50 = cp_als(T_log ,50);
output_file = ['C:\Ankit\Data\GUK\Weekly\CP_ALS_LOG_output_factors=50.mat'];
save(output_file,'P');
P_10 = cp_als(T_log ,10);
output_file = ['C:\Ankit\Data\GUK\Weekly\CP_ALS_LOG_output_factors=10.mat'];
save(output_file,'P');

% Testing - no need of factorization and 
file_friend = 'C:\Ankit\Data\GUK\Weekly\Housing Per month (charid,friendid)_';
file_trade = 'C:\Ankit\Data\GUK\Weekly\Trade Per month (buyer,seller)_';
file_mentor = 'C:\Ankit\Data\GUK\Weekly\Mentor Per month (mentor,mentee)_';
alpha = 0;
subt = build_N_way_tensor(file_friend,file_mentor, file_trade,alpha);
%subt = EQ2_N_way_tensor_test;
dim = [21540 21540 3 1];
final_tensor_test = sptensor(subt,ones(size(subt,1),1),dim); 
EQ2_N_way_tensor_test = final_tensor_test;
output_file = ['C:\Ankit\Data\GUK\Weekly\EQ2_N_way_tensor_TEST.mat'];
save(output_file,'EQ2_N_way_tensor_test');
Z_test = EQ2_N_way_tensor_test;
T_log_test = elemfun(Z_test, @(x) log(x)+1);

%{
P_test = cp_als(EQ2_N_way_tensor,50);
output_file = ['C:\Ankit\Data\GUK\Weekly\CP_ALS_output_factors=50_TEST.mat'];
save(output_file,'P_test');
P_test = cp_als(EQ2_N_way_tensor,10);
output_file = ['C:\Ankit\Data\GUK\Weekly\CP_ALS_output_factors=10_TEST.mat'];
save(output_file,'P_test');
%}

for t= 1:9
    f = P.U{4}(:,t); f(f<0)=0;
    subplot(9,1,t);plot([1:9]',f)
    title(num2str(t));
end

for t= 1:3
    f = P.U{3}(:,t); f(f<0)=0;
    subplot(3,1,t);plot([1:3]',f)
    title(num2str(t));
end

%{
S_k_n=[];
top = [];
r = 1000;
for n = 1:21540
    for k = 1:50
        gamma_k =  sum(P.U{4}(1:3,k))/3; % time dim as average of the last 3 time steps
        lambda_k = P.lambda(k);          % lambda of the factors
        alpha_k = P.U{3}(1,k);           % link dim for link type = 1
        S_k_n = S_k_n + lambda_k*gamma_k*alpha_k*P.U{1}(n,k)*P.U{2}(:,k);
    end
    idx = sort(S_k_n,1); val = sort(S_k_n,2); 
    S_k_n_top = [n*ones(r,1); idx(1:r); val(1:r)];
    S_k_n = [top; S_k_n_top];
    temp = sortrows(S_k_n,3);    
    top = temp(1:1000,:);
    S_k_n = []; S_k_n_top = [];
end
%}

P.U{3}(:,1)*P.U{2}(:,1)'

plot([1:9]',P.U{4}(:,1));
hold on;
plot([1:9]',P.U{4}(:,2));
hold on;
plot([1:9]',P.U{4}(:,3));
hold on;
plot([1:9]',P.U{4}(:,4));
hold on;
plot([1:9]',P.U{4}(:,5));
hold on;
plot([1:9]',P.U{4}(:,6));
hold on;
plot([1:9]',P.U{4}(:,7));
hold on;
plot([1:9]',P.U{4}(:,8));
hold on;
plot([1:9]',P.U{4}(:,9));
hold off;

results = [];

% Building result for # of factors = 10

r = 100; fac = 10; l_type = 1;
[pr1, pr2] = result_evaluation(1,10,100,P_10,T_log,T_log_test); 
results = [results; r fac l_type pr1 pr2];
r = 1000; fac = 10; l_type = 1;
[pr1, pr2] = result_evaluation(l_type,fac,r,P_10,T_log,T_log_test); 
results = [results; r fac l_type pr1 pr2];
r = 10000; fac = 10; l_type = 1;
[pr1, pr2] = result_evaluation(l_type,fac,r,P_10,T_log,T_log_test); 
results = [results; r fac l_type pr1 pr2];
r = 100; fac = 10; l_type = 2;
[pr1, pr2] = result_evaluation(l_type,fac,r,P_10,T_log,T_log_test); 
results = [results; r fac l_type pr1 pr2];
r = 1000; fac = 10; l_type = 2;
[pr1, pr2] = result_evaluation(l_type,fac,r,P_10,T_log,T_log_test); 
results = [results; r fac l_type pr1 pr2];
r = 10000; fac = 10; l_type = 2;
[pr1, pr2] = result_evaluation(l_type,fac,r,P_10,T_log,T_log_test); 
results = [results; r fac l_type pr1 pr2];
r = 100; fac = 10; l_type = 3;
[pr1, pr2] = result_evaluation(l_type,fac,r,P_10,T_log,T_log_test); 
results = [results; r fac l_type pr1 pr2];
r = 1000; fac = 10; l_type = 3;
[pr1, pr2] = result_evaluation(l_type,fac,r,P_10,T_log,T_log_test); 
results = [results; r fac l_type pr1 pr2];
r = 10000; fac = 10; l_type = 3;
[pr1, pr2] = result_evaluation(l_type,fac,r,P_10,T_log,T_log_test); 
results = [results; r fac l_type pr1 pr2];

r = 21540; fac = 10; l_type = 1;
[pr1, pr2, tot_all, tot_new] = result_evaluation(l_type,fac,r,P_10,T_log,T_log_test); 
results = [results; r fac l_type pr1 pr2 tot_all tot_new];
r = 21540; fac = 10; l_type = 2;
[pr1, pr2, tot_all, tot_new] = result_evaluation(l_type,fac,r,P_10,T_log,T_log_test); 
results = [results; r fac l_type pr1 pr2 tot_all tot_new];
r = 21540; fac = 10; l_type = 3;
[pr1, pr2, tot_all, tot_new] = result_evaluation(l_type,fac,r,P_10,T_log,T_log_test); 
results = [results; r fac l_type pr1 pr2 tot_all tot_new];

% Building result for # of factors = 50

r = 100; fac = 50; l_type = 1;
[pr1, pr2] = result_evaluation(1,10,100,P_50,T_log,T_log_test); 
results = [results; r fac l_type pr1 pr2];
r = 1000; fac = 50; l_type = 1;
[pr1, pr2] = result_evaluation(l_type,fac,r,P_50,T_log,T_log_test); 
results = [results; r fac l_type pr1 pr2];
r = 10000; fac = 50; l_type = 1;
[pr1, pr2] = result_evaluation(l_type,fac,r,P_50,T_log,T_log_test); 
results = [results; r fac l_type pr1 pr2];
r = 100; fac = 50; l_type = 2;
[pr1, pr2] = result_evaluation(l_type,fac,r,P_50,T_log,T_log_test); 
results = [results; r fac l_type pr1 pr2];
r = 1000; fac = 50; l_type = 2;
[pr1, pr2] = result_evaluation(l_type,fac,r,P_50,T_log,T_log_test); 
results = [results; r fac l_type pr1 pr2];
r = 10000; fac = 50; l_type = 2;
[pr1, pr2] = result_evaluation(l_type,fac,r,P_50,T_log,T_log_test); 
results = [results; r fac l_type pr1 pr2];
r = 100; fac = 50; l_type = 3;
[pr1, pr2] = result_evaluation(l_type,fac,r,P_50,T_log,T_log_test); 
results = [results; r fac l_type pr1 pr2];
r = 1000; fac = 50; l_type = 3;
[pr1, pr2] = result_evaluation(l_type,fac,r,P_50,T_log,T_log_test); 
results = [results; r fac l_type pr1 pr2];
r = 10000; fac = 50; l_type = 3;
[pr1, pr2] = result_evaluation(l_type,fac,r,P_50,T_log,T_log_test); 
results = [results; r fac l_type pr1 pr2];

