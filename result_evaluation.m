function [pred_all_link, pred_new_link, tot_all, tot_new] = result_evaluation(Link_Type,No_of_factors,r,P,T_log,T_log_test)
% Link_Type = 1; No_of_factors = 50; r = 100000; 
N_nodes = 21540; % Size of P_Tot
S_k_n=zeros(N_nodes,1);
top = [];
tic; % r is the top what # of results we want 
for n = 1:N_nodes
    
    for k = 1:No_of_factors
        gamma_k =  sum(P.U{4}(1:3,k))/3; % time dim as the average of last 3 steps
        lambda_k = P.lambda(k);          % lambda
        alpha_k = P.U{3}(Link_Type,k);   % link dim for link type = Link_Type
        temp1 = lambda_k*gamma_k*alpha_k*P.U{1}(n,k)*P.U{2}(:,k) ;
        S_k_n = S_k_n + temp1;
    end
    
    %idx = sort(S_k_n,1); val = sort(S_k_n,2); 
    [val,idx] = sort(S_k_n,'descend'); % sorted array
    S_k_n_top = [n*ones(r,1), idx(1:r), val(1:r)]; % top r from the sorted array
    S_k_n = [top; S_k_n_top]; % combining with the current top r array
    temp2 = flipdim(sortrows(S_k_n,3),1);    
    top = temp2(1:r,:);    % Top r 
    S_k_n = zeros(N_nodes,1); S_k_n_top = [];    
end
toc

% Get all the Link_Type=1 type records for all time snaps
train_log = T_log(:,:,Link_Type,:); 
display('train_log'); size(train_log)
% Get max on the mode-time(3) as train_log is 3-D
% Below is the n x n matrix for link = 1 that shall be having 
% all the links formed till T=9th snapshot.
train_log_collapse = collapse(train_log,3,@max); 
display('train_log_collapse.subs'); size(train_log_collapse.subs)
% Find the complete T=10 snapshot tensor 
% by combining test_tensor with the train_log_collapse calculated 
% above which shall be used for testing
test_log = T_log_test(:,:,Link_Type,:); % this subt are the subs that are produced by the 
display('test_log'); size(test_log)   % build_tensor function for T=10 test time snapshot
test_log_collapse = collapse(test_log,3,@max); 
display('test_log_collapse.subs'); size(test_log_collapse.subs)                    

complete_test = union(test_log_collapse.subs,train_log_collapse.subs,'rows');
display('complete_test'); tot_all = size(complete_test)
% # of new links formed
new_links = setdiff(complete_test,train_log_collapse.subs,'rows');
display('new_links'); tot_new = size(new_links,1)

%sum(test_tensor(:,1:2) == top(:,1:2),2)
% Correctly predicted link both existing and non existing
top_tuples = top(:,1:2);
TP = 0; 
for g = 1:size(complete_test,1)
    if(sum(ismember(top_tuples,complete_test(g,1:2),'rows'))>0)
       TP = TP + 1; 
    end
end
fprintf(['Correctly predicted link both existing and non existing for r = %d link_type= %d no_of_fact = %d'],r, Link_Type, No_of_factors);
TP
pred_all_link = TP;
% Correctly predicted links non existing
top_tuples = top(:,1:2);
TP = 0; 
for g = 1:size(new_links,1)
    if(sum(ismember(top_tuples,new_links(g,1:2),'rows'))>0)
       TP = TP + 1; 
    end
end
fprintf(['Correctly predicted links non existing for r = %d link_type= %d no_of_fact = %d'],r, Link_Type, No_of_factors);
TP
pred_new_link = TP;
end