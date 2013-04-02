%function [H,W,De,Dv,f_rank] = build_incidence_matrix(file_friend,file_mentor, file_group, file_item,alpha)
function [final_tensor] = build_N_way_tensor(file_friend,file_mentor, file_item, alpha)
tt2 = now; % Current Time
P_Tot = [];

% Below we build the superset P_Tot of all vertices spanning over the 
% whole training region as well as the test period.
% For each of the 10 months
% This loop is same for both testing and training runs of this code.
for i = 1:10 
    
    % Mentor Data (mentor_id, mentee_id)
    file_mentor_new = strcat(file_mentor, int2str(i),'.csv');    M1 = load(file_mentor_new);
        
    % Friend Data (char_id, friend_id)
    file_friend_new = strcat(file_friend, int2str(i),'.csv');    M2 = load(file_friend_new);
        
    % Trade Data (buyer_id, seller_id)
    file_item_new = strcat(file_item, int2str(i),'.csv');    M3 = load(file_item_new);
        
    % Vertices set
    P_Tot = unique([P_Tot unique(M1)'  unique(M2)' unique(M3)']);
    
end

disp('Time taken For node calculation......');
disp(datevec(now - tt2));
size(P_Tot)

% N X N x L X T (L = 3 link types; T= 9 weeks)
N = size(P_Tot,2); L = 3; T = 9;
dim = [N N L T]; subt = [];

% This loop is the main loop which shall fill the tensor with appropriate values.
% For each week
% for T = 1:9 % for training of 9 months
for T = 10 % for testing of the 10th month
    tt1 = now; % Current Time
    
    % Mentor Data (mentor_id, mentee_id)
    file_mentor_new = strcat(file_mentor, int2str(T),'.csv'); M1 = load(file_mentor_new);
    
    % Friend Data (char_id, friend_id)
    file_friend_new = strcat(file_friend, int2str(T),'.csv'); M2 = load(file_friend_new);
    
    % Trade Data (buyer_id, seller_id)
    file_item_new = strcat(file_item, int2str(T),'.csv'); M3 = load(file_item_new);
    
    % Vertices Set
    %P_Tot = unique([P_Tot unique(M1)'  unique(M2)' unique(M3)']);
    
    i=1; L = 1; % Mentor Link
    while(true)
        if(M1(i,1)~=M1(i,2))
            n1 = find(P_Tot==M1(i,1),1); n2 = find(P_Tot==M1(i,2),1);
            subt = [subt; n1 n2 L T];
        end
        i=i+1;
        if(i>size(M1,1))
            break;
        end
    end
    
    i=1; L = 2; % Friend Link
    while(true)
        if(M2(i,1)~=M2(i,2))
            n1 = find(P_Tot==M2(i,1),1); n2 = find(P_Tot==M2(i,2),1);
            subt = [subt; n1 n2 L T];
        end
        i=i+1;
        if(i>size(M2,1))
            break;
        end
    end
    
    i=1; L = 3; % Trade Link
    while(true)
        if(M3(3,1)~=M3(i,2))
            n1 = find(P_Tot==M3(i,1),1); n2 = find(P_Tot==M3(i,2),1);
            subt = [subt; n1 n2 L T];
        end
        i=i+1;
        if(i>size(M3,1))
            break;
        end
    end
    
    disp('Time taken for data loading...... ');
    disp(datevec(now - tt1));    
end

final_tensor = subt;
%final_tensor = sptensor(subt,ones(size(subt,1),1),dim);
display('Tensor making done.');
disp('Time taken for complete tesor building ......');
disp(datevec(now - tt2));
end