function [final_tensor] = build_tensor(file_trust, file_trade)

    trust_data = load(file_trust);
    [acc_1, acc_2, ~ , ~] = textread(file_trade,'%d %d %s %s','delimiter',',');    
    
    trade_data(:,1) = acc_1; trade_data(:,2) = acc_2;
    
    nodes  = unique([ trust_data(:,1)' trust_data(:,3)' trade_data(:,1)' trade_data(:,2)']);
    %nodes  = unique([ trust_data(:,1)' trust_data(:,3)]);
    [dim ,~] = size(nodes);
    
    subs = [];    
    
    %trust_network = zeros(dim,dim);
    [n ,~] = size(trust_data)            
    
    for i = 1 : n
        %trust_network( find(nodes==trust_data(i,1)) , find(nodes==trust_data(i,3)) )=1;
        %trust_network( find(nodes==trust_data(i,3)) , find(nodes==trust_data(i,1)) )=1;
        
        t1 = find(nodes==trust_data(i,1)); t2 = find(nodes==trust_data(i,3));
        
        if(size(t1,1)~=0 && size(t2,1)~=0)
            
            if(t1(1) ~= t2(1))
                %trust_network( t1(1),t2(1) )=1;
                %trust_network( t2(1),t1(1) )=1;
                subs = [subs; t1(1) t2(1) 1; t2(1) t1(1) 1];               
            end
        end        
    end
    
    %{
    [count, ~]=size(find(trust_network==1));
    disp('Number of Edges ......');
    disp(count/2);    
    disp('Matrix Size ......');
    disp(dim^2);
    
    disp('Saving the file ......');
    output_file = ['C:\Users\ankindian\Desktop\Work\Mohammed\Bi Graph IR work\code\Outputs\' file_trust '_trust_network.mat'];
    save(output_file,'trust_network');
    %}
      
    %trade_network = zeros(dim,dim);
    [n ,~] = size(trade_data)        
    
    for i = 1 : n
        
        t1 = find(nodes==trade_data(i,1)); t2 = find(nodes==trade_data(i,2));
        
        if(size(t1,1)~=0 && size(t2,1)~=0)
            
            if(t1(1) ~= t2(1))
                %trade_network( t1(1),t2(1) )=1;
                %trade_network( t2(1),t1(1) )=1;
                subs = [subs; t1(1) t2(1) 1; t2(1) t1(1) 2];               
            end
        end
    end
    
    %{
    [count, ~]=size(find(trade_network==1));
    disp('Number of Edges ......');
    disp(count/2);    
    disp('Matrix Size ......');
    disp(dim^2);
    
    disp('Saving the file ......');
    output_file = ['C:\Users\ankindian\Desktop\Work\Mohammed\Bi Graph IR work\code\Outputs\Trade\' file_trade '_trade_network.mat'];
    save(output_file,'trade_network');
    %}
    
    dim
    size(nodes,2)
    size(subs,1)
    final_tensor = sptensor(subs,ones(size(subs,1),1),[size(nodes,2) size(nodes,2) 2]);
    
end