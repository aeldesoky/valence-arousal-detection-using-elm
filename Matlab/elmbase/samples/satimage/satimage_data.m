function [P,T,TV] = GetPTforSatimagedata

    load sat_trn.dt;
    load sat_tst.dt;
    
    dataset(1:size(sat_trn,1),:)=sat_trn;
    dataset(size(sat_trn,1)+1:size(sat_trn,1)+size(sat_tst,1),:)=sat_tst;

    rand_sequence=randperm(size(dataset,1));
    temp_dataset=dataset;
    
    dataset=temp_dataset(rand_sequence, :);

    for i=1:size(dataset,2)-1
        dataset(:,i)=(dataset(:,i)-min(dataset(:,i)))/(max(dataset(:,i))-min(dataset(:,i)))*2-1;
    end             
    
    P1=dataset(1:floor(size(dataset,1)/2),1:size(dataset,2)-1);
    T1=dataset(1:floor(size(dataset,1)/2),size(dataset,2));

%    P=P1';
%    T=T1';
    
    %Obtain Random Validation Matrix
    
    X=dataset(floor(size(dataset,1)/2)+1:size(dataset,1),1:size(dataset,2)-1);
    Y=dataset(floor(size(dataset,1)/2)+1:size(dataset,1),size(dataset,2));
    
    
%    TV.P=X';
%    TV.T=Y';
    
    fid = fopen('sat_train','w');
    for i=1:size(P1,1)
        fprintf(fid,'%2.8f ',T1(i,1));
        for j=1:size(P1,2)
%            fprintf(fid,' %d:%2.8f',j, P1(i,j));    %   for SVM
            fprintf(fid,' %2.8f', P1(i,j));    %   for ELM
        end
            fprintf(fid,'\n');
        end
    fclose(fid);

    fid = fopen('sat_test','w');    
    for i=1:size(X,1)
        fprintf(fid,'%2.8f ',Y(i,1));
        for j=1:size(X,2)
%            fprintf(fid,' %d:%2.8f',j, X(i,j));     %   for SVM
            fprintf(fid,' %2.8f', X(i,j));     %   for ELM
        end
            fprintf(fid,'\n');
        end
    fclose(fid);