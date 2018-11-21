function [P,T,TV] = Segment_Data
    load segment1.dat
    
    dataset=segment1;
    
    rand_sequence=randperm(size(dataset,1));
    temp_dataset=dataset;
    
    dataset=temp_dataset(rand_sequence, :);

    for i=1:size(dataset,2)-1
        if max(dataset(:,i))~=min(dataset(:,i))
            dataset(:,i)=(dataset(:,i)-min(dataset(:,i)))/(max(dataset(:,i))-min(dataset(:,i)))*2-1;
        else
            dataset(:,i)=0;
        end
    end             
    
    P1=dataset(1:1500,1:size(dataset,2)-1);
    T1=dataset(1:1500,size(dataset,2));

%    P=P1';
%    T=T1';
    
    %Obtain Random Validation Matrix
    
    X=dataset(1501:size(dataset,1),1:size(dataset,2)-1);
    Y=dataset(1501:size(dataset,1),size(dataset,2));
    
    
%    TV.P=X';
%    TV.T=Y';
    
    fid = fopen('segment_train','w');
    for i=1:size(P1,1)
        fprintf(fid,'%2.8f ',T1(i,1));
        for j=1:size(P1,2)
%            fprintf(fid,' %d:%2.8f',j, P1(i,j));    %   for SVM
            fprintf(fid,' %2.8f', P1(i,j));    %   for ELM
        end
            fprintf(fid,'\n');
        end
    fclose(fid);

    fid = fopen('segment_test','w');    
    for i=1:size(X,1)
        fprintf(fid,'%2.8f ',Y(i,1));
        for j=1:size(X,2)
%            fprintf(fid,' %d:%2.8f',j, X(i,j));     %   for SVM
            fprintf(fid,' %2.8f', X(i,j));     %   for ELM
        end
            fprintf(fid,'\n');
        end
    fclose(fid);
    