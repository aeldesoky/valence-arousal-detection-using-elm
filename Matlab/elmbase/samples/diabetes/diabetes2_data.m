function [P, T, TV] = diabetes2_data
    %Obtain Random P, T
    
    load diabetes2.dt;
    
    rand_sequence=randperm(size(diabetes2,1));
    temp_dataset=diabetes2;
    
    diabetes2=temp_dataset(rand_sequence, :);
    
    P1=diabetes2(1:576,1:8);
    T1=diabetes2(1:576,9);

%    P=P1'.*2-1;
%    T=T1'.*2-1;
    
    %Obtain Random Test Matrix
    X=diabetes2(577:size(diabetes2,1),1:8);
    Y=diabetes2(577:size(diabetes2,1),9);
%    TV.P=X'.*2-1;
%    TV.T=Y'.*2-1;

    fid = fopen('diabetes_train','w');
    for i=1:size(P1,1)
        fprintf(fid,'%2.8f ',T1(i,1));
        for j=1:size(P1,2)
%            fprintf(fid,' %d:%2.8f',j, P1(i,j));    % for SVM
            fprintf(fid,' %2.8f', P1(i,j));    % for ELM
        end
            fprintf(fid,'\n');
        end
    fclose(fid);

    fid = fopen('diabetes_test','w');    
    for i=1:size(X,1)
        fprintf(fid,'%2.8f ',Y(i,1));
        for j=1:size(X,2)
%            fprintf(fid,' %d:%2.8f',j, X(i,j));     % for SVM
            fprintf(fid,' %2.8f', X(i,j));     % for ELM
        end
            fprintf(fid,'\n');
        end
    fclose(fid);
