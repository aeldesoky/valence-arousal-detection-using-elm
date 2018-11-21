function [P, T, TV, TP] = GetPTforSinc
    %Obtain Random P, T
    X1=20*rand(1,5000)-10;
    
    for i=1:size(X1,2)
        if X1(1,i)==0
            X1(1, i)=0.000001;
        end
    end    
    Y1 = [sin(X1).*X1.^(-1); X1];

    for i=1:size(X1,2)
        Y1(1, i)=Y1(1, i)+0.4*rand(1,1)-0.2;
    end        
    
    fid = fopen('sinc_train','w');
    fprintf(fid,'%2.8f %2.8f\n',Y1);
    fclose(fid);

    X2=sort(20*rand(1,5000)-10);
    
    for i=1:size(X2,2)
        if X2(1,i)==0
            X2(1, i)=0.000001;
        end
    end    
    Y2 = [sin(X2).*X2.^(-1); X2];
    
    fid = fopen('sinc_test','w');
    fprintf(fid,'%2.8f %2.8f\n',Y2);
    fclose(fid);
    