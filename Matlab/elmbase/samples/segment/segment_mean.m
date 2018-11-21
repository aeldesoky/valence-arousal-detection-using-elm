number_of_trial=50;
test=zeros(number_of_trial,1);
train=zeros(number_of_trial,1);
train_time=zeros(number_of_trial,1);
testing_time=zeros(number_of_trial,1);

wb=waitbar(0,'Please waiting...');

rnd=0;
waitbar(rnd/number_of_trial,wb);

for rnd = 1 : number_of_trial
    
    disp('Randomly generating training and testing dataset ... ...');
    segment_data;

    disp('Start training ... ...');
    [learn_time, test_time, train_accuracy, test_accuracy]=ELM('segment_train','segment_test',1,200,'sig');
    test(rnd,1)=test_accuracy;
    train(rnd,1)=train_accuracy;
    train_time(rnd,1)=learn_time;
    testing_time(rnd,1)=test_time;

    waitbar(rnd/number_of_trial,wb);
    
end
close(wb);

AverageTrainingTime=mean(train_time)
StandardDeviationofTrainingTime=std(train_time)
AvergeTestingTime=mean(testing_time)
StandardDeviationofTestingTime=std(testing_time)
AverageTrainingAccuracy=mean(train)
StandardDeviationofTrainingAccuracy=std(train)
AverageTestingAccuracy=mean(test)
StandardDeviationofTestingAccuracy=std(test)
