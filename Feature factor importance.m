
X=[TrainSum(:,1:4),JiangShui,TrainSum(:,13),TrainSum(:,3)./TrainSum(:,2),TrainSum(:,14)];


train_data = X(:,1:end-1);train_label = X(:,end);




%%
b = TreeBagger(100,train_data, train_label,'Method','regression','Surrogate','on','OOBPrediction','On','OOBPredictorImportance','On',...
   'PredictorSelection','curvature');
%%
% Feature factor importance
idxvar = find(b.OOBPermutedPredictorDeltaError>0.7)
idxCategorical = find(iscategorical(idxvar)==1);
finbag = zeros(1,b.NTrees);figure
bar(b.OOBPermutedPredictorDeltaError)
xlabel('Feature') 
ylabel('Out-of-Bag Feature Importance')




