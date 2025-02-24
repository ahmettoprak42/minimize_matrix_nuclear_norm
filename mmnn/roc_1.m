function [auc] = roc_1(deci,label_y,colour)
[threshold,ind] = sort(deci,'descend');
roc_y = label_y(ind);
stack_x = cumsum(roc_y == 0)/sum(roc_y == 0);
stack_y = cumsum(roc_y == 1)/sum(roc_y == 1);
auc=sum((stack_x(2:length(roc_y))-stack_x(1:length(roc_y)-1)).*stack_y(2:length(roc_y)));
plot(stack_x,stack_y,colour);
xlabel('False Positive Rate');
ylabel('True Positive Rate');
% hold on;
title(['ROC curve of (AUC = ' num2str(auc) ' )']);
grid on
end
