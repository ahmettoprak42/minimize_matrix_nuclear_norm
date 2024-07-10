function W=Label_Propagation(feature_matrix,tag,neighbor_num,regulation)
    distance_matrix=calculate_instances(feature_matrix);
    nearst_neighbor_matrix=calculate_neighbors(distance_matrix,neighbor_num);
    W=optimization_similairty_matrix(feature_matrix,nearst_neighbor_matrix,tag,regulation);
end

function W=optimization_similairty_matrix(feature_matrix,nearst_neighbor_matrix,tag,regulation)
   row_num=size(feature_matrix,1);
   W=zeros(1,row_num);
   if tag==1
       row_num=1;
   end
   for i=1:row_num
       nearst_neighbors=feature_matrix(logical(nearst_neighbor_matrix(i,:)'),:);   
       neighbors_num=size(nearst_neighbors,1);
       G1=repmat(feature_matrix(i,:),neighbors_num,1)-nearst_neighbors;
       G2=repmat(feature_matrix(i,:),neighbors_num,1)'-nearst_neighbors';
       if regulation=='regulation2'
         G_i=G1*G2+eye(neighbors_num);
       end
       if regulation=='regulation1'
         G_i=G1*G2;
       end
       H=2*G_i;
       f=[];
       A=[];
       if isempty(H)
           A;
       end
       
       b=[];
       Aeq=ones(neighbors_num,1)';
       beq=1;
       lb=zeros(neighbors_num,1);
       ub=[];
       options=optimset('Display','off');
       [w,fval]= quadprog(H,f,A,b,Aeq,beq,lb,ub,[],options);
       w=w';
       W(i,logical(nearst_neighbor_matrix(i,:)))=w;     
   end
end

function distance_matrix=calculate_instances(feature_matrix)
    [row_num,col_num]=size(feature_matrix);
    distance_matrix=zeros(row_num,row_num);
    for i=1:row_num
        for j=i+1:row_num
            distance_matrix(i,j)=sqrt(sum((feature_matrix(i,:)-feature_matrix(j,:)).^2));
            distance_matrix(j,i)=distance_matrix(i,j);
        end
        distance_matrix(i,i)=col_num;
    end
end

function nearst_neighbor_matrix=calculate_neighbors(distance_matrix,neighbor_num)
  [sv si]=sort(distance_matrix,2,'ascend');
  [row_num,col_num]=size(distance_matrix);
  nearst_neighbor_matrix=zeros(row_num,col_num);
  index=si(:,1:neighbor_num);
  for i=1:row_num
       nearst_neighbor_matrix(i,index(i,:))=1;
  end
end

